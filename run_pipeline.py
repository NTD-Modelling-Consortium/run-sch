#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import enum
from pathlib import Path
from typing import Dict, List, Set, Optional
from dataclasses import dataclass


class Stage(enum.Enum):
    FITTING_PREP = "fitting-prep"
    FITTING = "fitting"
    PROJECTIONS_PREP = "projections-prep"
    NEARTERM_PROJECTIONS = "nearterm-projections"
    ALL = "all"
    SKIP_FITTING_PREP = "skip-fitting-prep"


class Disease(enum.Enum):
    STH = "sth"
    SCH = "sch"
    ALL = "all"


class STHSpecies(enum.Enum):
    ASCARIS = "ascaris"
    HOOKWORM = "hookworm"
    TRICHURIS = "trichuris"
    ALL = "all"


class SCHSpecies(enum.Enum):
    HAEMATOBIUM = "haematobium"
    MANSONI_HIGH = "mansoni_high_burden"
    MANSONI_LOW = "mansoni_low_burden"
    ALL = "all"


@dataclass
class StageArtifact:
    """Represents an artifact required by a stage."""

    path_pattern: str  # Path pattern with {species} and {batch_id} placeholders
    required_for_species: Set[str]  # Which species require this artifact
    is_batch_specific: bool = False  # Whether artifact depends on batch_id
    description: str = ""


@dataclass
class StageDependency:
    """Represents dependencies for a pipeline stage."""

    stage: Stage
    artifacts: List[StageArtifact]
    upstream_stages: Set[Stage]  # Which stages must complete before this one


def _get_base_path() -> Path:
    return Path(os.getenv("STH_SCH_AMIS_DIR", "/ntdmc/sth-sch-amis-integration"))


def _get_species_list(species_enums: List) -> List[str]:
    """Convert species enum list to string list."""
    species_strings = []

    for species_enum in species_enums:
        if species_enum == STHSpecies.ALL:
            species_strings.extend(["ascaris", "hookworm", "trichuris"])
        elif species_enum == SCHSpecies.ALL:
            species_strings.extend(
                ["haematobium", "mansoni_low_burden", "mansoni_high_burden"]
            )
        elif isinstance(species_enum, (STHSpecies, SCHSpecies)):
            species_strings.append(species_enum.value)

    return list(set(species_strings))  # Remove duplicates


def parse_species_arguments(args, disease):
    """Parse and validate species arguments for both STH and SCH diseases."""
    sth_species = []
    sch_species = []

    def parse_species_string(species_string):
        """Parse comma-separated or space-separated species string."""
        if not species_string:
            return []
        
        if "," in species_string:
            return [s.strip() for s in species_string.split(",")]
        else:
            return species_string.split()

    def validate_species_choices(species_list, valid_species, disease_name):
        """Validate species choices against valid options."""
        for species in species_list:
            if species not in valid_species:
                print(
                    f"Error: Invalid {disease_name} species '{species}'. Valid choices: {valid_species}",
                    file=sys.stderr,
                )
                sys.exit(1)

    # Parse STH species
    if disease in [Disease.STH, Disease.ALL]:
        if args.sth_species:
            species_list = parse_species_string(args.sth_species)
            valid_sth_species = [s.value for s in STHSpecies]
            validate_species_choices(species_list, valid_sth_species, "STH")
            sth_species = [STHSpecies(s) for s in species_list]
        else:
            sth_species = [STHSpecies.ALL]

    # Parse SCH species  
    if disease in [Disease.SCH, Disease.ALL]:
        if args.sch_species:
            species_list = parse_species_string(args.sch_species)
            valid_sch_species = [s.value for s in SCHSpecies]
            validate_species_choices(species_list, valid_sch_species, "SCH")
            sch_species = [SCHSpecies(s) for s in species_list]
        else:
            sch_species = [SCHSpecies.ALL]

    return sth_species, sch_species


def _get_pipeline_dependencies() -> Dict[Stage, StageDependency]:
    """Define the complete pipeline dependency graph."""
    return {
        Stage.FITTING_PREP: StageDependency(
            stage=Stage.FITTING_PREP,
            artifacts=[],  # No dependencies - generates base artifacts
            upstream_stages=set(),
        ),
        Stage.FITTING: StageDependency(
            stage=Stage.FITTING,
            artifacts=[
                # STH fitting-prep artifacts
                StageArtifact(
                    path_pattern="fitting-prep/artefacts/Maps/iu_task_lookup_sth.rds",
                    required_for_species={"ascaris", "hookworm"},
                    description="STH task lookup table",
                ),
                StageArtifact(
                    path_pattern="fitting-prep/artefacts/Maps/iu_task_lookup_trichuris.rds",
                    required_for_species={"trichuris"},
                    description="Trichuris task lookup table",
                ),
                StageArtifact(
                    path_pattern="fitting-prep/artefacts/endgame_inputs/STH/InputMDA_MTP_{batch_id}.xlsx",
                    required_for_species={"ascaris", "hookworm", "trichuris"},
                    is_batch_specific=True,
                    description="STH coverage data for batch",
                ),
                # SCH fitting-prep artifacts
                StageArtifact(
                    path_pattern="fitting-prep/artefacts/Maps/iu_task_lookup_haema.rds",
                    required_for_species={"haematobium"},
                    description="Haematobium task lookup table",
                ),
                StageArtifact(
                    path_pattern="fitting-prep/artefacts/Maps/iu_task_lookup_{species}.rds",
                    required_for_species={"mansoni_low_burden", "mansoni_high_burden"},
                    description="Mansoni species-specific lookup table",
                ),
                StageArtifact(
                    path_pattern="fitting-prep/artefacts/endgame_inputs/sch-haematobium/InputMDA_MTP_{batch_id}.xlsx",
                    required_for_species={"haematobium"},
                    is_batch_specific=True,
                    description="Haematobium coverage data for batch",
                ),
                StageArtifact(
                    path_pattern="fitting-prep/artefacts/endgame_inputs/sch-mansoni/InputMDA_MTP_{batch_id}.xlsx",
                    required_for_species={"mansoni_low_burden", "mansoni_high_burden"},
                    is_batch_specific=True,
                    description="Mansoni coverage data for batch",
                ),
            ],
            upstream_stages={Stage.FITTING_PREP},
        ),
        Stage.PROJECTIONS_PREP: StageDependency(
            stage=Stage.PROJECTIONS_PREP,
            artifacts=[
                # Fitting outputs - actual file location
                StageArtifact(
                    path_pattern="fitting/artefacts/AMIS_output/{species}_amis_output{batch_id}{sigma_suffix}.Rdata",
                    required_for_species={
                        "ascaris",
                        "hookworm",
                        "trichuris",
                        "haematobium",
                        "mansoni_low_burden",
                        "mansoni_high_burden",
                    },
                    is_batch_specific=True,
                    description="AMIS fitting results",
                ),
            ],
            upstream_stages={Stage.FITTING},
        ),
        Stage.NEARTERM_PROJECTIONS: StageDependency(
            stage=Stage.NEARTERM_PROJECTIONS,
            artifacts=[
                # Note: Nearterm-projections uses IU-specific files, not batch-specific
                # The validation for this stage should check for the existence of lookup tables
                # and let the projection script handle IU-to-file mapping
                StageArtifact(
                    path_pattern="fitting-prep/artefacts/Maps/table_iu_idx_STH.csv",
                    required_for_species={"ascaris", "hookworm"},
                    is_batch_specific=False,
                    description="STH IU lookup table for projections",
                ),
                StageArtifact(
                    path_pattern="fitting-prep/artefacts/Maps/table_iu_idx_trichuris.csv",
                    required_for_species={"trichuris"},
                    is_batch_specific=False,
                    description="Trichuris IU lookup table for projections",
                ),
                StageArtifact(
                    path_pattern="fitting-prep/artefacts/Maps/table_iu_idx_haematobium.csv",
                    required_for_species={"haematobium"},
                    is_batch_specific=False,
                    description="Haematobium IU lookup table for projections",
                ),
            ],
            upstream_stages={
                Stage.PROJECTIONS_PREP,
                Stage.FITTING_PREP,
            },  # Cross-stage dependency!
        ),
    }


def validate_stage_dependencies(
    stage: Stage,
    species_list: List,
    batch_id: Optional[int] = None,
    amis_sigma: Optional[float] = None,
) -> bool:
    """Validate that all dependencies for a given stage are available."""
    base_path = _get_base_path()
    dependencies = _get_pipeline_dependencies()

    if stage not in dependencies:
        print(f"   ⚠️  Unknown stage: {stage.value}")
        return False

    print(f"\n🔍 Validating dependencies for stage: {stage.value}")

    stage_dep = dependencies[stage]
    if not stage_dep.artifacts:
        print("   ✅ No dependencies required")
        return True

    # Convert species enums to strings
    requested_species = _get_species_list(species_list)
    missing_artifacts = []

    # Check each required artifact
    for artifact in stage_dep.artifacts:
        # Check if this artifact is needed for any of the requested species
        needed_species = artifact.required_for_species.intersection(
            set(requested_species)
        )
        if not needed_species:
            continue  # Skip artifacts not needed for requested species

        # Check if batch_id is required but not provided
        if artifact.is_batch_specific and batch_id is None:
            missing_artifacts.append(f"{artifact.description}: requires --id parameter")
            continue

        # Check artifact existence for each needed species
        for species in needed_species:
            # Build file path from pattern
            path_str = artifact.path_pattern

            # Handle sigma suffix for AMIS output files
            if "{sigma_suffix}" in path_str:
                # Check for files with any sigma value
                # First try default (no suffix)
                default_path = path_str.replace("{sigma_suffix}", "")
                # Also check with specific sigma if provided
                sigma_path = None
                if amis_sigma and amis_sigma != 0.0025:
                    sigma_path = path_str.replace(
                        "{sigma_suffix}", f"_sigma{amis_sigma}"
                    )

                # Check both possible locations
                found = False
                for check_path in [default_path, sigma_path]:
                    if check_path:
                        full_path = check_path.format(
                            species=species, batch_id=batch_id or ""
                        )
                        artifact_path = base_path / full_path
                        if artifact_path.exists():
                            found = True
                            break

                # If not found with specific sigma values, check for any sigma file
                if not found:
                    # Try pattern matching for any sigma value
                    pattern_path = path_str.replace("{sigma_suffix}", "_sigma*")
                    pattern_full = pattern_path.format(
                        species=species, batch_id=batch_id or ""
                    )
                    matching_files = list(base_path.glob(pattern_full))
                    if matching_files:
                        found = True
                        print(
                            f"   ℹ️  Found {artifact.description} for {species} with different sigma value"
                        )

                if not found:
                    missing_artifacts.append(f"{artifact.description} for {species}")
                continue

            # Handle other placeholders
            if "{species}" in path_str:
                path_str = path_str.format(species=species, batch_id=batch_id or "")
            elif "{batch_id}" in path_str:
                path_str = path_str.format(batch_id=batch_id)

            artifact_path = base_path / path_str

            if not artifact_path.exists():
                missing_artifacts.append(
                    f"{artifact.description} for {species}: {artifact_path}"
                )

    if missing_artifacts:
        print("   ❌ Missing dependencies:")
        for artifact in missing_artifacts:
            print(f"      - {artifact}")

        # Suggest upstream stages to run
        upstream_stages = stage_dep.upstream_stages
        if upstream_stages:
            stage_names = [s.value for s in upstream_stages]
            if batch_id:
                print(
                    f"   💡 Try running: --stage={' or --stage='.join(stage_names)} --id={batch_id}"
                )
            else:
                print(f"   💡 Try running: --stage={' or --stage='.join(stage_names)}")

        return False

    print("   ✅ All dependencies satisfied")
    return True


def validate_environment():
    """Validate required environment variables are set."""
    base_path = _get_base_path()

    if not base_path.exists():
        raise ValueError(f"Base path '{base_path}' does not exist")

    # Validate required paths exist
    required_paths = [
        base_path / "fitting-prep",
        base_path / "fitting",
        base_path / "projections-prep",
        base_path / "projections",
    ]

    for path in required_paths:
        if not path.exists():
            raise ValueError(f"Required path '{path}' does not exist")


def run_command(command, description=None, cwd=None):
    """Run a command and handle any errors."""
    if description:
        print(f"\n{description}...")
        print(f"Running: {command}")

    try:
        subprocess.run(command, shell=True, check=True, cwd=cwd)
        return True
    except subprocess.CalledProcessError:
        print(f"Error executing command: {command}", file=sys.stderr)
        return False


def generate_fitting_manifest(species_list, args):
    """Generate comprehensive fitting manifest with failure detection and ESS analysis."""
    base_path = _get_base_path()
    fitting_artefacts = base_path / "fitting" / "artefacts"
    
    species_to_analyze = _get_species_list(species_list)
    failed_batches = []
    successful_batches = []
    
    print(f"\n📋 FITTING MANIFEST - Batch {args.id}")
    print("=" * 50)
    
    for species in species_to_analyze:
        # Use correct file pattern matching our Docker scripts
        sigma_suffix = f"_sigma{args.amis_sigma}" if args.amis_sigma and args.amis_sigma != 0.0025 else ""
        amis_file = fitting_artefacts / "AMIS_output" / f"{species}_amis_output{args.id}{sigma_suffix}.Rdata"
        
        print(f"\n🔬 {species.upper()}")
        
        if amis_file.exists():
            # Try to analyze ESS for operational insights
            try:
                import subprocess
                r_check = f"""
                load("{amis_file}")
                if(exists("amis_output") && !is.null(amis_output$ess)) {{
                    ess_values <- amis_output$ess
                    cat("ESS_SUMMARY:", length(ess_values), "IUs,", 
                        "min:", round(min(ess_values), 1), 
                        "max:", round(max(ess_values), 1), 
                        "mean:", round(mean(ess_values), 1),
                        "below_threshold:", sum(ess_values < {args.ess_threshold}))
                }} else {{
                    cat("ESS_UNAVAILABLE")
                }}
                """
                
                result = subprocess.run(
                    ["Rscript", "-e", r_check], 
                    capture_output=True, text=True, timeout=10
                )
                
                if result.returncode == 0 and "ESS_SUMMARY:" in result.stdout:
                    ess_info = result.stdout.split("ESS_SUMMARY:")[1].strip()
                    print(f"   ✅ AMIS output: {amis_file.name}")
                    print(f"   📊 ESS: {ess_info}")
                    
                    # Check if any IUs need retry with higher sigma
                    if "below_threshold:" in ess_info and int(ess_info.split("below_threshold:")[-1].strip()) > 0:
                        below_count = int(ess_info.split("below_threshold:")[-1].strip())
                        print(f"   ⚠️  {below_count} IUs with ESS < {args.ess_threshold} - consider --amis-sigma=0.025")
                else:
                    print(f"   ✅ AMIS output: {amis_file.name}")
                    print(f"   ℹ️  ESS analysis unavailable")
                    
                successful_batches.append(f"{species}_{args.id}")
                
            except Exception:
                print(f"   ✅ AMIS output: {amis_file.name}")
                print(f"   ℹ️  ESS analysis skipped")
                successful_batches.append(f"{species}_{args.id}")
        else:
            print(f"   ❌ Missing: {amis_file.name}")
            print(f"   💡 Retry: --stage=fitting --id={args.id} --amis-sigma=0.025")
            failed_batches.append(f"{species}_{args.id}")
    
    # Provide actionable summary
    print(f"\n📊 SUMMARY")
    print(f"   ✅ Successful: {len(successful_batches)}")  
    print(f"   ❌ Failed: {len(failed_batches)}")
    
    if failed_batches:
        print(f"\n🔄 For find_lowESS_ids.R script:")
        failed_ids = [b.split('_')[-1] for b in failed_batches]
        print(f"   --failed-ids={','.join(failed_ids)}")
    
    return failed_batches


def generate_task_iu_mapping(base_path: Path, output_file="task_iu_mapping.csv"):
    """Generate a simple task mapping by calling R script directly."""
    print("\n🗺️  Generating Task ID to IU mapping...")
    
    # Simple R script to combine all lookup tables
    r_script = f"""
    library(dplyr)
    base_path <- "{base_path}"
    maps_path <- file.path(base_path, "fitting-prep/artefacts/Maps")
    
    all_data <- list()
    species_info <- list(
        list(name = "ascaris", file = "iu_task_lookup_sth.rds", iu_col = "IU_2021"),
        list(name = "hookworm", file = "iu_task_lookup_sth.rds", iu_col = "IU_2021"),
        list(name = "trichuris", file = "iu_task_lookup_trichuris.rds", iu_col = "IU_2021"),
        list(name = "haematobium", file = "iu_task_lookup_haema.rds", iu_col = "IU_ID"),
        list(name = "mansoni_low_burden", file = "iu_task_lookup_mansoni_low_burden.rds", iu_col = "IU_ID"),
        list(name = "mansoni_high_burden", file = "iu_task_lookup_mansoni_high_burden.rds", iu_col = "IU_ID")
    )
    
    for (info in species_info) {{
        lookup_file <- file.path(maps_path, info$file)
        if (file.exists(lookup_file)) {{
            load(lookup_file)
            df <- data.frame(
                TaskID = iu_task_lookup$TaskID,
                IU_CODE = iu_task_lookup[[info$iu_col]],
                Species = info$name,
                stringsAsFactors = FALSE
            )
            all_data[[info$name]] <- df
            cat("   ✓ Found", nrow(df), "entries for", info$name, "\\n")
        }} else {{
            cat("   ⚠️ Could not find", lookup_file, "\\n")
        }}
    }}
    
    if (length(all_data) > 0) {{
        combined_df <- do.call(rbind, all_data)
        output_path <- file.path(base_path, "fitting-prep/artefacts", "{output_file}")
        write.csv(combined_df, output_path, row.names = FALSE)
        cat("   📄 Mapping saved to:", "{output_file}", "\\n")
        cat("   📊 Total entries:", nrow(combined_df), "\\n")
        cat("   📈 Species coverage:", length(unique(combined_df$Species)), "species\\n")
        cat("SUCCESS\\n")
    }} else {{
        cat("ERROR: No lookup tables found\\n")
    }}
    """
    
    # Execute R script
    result = subprocess.run(
        ["Rscript", "-e", r_script], capture_output=True, text=True, cwd=str(base_path)
    )
    
    return "SUCCESS" in result.stdout


def _run_fitting_prep_sth_impl(species_list, batch_id, scripts_path):
    """Run fitting preparation for STH species."""

    commands = []

    # For STH maps preparation, we need BOTH ascaris/hookworm AND trichuris histories
    # because prepare_maps_allspecies.R reads both files
    needs_ascaris_hookworm = (
        STHSpecies.ASCARIS in species_list
        or STHSpecies.HOOKWORM in species_list
        or STHSpecies.ALL in species_list
    )
    needs_trichuris = (
        STHSpecies.TRICHURIS in species_list or STHSpecies.ALL in species_list
    )

    # Build batch argument for projection scripts
    batch_arg = f" --id {batch_id}" if batch_id else ""

    # First, prepare histories for fitting
    if needs_ascaris_hookworm:
        commands.append(
            ("Rscript prepare_histories.R", "Preparing histories for ascaris/hookworm")
        )

    if (
        needs_trichuris or needs_ascaris_hookworm
    ):  # Always run trichuris if any STH species needed
        commands.append(
            (
                "Rscript prepare_histories_trichuris.R",
                "Preparing histories for trichuris",
            )
        )

    # Always run maps preparation after histories
    if commands:
        commands.append(
            ("Rscript prepare_maps_allspecies.R", "Preparing maps for all STH species")
        )

    # Also prepare projection lookup tables (table_iu_idx files)
    if needs_ascaris_hookworm:
        commands.append(
            (
                f"Rscript prepare_histories_projections.R --species ascaris{batch_arg}",
                "Preparing projection lookup tables for ascaris",
            )
        )
        commands.append(
            (
                f"Rscript prepare_histories_projections.R --species hookworm{batch_arg}",
                "Preparing projection lookup tables for hookworm",
            )
        )

    if needs_trichuris:
        commands.append(
            (
                f"Rscript prepare_histories_trichuris_projections.R --species trichuris{batch_arg}",
                "Preparing projection lookup tables for trichuris",
            )
        )

    for cmd, desc in commands:
        if not run_command(cmd, desc, cwd=str(scripts_path)):
            return False

    return True


def _run_fitting_prep_sch_impl(species_list, batch_id, scripts_path):
    """Run fitting preparation for SCH species."""

    commands = []

    batch_arg = f" --id {batch_id}" if batch_id else ""

    # Prepare histories and maps for fitting
    if SCHSpecies.HAEMATOBIUM in species_list or SCHSpecies.ALL in species_list:
        commands.append(
            (
                f"Rscript prepare_histories_and_maps_haematobium.R{batch_arg}",
                "Preparing histories and maps for haematobium",
            )
        )

    if (
        SCHSpecies.MANSONI_HIGH in species_list
        or SCHSpecies.MANSONI_LOW in species_list
        or SCHSpecies.ALL in species_list
    ):
        commands.append(
            (
                f"Rscript prepare_histories_and_maps_mansoni.R{batch_arg}",
                "Preparing histories and maps for mansoni",
            )
        )

    # Also prepare projection lookup tables (table_iu_idx files)
    if species_list:  # If any SCH species requested
        # Determine species argument for SCH projection script
        if SCHSpecies.HAEMATOBIUM in species_list or SCHSpecies.ALL in species_list:
            species_arg = "haematobium"
        elif SCHSpecies.MANSONI_HIGH in species_list:
            species_arg = "mansoni_high_burden"
        elif SCHSpecies.MANSONI_LOW in species_list:
            species_arg = "mansoni_low_burden"
        else:
            species_arg = "haematobium"  # default

        commands.append(
            (
                f"Rscript prepare_histories_projections_sch.R --species {species_arg}{batch_arg}",
                f"Preparing projection lookup tables for SCH species: {species_arg}",
            )
        )

    for cmd, desc in commands:
        if not run_command(cmd, desc, cwd=str(scripts_path)):
            return False

    return True


def build_amis_args(args):
    """Build AMIS parameter arguments from command line args."""
    amis_args = []

    if args.amis_sigma is not None:
        amis_args.extend(["--amis-sigma", str(args.amis_sigma)])

    if args.amis_n_samples is not None:
        amis_args.extend(["--amis-n-samples", str(args.amis_n_samples)])

    if args.amis_target_ess is not None:
        amis_args.extend(["--amis-target-ess", str(args.amis_target_ess)])

    if args.amis_n_iters is not None:
        amis_args.extend(["--amis-n-iters", str(args.amis_n_iters)])

    if args.num_cores is not None:
        amis_args.extend(["--num-cores", str(args.num_cores)])

    return amis_args


def run_fitting(sth_species_list, sch_species_list, args):
    """Run fitting for both STH and SCH species in a unified way."""
    base_path = _get_base_path()
    scripts_path = base_path / "fitting" / "scripts"

    # Process both disease types with their respective scripts
    disease_configs = []

    if sth_species_list:
        disease_configs.append(
            {
                "species_list": sth_species_list,
                "script": "sth_fitting.R",
                "disease_type": "STH",
                "is_sth": True,
            }
        )

    if sch_species_list:
        disease_configs.append(
            {
                "species_list": sch_species_list,
                "script": "sch_fitting.R",
                "disease_type": "SCH",
                "is_sth": False,
            }
        )

    for config in disease_configs:
        # Get species strings from enum list
        species_to_fit = _get_species_list(config["species_list"])

        # Fit each species
        for species in species_to_fit:
            print(f"\nFitting {config['disease_type']} species: {species}")

            # Build command with base arguments
            cmd_parts = [
                "Rscript",
                config["script"],
                "--id",
                str(args.id),
                "--species",
                species,
            ]

            # Add AMIS parameters
            cmd_parts.extend(build_amis_args(args))
            cmd = " ".join(cmd_parts)

            # Run fitting command
            if not run_command(
                cmd, f"Running AMIS fitting for {species}", cwd=str(scripts_path)
            ):
                print(
                    f"✗ Fitting failed for {species}. Check ESS_NOT_REACHED_{species}.txt for details.",
                    file=sys.stderr,
                )
                print(
                    f"  Consider re-running with --amis-sigma=0.025 for failed batches",
                    file=sys.stderr,
                )
                return False

            print(f"✓ Fitting completed for {species}")

        # Generate fitting manifests
        if species_to_fit:  # Only generate if we processed species
            generate_fitting_manifest(config["species_list"], args)

    return True


def run_projections_prep(sth_species_list, sch_species_list, args):
    """Run projections preparation for both STH and SCH species in a unified way."""
    base_path = _get_base_path()
    scripts_path = base_path / "projections-prep" / "scripts"

    # Process both disease types with the same script
    disease_configs = []

    if sth_species_list:
        disease_configs.append(
            {"species_list": sth_species_list, "disease_type": "STH"}
        )

    if sch_species_list:
        disease_configs.append(
            {"species_list": sch_species_list, "disease_type": "SCH"}
        )

    for config in disease_configs:
        # Get species strings from enum list
        species_to_process = _get_species_list(config["species_list"])

        # Process each species
        for species in species_to_process:
            print(
                f"\nProcessing projections prep for {config['disease_type']} species: {species}"
            )

            # Build command with base arguments
            cmd_parts = [
                "Rscript",
                "preprocess_for_projections.R",
                "--species",
                species,
            ]

            # Add single batch ID if provided
            if args.id:
                cmd_parts.extend(["--id", str(args.id)])

            # Add failed IDs if provided
            if args.failed_ids:
                cmd_parts.extend(["--failed-ids", args.failed_ids])

            # Add ESS threshold if provided
            if args.ess_threshold:
                cmd_parts.extend(["--ess-threshold", str(args.ess_threshold)])

            # Add AMIS sigma if provided
            if args.amis_sigma:
                cmd_parts.extend(["--amis-sigma", str(args.amis_sigma)])

            cmd = " ".join(cmd_parts)

            # Run preprocessing command
            if not run_command(
                cmd,
                f"Running projections preprocessing for {species}",
                cwd=str(scripts_path),
            ):
                print(f"✗ Projections prep failed for {species}", file=sys.stderr)
                return False

            print(f"✓ Projections prep completed for {species}")

    return True


def run_nearterm_projections(sth_species_list, sch_species_list, args):
    """Run near-term projections for both STH and SCH species in a unified way."""
    base_path = _get_base_path()
    scripts_path = base_path / "projections" / "scripts"

    # Validate that batch ID is provided (required for projections)
    if not args.id:
        print(
            "Error: --id is required for nearterm-projections stage",
            file=sys.stderr,
        )
        return False

    # Process both disease types with their respective scripts
    disease_configs = []

    if sth_species_list:
        disease_configs.append(
            {
                "species_list": sth_species_list,
                "script": "sth_projections_per_IU.py",
                "disease_type": "STH",
            }
        )

    if sch_species_list:
        disease_configs.append(
            {
                "species_list": sch_species_list,
                "script": "sch_projections_per_IU.py",
                "disease_type": "SCH",
            }
        )

    for config in disease_configs:
        # Get species strings from enum list
        species_to_process = _get_species_list(config["species_list"])

        # Process each species
        for species in species_to_process:
            print(
                f"\nRunning near-term projections for {config['disease_type']} species: {species}"
            )

            # Build command with base arguments
            cmd_parts = ["python", config["script"], "--species", species]

            # Add batch ID - required for projections
            cmd_parts.extend(["--id", str(args.id)])

            # Add number of cores if provided
            if args.num_cores:
                cmd_parts.extend(["--num-cores", str(args.num_cores)])

            cmd = " ".join(cmd_parts)

            # Run projections command
            if not run_command(
                cmd,
                f"Running near-term projections for {species}",
                cwd=str(scripts_path),
            ):
                print(f"✗ Near-term projections failed for {species}", file=sys.stderr)
                return False

            print(f"✓ Near-term projections completed for {species}")

    return True


def execute_single_stage(stage, sth_species, sch_species, args):
    """Execute a single pipeline stage with unified validation and error handling."""
    # Stage function mapping
    stage_functions = {
        Stage.FITTING_PREP: lambda: run_fitting_prep(sth_species, sch_species, args),
        Stage.FITTING: lambda: run_fitting(sth_species, sch_species, args),
        Stage.PROJECTIONS_PREP: lambda: run_projections_prep(sth_species, sch_species, args),
        Stage.NEARTERM_PROJECTIONS: lambda: run_nearterm_projections(sth_species, sch_species, args),
    }
    
    # Success messages
    success_messages = {
        Stage.FITTING_PREP: "Fitting preparation completed successfully!",
        Stage.FITTING: "Fitting completed successfully!",
        Stage.PROJECTIONS_PREP: "Projections preparation completed successfully!",
        Stage.NEARTERM_PROJECTIONS: "Near-term projections completed successfully!",
    }
    
    # Validate dependencies (except for fitting-prep)
    if stage != Stage.FITTING_PREP:
        all_species = sth_species + sch_species
        if not validate_stage_dependencies(stage, all_species, args.id, args.amis_sigma):
            print(f"\n❌ Dependencies not satisfied for stage: {stage.value}", file=sys.stderr)
            return False
    
    # Validate required arguments
    if stage in [Stage.FITTING, Stage.NEARTERM_PROJECTIONS] and args.id is None:
        print(f"Error: --id is required for {stage.value} stage", file=sys.stderr)
        return False
    
    # Execute the stage
    if stage_functions[stage]():
        print(f"\n{success_messages[stage]}")
        return True
    else:
        print(f"\n{stage.value.title()} failed!", file=sys.stderr)
        return False


def execute_pipeline_stages(stages, sth_species, sch_species, args):
    """Execute a sequence of pipeline stages with unified logic."""
    stage_names = {
        Stage.FITTING_PREP: "FITTING-PREP",
        Stage.FITTING: "FITTING", 
        Stage.PROJECTIONS_PREP: "PROJECTIONS-PREP",
        Stage.NEARTERM_PROJECTIONS: "NEAR-TERM PROJECTIONS",
    }

    for i, stage in enumerate(stages, 1):
        print(f"\n{'='*60}")
        print(f"STAGE {i}: {stage_names[stage]}")
        print(f"{'='*60}")

        if not execute_single_stage(stage, sth_species, sch_species, args):
            return False

    return True


def run_fitting_prep(sth_species_list, sch_species_list, args):
    """Run fitting preparation for both STH and SCH species."""
    base_path = _get_base_path()
    scripts_path = base_path / "fitting-prep" / "scripts"
    success = True

    if sth_species_list:
        print(
            f"\nRunning fitting-prep for STH species: {[s.value for s in sth_species_list]}"
        )
        success = success and _run_fitting_prep_sth_impl(
            sth_species_list, args.id, scripts_path
        )

    if sch_species_list:
        print(
            f"\nRunning fitting-prep for SCH species: {[s.value for s in sch_species_list]}"
        )
        success = success and _run_fitting_prep_sch_impl(
            sch_species_list, args.id, scripts_path
        )

    if success:
        # Generate comprehensive Task ID to IU mapping
        generate_task_iu_mapping(base_path)

    return success


def main():
    parser = argparse.ArgumentParser(
        description="Run STH/SCH fitting and projections pipeline stages",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run fitting-prep for all STH species
  docker run --rm -v $(pwd):/ntdmc/sth-sch-amis-integration sch-sth-amis-pipeline --stage=fitting-prep --disease=sth
  
  # Run fitting-prep for specific SCH species
  docker run --rm -v $(pwd):/ntdmc/sth-sch-amis-integration sch-sth-amis-pipeline --stage=fitting-prep --disease=sch --sch-species=haematobium
  
  # Run fitting-prep for multiple STH species
  docker run --rm -v $(pwd):/ntdmc/sth-sch-amis-integration sch-sth-amis-pipeline --stage=fitting-prep --disease=sth --sth-species=ascaris,hookworm
""",
    )

    # Basic arguments
    parser.add_argument(
        "--stage",
        type=str,
        choices=[s.value for s in Stage],
        required=False,
        default=Stage.ALL.value,
        help="Pipeline stage to run",
    )

    parser.add_argument(
        "--disease",
        type=str,
        choices=[d.value for d in Disease],
        default=Disease.ALL.value,
        help="Disease type to process (default: all)",
    )

    parser.add_argument(
        "--sth-species",
        type=str,
        help="STH species to process - can be comma-separated (e.g., 'ascaris,hookworm') or space-separated (default: all)",
    )

    parser.add_argument(
        "--sch-species",
        type=str,
        help="SCH species to process - can be comma-separated (e.g., 'haematobium,mansoni_high_burden') or space-separated (default: all)",
    )

    # Stage-specific arguments
    parser.add_argument("--id", type=int, help="Batch/Task ID")
    parser.add_argument(
        "--failed-ids", type=str, help="Comma-separated list of failed IDs"
    )
    parser.add_argument("--amis-sigma", type=float, help="AMIS sigma parameter")
    parser.add_argument("--amis-n-samples", type=int, help="AMIS number of samples")
    parser.add_argument("--amis-target-ess", type=int, help="AMIS target ESS")
    parser.add_argument("--amis-n-iters", type=int, help="AMIS maximum iterations")
    parser.add_argument("--num-cores", type=int, help="Number of cores to use")
    parser.add_argument("--ess-threshold", type=int, default=200, help="ESS threshold")

    args = parser.parse_args()

    try:
        validate_environment()
    except ValueError as e:
        print(f"Environment validation failed: {e}", file=sys.stderr)
        sys.exit(1)

    # Parse disease and species selections
    disease = Disease(args.disease)
    sth_species, sch_species = parse_species_arguments(args, disease)

    # Execute the requested stage
    stage = Stage(args.stage)

    # Execute single stages or multi-stage pipelines
    if stage in [Stage.FITTING_PREP, Stage.FITTING, Stage.PROJECTIONS_PREP, Stage.NEARTERM_PROJECTIONS]:
        if not execute_single_stage(stage, sth_species, sch_species, args):
            sys.exit(1)

    elif stage == Stage.ALL:
        print(
            "Running complete pipeline: fitting-prep → fitting → projections-prep → nearterm-projections"
        )

        all_stages = [
            Stage.FITTING_PREP,
            Stage.FITTING,
            Stage.PROJECTIONS_PREP,
            Stage.NEARTERM_PROJECTIONS,
        ]

        if execute_pipeline_stages(all_stages, sth_species, sch_species, args):
            print("✅ Complete pipeline executed successfully!")
            print(f"\n{'='*60}")
            print("COMPLETE PIPELINE SUMMARY")
            print(f"{'='*60}")
            print("✅ fitting-prep: COMPLETED")
            print("✅ fitting: COMPLETED")
            print("✅ projections-prep: COMPLETED")
            print("✅ nearterm-projections: COMPLETED")
            print(f"{'='*60}")
        else:
            sys.exit(1)

    elif stage == Stage.SKIP_FITTING_PREP:
        print(
            "Running pipeline (skipping fitting-prep): fitting → projections-prep → nearterm-projections"
        )

        # Validate that fitting-prep artifacts exist before proceeding
        all_species = sth_species + sch_species
        if not validate_stage_dependencies(Stage.FITTING, all_species, args.id):
            print(
                "\n❌ Cannot skip fitting-prep: required artifacts are missing",
                file=sys.stderr,
            )
            print(
                "   Either run --stage=fitting-prep first, or use pre-built artifacts",
                file=sys.stderr,
            )
            sys.exit(1)

        skip_stages = [
            Stage.FITTING,
            Stage.PROJECTIONS_PREP,
            Stage.NEARTERM_PROJECTIONS,
        ]

        if execute_pipeline_stages(skip_stages, sth_species, sch_species, args):
            print("✅ Pipeline executed successfully (fitting-prep skipped)!")
            print(f"\n{'='*60}")
            print("PIPELINE SUMMARY (fitting-prep SKIPPED)")
            print(f"{'='*60}")
            print("⏭️  fitting-prep: SKIPPED (using pre-built artifacts)")
            print("✅ fitting: COMPLETED")
            print("✅ projections-prep: COMPLETED")
            print("✅ nearterm-projections: COMPLETED")
            print(f"{'='*60}")
        else:
            sys.exit(1)


if __name__ == "__main__":
    main()
