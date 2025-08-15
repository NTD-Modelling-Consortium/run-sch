#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import enum
import json
import glob
import tempfile
from datetime import datetime
from pathlib import Path
from io import StringIO
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


def _get_pipeline_dependencies() -> Dict[Stage, StageDependency]:
    """Define the complete pipeline dependency graph."""
    base_path = _get_base_path()

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
        result = subprocess.run(command, shell=True, check=True, cwd=cwd)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {command}", file=sys.stderr)
        return False


def generate_fitting_manifest(species_list, args, is_sth=True):
    """Generate a manifest file after fitting stage completion."""
    base_path = _get_base_path()
    fitting_artefacts = base_path / "fitting" / "artefacts"

    # Determine species to analyze
    if is_sth:
        if STHSpecies.ALL in species_list:
            species_to_analyze = ["ascaris", "hookworm", "trichuris"]
        else:
            species_to_analyze = [s.value for s in species_list]
    else:
        if SCHSpecies.ALL in species_list:
            species_to_analyze = [
                "haematobium",
                "mansoni_low_burden",
                "mansoni_high_burden",
            ]
        else:
            species_to_analyze = [s.value for s in species_list]

    all_manifests = []

    for species in species_to_analyze:
        print(f"\nGenerating fitting manifest for {species}...")

        manifest = {
            "species": species,
            "sigma": str(args.amis_sigma) if args.amis_sigma else "0.0025",
            "timestamp": datetime.now().isoformat(),
            "batch_id": args.id,
            "ius_in_batch": [],
            "batches_processed": [],
            "successful_batches": [],
            "failed_batches": [],
            "low_ess_batches": [],
            "missing_output_files": [],
            "recommendations": {},
            "summary": {},
        }

        # Look for AMIS output files for this species and batch
        sigma_str = str(args.amis_sigma) if args.amis_sigma else "0.0025"
        amis_pattern = (
            fitting_artefacts / f"fit_amis_{species}_{args.id}_sigma{sigma_str}.RData"
        )
        amis_files = glob.glob(str(amis_pattern))

        if amis_files:
            manifest["batches_processed"].append(args.id)
            manifest["successful_batches"].append(args.id)

            # Try to read ESS information from R using subprocess
            try:
                # Create a temporary R script to read AMIS output and extract ESS
                r_script = f"""
                library("AMISforInfectiousDiseases")
                
                # Load lookup table to get IU count expectation
                base_path <- "{base_path}"
                kPathToMapsArtefacts <- file.path(base_path, "fitting-prep/artefacts/Maps")
                
                if("{species}" == "trichuris") {{
                  load(file.path(kPathToMapsArtefacts, "iu_task_lookup_trichuris.rds"))
                }} else if ("{species}" %in% c("ascaris","hookworm")) {{
                  load(file.path(kPathToMapsArtefacts, "iu_task_lookup_sth.rds"))  
                }} else if ("{species}" == "haematobium") {{
                  load(file.path(kPathToMapsArtefacts, "iu_task_lookup_haema.rds"))
                }} else {{
                  lookup_file <- paste0("iu_task_lookup_", "{species}", ".rds")
                  load(file.path(kPathToMapsArtefacts, lookup_file))
                }}
                
                # Get expected IU count for this batch
                if("{species}" %in% c("ascaris", "hookworm", "trichuris")) {{
                  expected_ius <- iu_task_lookup[iu_task_lookup$TaskID == {args.id}, "IU_2021"]
                }} else {{
                  expected_ius <- iu_task_lookup[iu_task_lookup$TaskID == {args.id}, "IU_ID"]
                }}
                expected_count <- length(expected_ius)
                
                # Load AMIS output
                amis_file <- "{amis_files[0]}"
                if(file.exists(amis_file)) {{
                  load(amis_file)  # loads amis_output
                  
                  if(exists("amis_output") && !is.null(amis_output$ess)) {{
                    ess_values <- amis_output$ess
                    min_ess <- min(ess_values, na.rm=TRUE)
                    max_ess <- max(ess_values, na.rm=TRUE)
                    avg_ess <- mean(ess_values, na.rm=TRUE)
                    low_ess_count <- sum(ess_values < 200, na.rm=TRUE)
                    total_ius <- length(ess_values)
                    
                    cat("ESS_ANALYSIS:", min_ess, max_ess, avg_ess, low_ess_count, total_ius, expected_count, "\\n")
                  }} else {{
                    cat("ESS_ANALYSIS: ERROR - No ESS data found\\n")
                  }}
                }} else {{
                  cat("ESS_ANALYSIS: ERROR - File not found\\n")
                }}
                """

                # Write and execute R script using temporary file
                with tempfile.NamedTemporaryFile(
                    mode="w",
                    suffix=".R",
                    prefix=f"temp_ess_analysis_{species}_",
                ) as temp_file:
                    temp_file.write(r_script)
                    temp_file.flush()

                    result = subprocess.run(
                        ["Rscript", temp_file.name],
                        capture_output=True,
                        text=True,
                        cwd=str(fitting_artefacts),
                    )

                if result.returncode == 0:
                    # Parse ESS analysis output
                    for line in result.stdout.split("\n"):
                        if line.startswith("ESS_ANALYSIS:"):
                            parts = line.split()
                            if len(parts) >= 7 and parts[1] != "ERROR":
                                min_ess = float(parts[1])
                                max_ess = float(parts[2])
                                avg_ess = float(parts[3])
                                low_ess_count = int(parts[4])
                                total_ius = int(parts[5])
                                expected_count = int(parts[6])

                                manifest["summary"] = {
                                    "total_ius_processed": total_ius,
                                    "expected_ius": expected_count,
                                    "processing_complete": total_ius == expected_count,
                                    "min_ess": round(min_ess, 2),
                                    "max_ess": round(max_ess, 2),
                                    "avg_ess": round(avg_ess, 2),
                                    "low_ess_ius": low_ess_count,
                                    "low_ess_percentage": (
                                        round(100 * low_ess_count / total_ius, 1)
                                        if total_ius > 0
                                        else 0
                                    ),
                                }

                                if low_ess_count > 0:
                                    manifest["low_ess_batches"].append(args.id)
                                    manifest["recommendations"][
                                        "rerun_with_higher_sigma"
                                    ] = [args.id]
                                    manifest["recommendations"][
                                        "suggested_sigma"
                                    ] = "0.025"
                                    manifest["recommendations"][
                                        "reason"
                                    ] = f"{low_ess_count} IUs have ESS < 200"

                                break
                            elif "ERROR" in parts:
                                manifest["summary"][
                                    "error"
                                ] = "Could not analyze ESS data"
                                break

            except Exception as e:
                manifest["summary"]["error"] = f"ESS analysis failed: {str(e)}"
        else:
            # No output files found - batch failed
            manifest["failed_batches"].append(args.id)
            manifest["missing_output_files"].append(
                f"fit_amis_{species}_{args.id}_sigma{sigma_str}.RData"
            )
            manifest["recommendations"]["rerun_batch"] = [args.id]
            manifest["recommendations"][
                "reason"
            ] = "No AMIS output file found - batch likely failed"

        # Set final recommendations
        if manifest["failed_batches"] or manifest["low_ess_batches"]:
            if manifest["failed_batches"]:
                manifest["recommendations"]["action"] = "rerun_failed_batches"
                manifest["recommendations"]["suggested_sigma"] = "0.025"
            elif manifest["low_ess_batches"]:
                manifest["recommendations"]["action"] = "rerun_low_ess_batches"
                manifest["recommendations"]["suggested_sigma"] = "0.025"
        else:
            manifest["recommendations"]["action"] = "none_needed"
            manifest["recommendations"][
                "message"
            ] = "All batches completed successfully with adequate ESS"

        # Save manifest
        manifest_file = (
            fitting_artefacts
            / f"fitting_manifest_{species}_{args.id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        )
        manifest_file.write_text(json.dumps(manifest, indent=2))

        all_manifests.append(manifest)

        # Print summary to user
        print(f"📊 Fitting Summary for {species} (Batch {args.id}):")
        if manifest["successful_batches"]:
            if "total_ius_processed" in manifest["summary"]:
                summary = manifest["summary"]
                print(
                    f"   ✓ Processed {summary['total_ius_processed']} IUs successfully"
                )
                print(
                    f"   📈 ESS: min={summary['min_ess']}, max={summary['max_ess']}, avg={summary['avg_ess']}"
                )
                if summary["low_ess_ius"] > 0:
                    print(
                        f"   ⚠️  {summary['low_ess_ius']} IUs ({summary['low_ess_percentage']}%) have ESS < 200"
                    )
                else:
                    print(f"   ✅ All IUs have adequate ESS (≥200)")
            else:
                print(f"   ✓ AMIS output file generated")

        if manifest["failed_batches"]:
            print(f"   ❌ Batch failed - no output file generated")

        if manifest["recommendations"]["action"] != "none_needed":
            print(f"   💡 Recommendation: {manifest['recommendations']['action']}")
            if "suggested_sigma" in manifest["recommendations"]:
                print(
                    f"      Retry with: --amis-sigma={manifest['recommendations']['suggested_sigma']}"
                )

        print(f"   📄 Manifest saved: {manifest_file.name}")

    return all_manifests


def load_species_lookup_table(species, base_path):
    """Load the appropriate lookup table for the given species."""
    base_path = Path(base_path)
    maps_path = base_path / "fitting-prep" / "artefacts" / "Maps"

    # Determine which lookup file to use based on species
    if species == "trichuris":
        lookup_file = maps_path / "iu_task_lookup_trichuris.rds"
    elif species in ["ascaris", "hookworm"]:
        lookup_file = maps_path / "iu_task_lookup_sth.rds"
    elif species == "haematobium":
        lookup_file = maps_path / "iu_task_lookup_haema.rds"
    elif species in ["mansoni_low_burden", "mansoni_high_burden"]:
        lookup_file = maps_path / f"iu_task_lookup_{species}.rds"
    else:
        raise ValueError(f"Unsupported species: {species}")

    if not lookup_file.exists():
        raise FileNotFoundError(f"Lookup file not found: {lookup_file}")

    # Use R to read the RDS file and convert to CSV
    r_script = f"""
    load("{lookup_file}")
    
    # Standardize column names
    if("{species}" %in% c("ascaris", "hookworm", "trichuris")) {{
        # STH species use IU_2021 column
        if(!"IU_CODE" %in% colnames(iu_task_lookup)) {{
            iu_task_lookup$IU_CODE <- iu_task_lookup$IU_2021
        }}
    }} else {{
        # SCH species use IU_ID column  
        if(!"IU_CODE" %in% colnames(iu_task_lookup)) {{
            iu_task_lookup$IU_CODE <- iu_task_lookup$IU_ID
        }}
    }}
    
    # Output as CSV to stdout
    write.csv(iu_task_lookup, stdout(), row.names=FALSE)
    """

    # Execute R script
    result = subprocess.run(
        ["Rscript", "-e", r_script], capture_output=True, text=True, cwd=str(base_path)
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to read lookup table for {species}: {result.stderr}"
        )

    # Parse CSV output
    import pandas as pd

    df = pd.read_csv(StringIO(result.stdout))
    # Add species column
    df["Species"] = species
    return df


def generate_task_iu_mapping(base_path: Path, output_file="task_iu_mapping.csv"):
    """Generate comprehensive TaskID to IU mapping for all species."""

    all_species = [
        "ascaris",
        "hookworm",
        "trichuris",  # STH species
        "haematobium",
        "mansoni_low_burden",
        "mansoni_high_burden",  # SCH species
    ]

    print("\n🗺️  Generating comprehensive Task ID to IU mapping...")

    all_mappings = []

    for species in all_species:
        try:
            print(f"   Processing {species}...")
            df = load_species_lookup_table(species, str(base_path))
            all_mappings.append(df)
            print(f"     ✓ Found {len(df)} entries for {species}")
        except (FileNotFoundError, RuntimeError) as e:
            print(f"     ⚠️ Could not process {species}: {e}")
            continue

    if not all_mappings:
        print("   ❌ No lookup tables could be processed")
        return False

    try:
        import pandas as pd

        # Combine all mappings using pandas
        combined_df = pd.concat(all_mappings, ignore_index=True)

        # Standardize column names and order
        columns_to_keep = ["TaskID", "IU_CODE", "Species"]

        # Add additional columns if they exist
        if "Country" in combined_df.columns:
            columns_to_keep.append("Country")
        elif "country" in combined_df.columns:
            combined_df["Country"] = combined_df["country"]
            columns_to_keep.append("Country")

        # Keep only relevant columns
        available_columns = [
            col for col in columns_to_keep if col in combined_df.columns
        ]
        final_df = combined_df[available_columns].copy()

        # Sort by Species, then TaskID, then IU_CODE
        final_df = final_df.sort_values(["Species", "TaskID", "IU_CODE"]).reset_index(
            drop=True
        )

        # Save to CSV
        output_path = base_path / "fitting-prep" / "artefacts" / output_file
        final_df.to_csv(str(output_path), index=False)

        print(f"   📄 Mapping saved to: {output_file}")
        print(f"   📊 Total entries: {len(final_df)}")
        print(
            f"   📈 Species coverage: {len(final_df['Species'].unique())} species, {final_df['TaskID'].nunique()} total batches"
        )

        return True

    except Exception as e:
        print(f"   ❌ Error generating mapping: {e}")
        return False


def _run_fitting_prep_sth_impl(species_list, batch_id, base_path, scripts_path):
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


def _run_fitting_prep_sch_impl(species_list, batch_id, base_path, scripts_path):
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
            generate_fitting_manifest(
                config["species_list"], args, is_sth=config["is_sth"]
            )

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


def execute_pipeline_stages(stages, sth_species, sch_species, args):
    """Execute a sequence of pipeline stages with unified logic."""
    # Stage function mapping
    stage_functions = {
        Stage.FITTING_PREP: lambda: run_fitting_prep(sth_species, sch_species, args),
        Stage.FITTING: lambda: run_fitting(sth_species, sch_species, args),
        Stage.PROJECTIONS_PREP: lambda: run_projections_prep(
            sth_species, sch_species, args
        ),
        Stage.NEARTERM_PROJECTIONS: lambda: run_nearterm_projections(
            sth_species, sch_species, args
        ),
    }

    # Stage display names
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

        # Validate dependencies (except for fitting-prep which has no deps)
        if stage != Stage.FITTING_PREP:
            all_species = sth_species + sch_species
            if not validate_stage_dependencies(
                stage, all_species, args.id, args.amis_sigma
            ):
                print(
                    f"\n❌ Dependencies not satisfied for stage: {stage.value}",
                    file=sys.stderr,
                )
                return False

        # Validate required arguments for specific stages
        if stage in [Stage.FITTING, Stage.NEARTERM_PROJECTIONS] and args.id is None:
            print(f"Error: --id is required for {stage.value} stage", file=sys.stderr)
            return False

        # Execute the stage
        if not stage_functions[stage]():
            print(f"\n{stage_names[stage]} failed!", file=sys.stderr)
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
            sth_species_list, args.id, base_path, scripts_path
        )

    if sch_species_list:
        print(
            f"\nRunning fitting-prep for SCH species: {[s.value for s in sch_species_list]}"
        )
        success = success and _run_fitting_prep_sch_impl(
            sch_species_list, args.id, base_path, scripts_path
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

    # Determine which species to process
    sth_species = []
    sch_species = []

    if disease in [Disease.STH, Disease.ALL]:
        if args.sth_species:
            # Handle comma-separated or space-separated species
            species_list = []
            if "," in args.sth_species:
                species_list = [s.strip() for s in args.sth_species.split(",")]
            else:
                species_list = args.sth_species.split()

            # Validate species choices
            valid_sth_species = [s.value for s in STHSpecies]
            for species in species_list:
                if species not in valid_sth_species:
                    print(
                        f"Error: Invalid STH species '{species}'. Valid choices: {valid_sth_species}",
                        file=sys.stderr,
                    )
                    sys.exit(1)

            sth_species = [STHSpecies(s) for s in species_list]
        else:
            sth_species = [STHSpecies.ALL]

    if disease in [Disease.SCH, Disease.ALL]:
        if args.sch_species:
            # Handle comma-separated or space-separated species
            species_list = []
            if "," in args.sch_species:
                species_list = [s.strip() for s in args.sch_species.split(",")]
            else:
                species_list = args.sch_species.split()

            # Validate species choices
            valid_sch_species = [s.value for s in SCHSpecies]
            for species in species_list:
                if species not in valid_sch_species:
                    print(
                        f"Error: Invalid SCH species '{species}'. Valid choices: {valid_sch_species}",
                        file=sys.stderr,
                    )
                    sys.exit(1)

            sch_species = [SCHSpecies(s) for s in species_list]
        else:
            sch_species = [SCHSpecies.ALL]

    # Execute the requested stage
    stage = Stage(args.stage)

    # Execute stages based on the selected stage
    if stage == Stage.FITTING_PREP:
        success = run_fitting_prep(sth_species, sch_species, args)
        if success:
            print("\nFitting preparation completed successfully!")
        else:
            print("\nFitting preparation failed!", file=sys.stderr)
            sys.exit(1)

    elif stage == Stage.FITTING:
        # Validate dependencies and required arguments
        all_species = sth_species + sch_species
        if not validate_stage_dependencies(
            stage, all_species, args.id, args.amis_sigma
        ):
            print(
                f"\n❌ Dependencies not satisfied for stage: {stage.value}",
                file=sys.stderr,
            )
            sys.exit(1)

        if args.id is None:
            print("Error: --id is required for fitting stage", file=sys.stderr)
            sys.exit(1)

        success = run_fitting(sth_species, sch_species, args)
        if success:
            print("\nFitting completed successfully!")
        else:
            print("\nFitting failed!", file=sys.stderr)
            sys.exit(1)

    elif stage == Stage.PROJECTIONS_PREP:
        # Validate dependencies
        all_species = sth_species + sch_species
        if not validate_stage_dependencies(
            stage, all_species, args.id, args.amis_sigma
        ):
            print(
                f"\n❌ Dependencies not satisfied for stage: {stage.value}",
                file=sys.stderr,
            )
            sys.exit(1)

        success = run_projections_prep(sth_species, sch_species, args)
        if success:
            print("\nProjections preparation completed successfully!")
        else:
            print("\nProjections preparation failed!", file=sys.stderr)
            sys.exit(1)

    elif stage == Stage.NEARTERM_PROJECTIONS:
        # Validate dependencies and required arguments
        all_species = sth_species + sch_species
        if not validate_stage_dependencies(
            stage, all_species, args.id, args.amis_sigma
        ):
            print(
                f"\n❌ Dependencies not satisfied for stage: {stage.value}",
                file=sys.stderr,
            )
            sys.exit(1)

        if args.id is None:
            print(
                "Error: --id is required for nearterm-projections stage",
                file=sys.stderr,
            )
            sys.exit(1)

        success = run_nearterm_projections(sth_species, sch_species, args)
        if success:
            print("\nNear-term projections completed successfully!")
        else:
            print("\nNear-term projections failed!", file=sys.stderr)
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
