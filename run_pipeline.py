#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import enum
from pathlib import Path


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


def validate_environment():
    """Validate required environment variables are set."""
    # Check if we're in a Docker container
    if os.path.exists('/.dockerenv'):
        # In Docker, validate paths exist
        required_paths = [
            "/ntdmc/sth-sch-amis-integration/fitting-prep",
            "/ntdmc/sth-sch-amis-integration/fitting",
            "/ntdmc/sth-sch-amis-integration/projections-prep",
            "/ntdmc/sth-sch-amis-integration/projections"
        ]
        for path in required_paths:
            if not os.path.exists(path):
                raise ValueError(f"Required path '{path}' does not exist in container")
    else:
        # Outside Docker, just check we're in the right directory
        if not os.path.exists("fitting-prep"):
            raise ValueError("This script should be run from the project root directory")


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


def run_fitting_prep_sth(species_list, batch_id=None):
    """Run fitting preparation for STH species."""
    base_path = "/ntdmc/sth-sch-amis-integration" if os.path.exists('/.dockerenv') else "."
    scripts_path = f"{base_path}/fitting-prep/scripts"
    
    commands = []
    
    # For STH maps preparation, we need BOTH ascaris/hookworm AND trichuris histories
    # because prepare_maps_allspecies.R reads both files
    needs_ascaris_hookworm = (STHSpecies.ASCARIS in species_list or 
                             STHSpecies.HOOKWORM in species_list or 
                             STHSpecies.ALL in species_list)
    needs_trichuris = (STHSpecies.TRICHURIS in species_list or 
                      STHSpecies.ALL in species_list)
    
    # Build batch argument for projection scripts
    batch_arg = f" --id {batch_id}" if batch_id else ""
    
    # First, prepare histories for fitting
    if needs_ascaris_hookworm:
        commands.append(("Rscript prepare_histories.R", "Preparing histories for ascaris/hookworm"))
    
    if needs_trichuris or needs_ascaris_hookworm:  # Always run trichuris if any STH species needed
        commands.append(("Rscript prepare_histories_trichuris.R", "Preparing histories for trichuris"))
    
    # Always run maps preparation after histories
    if commands:
        commands.append(("Rscript prepare_maps_allspecies.R", "Preparing maps for all STH species"))
    
    # Also prepare projection lookup tables (table_iu_idx files)
    if needs_ascaris_hookworm:
        commands.append((f"Rscript prepare_histories_projections.R --species ascaris{batch_arg}", "Preparing projection lookup tables for ascaris"))
        commands.append((f"Rscript prepare_histories_projections.R --species hookworm{batch_arg}", "Preparing projection lookup tables for hookworm"))
    
    if needs_trichuris:
        commands.append((f"Rscript prepare_histories_trichuris_projections.R --species trichuris{batch_arg}", "Preparing projection lookup tables for trichuris"))
    
    for cmd, desc in commands:
        if not run_command(cmd, desc, cwd=scripts_path):
            return False
    
    return True


def run_fitting_prep_sch(species_list, batch_id=None):
    """Run fitting preparation for SCH species."""
    base_path = "/ntdmc/sth-sch-amis-integration" if os.path.exists('/.dockerenv') else "."
    scripts_path = f"{base_path}/fitting-prep/scripts"
    
    commands = []
    
    batch_arg = f" --id {batch_id}" if batch_id else ""
    
    # Prepare histories and maps for fitting
    if SCHSpecies.HAEMATOBIUM in species_list or SCHSpecies.ALL in species_list:
        commands.append((f"Rscript prepare_histories_and_maps_haematobium.R{batch_arg}", "Preparing histories and maps for haematobium"))
    
    if SCHSpecies.MANSONI_HIGH in species_list or SCHSpecies.MANSONI_LOW in species_list or SCHSpecies.ALL in species_list:
        commands.append((f"Rscript prepare_histories_and_maps_mansoni.R{batch_arg}", "Preparing histories and maps for mansoni"))
    
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
        
        commands.append((f"Rscript prepare_histories_projections_sch.R --species {species_arg}{batch_arg}", 
                        f"Preparing projection lookup tables for SCH species: {species_arg}"))
    
    for cmd, desc in commands:
        if not run_command(cmd, desc, cwd=scripts_path):
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


def run_fitting_sth(species_list, args):
    """Run fitting for STH species."""
    base_path = "/ntdmc/sth-sch-amis-integration" if os.path.exists('/.dockerenv') else "."
    scripts_path = f"{base_path}/fitting/scripts"
    
    # Determine which species to fit
    species_to_fit = []
    if STHSpecies.ALL in species_list:
        species_to_fit = ["ascaris", "hookworm", "trichuris"]
    else:
        species_to_fit = [s.value for s in species_list]
    
    # Fit each species
    for species in species_to_fit:
        print(f"\nFitting STH species: {species}")
        
        # Build command with base arguments
        cmd_parts = [
            "Rscript", "sth_fitting.R",
            "--id", str(args.id),
            "--species", species
        ]
        
        # Add AMIS parameters
        cmd_parts.extend(build_amis_args(args))
        cmd = " ".join(cmd_parts)
        
        # Run fitting command
        if not run_command(cmd, f"Running AMIS fitting for {species}", cwd=scripts_path):
            print(f"✗ Fitting failed for {species}. Check ESS_NOT_REACHED_{species}.txt for details.", file=sys.stderr)
            print(f"  Consider re-running with --amis-sigma=0.025 for failed batches", file=sys.stderr)
            return False
        
        print(f"✓ Fitting completed for {species}")
    
    return True


def run_fitting_sch(species_list, args):
    """Run fitting for SCH species."""
    base_path = "/ntdmc/sth-sch-amis-integration" if os.path.exists('/.dockerenv') else "."
    scripts_path = f"{base_path}/fitting/scripts"
    
    # Determine which species to fit
    species_to_fit = []
    if SCHSpecies.ALL in species_list:
        species_to_fit = ["haematobium", "mansoni_low_burden", "mansoni_high_burden"]
    else:
        species_to_fit = [s.value for s in species_list]
    
    # Fit each species
    for species in species_to_fit:
        print(f"\nFitting SCH species: {species}")
        
        # Build command with base arguments
        cmd_parts = [
            "Rscript", "sch_fitting.R",
            "--id", str(args.id),
            "--species", species
        ]
        
        # Add AMIS parameters
        cmd_parts.extend(build_amis_args(args))
        cmd = " ".join(cmd_parts)
        
        # Run fitting command
        if not run_command(cmd, f"Running AMIS fitting for {species}", cwd=scripts_path):
            print(f"✗ Fitting failed for {species}. Check ESS_NOT_REACHED_{species}.txt for details.", file=sys.stderr)
            print(f"  Consider re-running with --amis-sigma=0.025 for failed batches", file=sys.stderr)
            return False
        
        print(f"✓ Fitting completed for {species}")
    
    return True


def run_projections_prep_sth(species_list, args):
    """Run projections preparation for STH species."""
    base_path = "/ntdmc/sth-sch-amis-integration" if os.path.exists('/.dockerenv') else "."
    scripts_path = f"{base_path}/projections-prep/scripts"
    
    # Determine which species to process
    species_to_process = []
    if STHSpecies.ALL in species_list:
        species_to_process = ["ascaris", "hookworm", "trichuris"]
    else:
        species_to_process = [s.value for s in species_list]
    
    # Process each species
    for species in species_to_process:
        print(f"\nProcessing projections prep for STH species: {species}")
        
        # Build command with base arguments
        cmd_parts = [
            "Rscript", "preprocess_for_projections.R",
            "--species", species
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
        if not run_command(cmd, f"Running projections preprocessing for {species}", cwd=scripts_path):
            print(f"✗ Projections prep failed for {species}", file=sys.stderr)
            return False
        
        print(f"✓ Projections prep completed for {species}")
    
    return True


def run_projections_prep_sch(species_list, args):
    """Run projections preparation for SCH species."""
    base_path = "/ntdmc/sth-sch-amis-integration" if os.path.exists('/.dockerenv') else "."
    scripts_path = f"{base_path}/projections-prep/scripts"
    
    # Determine which species to process
    species_to_process = []
    if SCHSpecies.ALL in species_list:
        species_to_process = ["haematobium", "mansoni_low_burden", "mansoni_high_burden"]
    else:
        species_to_process = [s.value for s in species_list]
    
    # Process each species
    for species in species_to_process:
        print(f"\nProcessing projections prep for SCH species: {species}")
        
        # Build command with base arguments
        cmd_parts = [
            "Rscript", "preprocess_for_projections.R",
            "--species", species
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
        if not run_command(cmd, f"Running projections preprocessing for {species}", cwd=scripts_path):
            print(f"✗ Projections prep failed for {species}", file=sys.stderr)
            return False
        
        print(f"✓ Projections prep completed for {species}")
    
    return True


def run_nearterm_projections_sth(species_list, args):
    """Run near-term projections for STH species."""
    base_path = "/ntdmc/sth-sch-amis-integration" if os.path.exists('/.dockerenv') else "."
    scripts_path = f"{base_path}/projections/scripts"
    
    # Determine which species to process
    species_to_process = []
    if STHSpecies.ALL in species_list:
        species_to_process = ["ascaris", "hookworm", "trichuris"]
    else:
        species_to_process = [s.value for s in species_list]
    
    # Process each species
    for species in species_to_process:
        print(f"\nRunning near-term projections for STH species: {species}")
        
        # Build command with base arguments
        cmd_parts = [
            "python", "sth_projections_per_IU.py",
            "--species", species
        ]
        
        # Add batch ID - required for projections
        if args.id:
            cmd_parts.extend(["--id", str(args.id)])
        else:
            print("Error: --id is required for nearterm-projections stage", file=sys.stderr)
            return False
        
        # Add number of cores if provided
        if args.num_cores:
            cmd_parts.extend(["--num-cores", str(args.num_cores)])
        
        cmd = " ".join(cmd_parts)
        
        # Run projections command
        if not run_command(cmd, f"Running near-term projections for {species}", cwd=scripts_path):
            print(f"✗ Near-term projections failed for {species}", file=sys.stderr)
            return False
        
        print(f"✓ Near-term projections completed for {species}")
    
    return True


def run_nearterm_projections_sch(species_list, args):
    """Run near-term projections for SCH species."""
    base_path = "/ntdmc/sth-sch-amis-integration" if os.path.exists('/.dockerenv') else "."
    scripts_path = f"{base_path}/projections/scripts"
    
    # Determine which species to process
    species_to_process = []
    if SCHSpecies.ALL in species_list:
        species_to_process = ["haematobium", "mansoni_low_burden", "mansoni_high_burden"]
    else:
        species_to_process = [s.value for s in species_list]
    
    # Process each species
    for species in species_to_process:
        print(f"\nRunning near-term projections for SCH species: {species}")
        
        # Build command with base arguments
        cmd_parts = [
            "python", "sch_projections_per_IU.py",
            "--species", species
        ]
        
        # Add batch ID - required for projections
        if args.id:
            cmd_parts.extend(["--id", str(args.id)])
        else:
            print("Error: --id is required for nearterm-projections stage", file=sys.stderr)
            return False
        
        # Add number of cores if provided
        if args.num_cores:
            cmd_parts.extend(["--num-cores", str(args.num_cores)])
        
        cmd = " ".join(cmd_parts)
        
        # Run projections command
        if not run_command(cmd, f"Running near-term projections for {species}", cwd=scripts_path):
            print(f"✗ Near-term projections failed for {species}", file=sys.stderr)
            return False
        
        print(f"✓ Near-term projections completed for {species}")
    
    return True


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
"""
    )
    
    # Basic arguments
    parser.add_argument(
        "--stage",
        type=str,
        choices=[s.value for s in Stage],
        required=False,
        default=Stage.ALL.value,
        help="Pipeline stage to run"
    )
    
    parser.add_argument(
        "--disease",
        type=str,
        choices=[d.value for d in Disease],
        default=Disease.ALL.value,
        help="Disease type to process (default: all)"
    )
    
    parser.add_argument(
        "--sth-species",
        type=str,
        help="STH species to process - can be comma-separated (e.g., 'ascaris,hookworm') or space-separated (default: all)"
    )
    
    parser.add_argument(
        "--sch-species",
        type=str,
        help="SCH species to process - can be comma-separated (e.g., 'haematobium,mansoni_high_burden') or space-separated (default: all)"
    )
    
    
    # Stage-specific arguments
    parser.add_argument("--id", type=int, help="Batch/Task ID")
    parser.add_argument("--failed-ids", type=str, help="Comma-separated list of failed IDs")
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
            if ',' in args.sth_species:
                species_list = [s.strip() for s in args.sth_species.split(',')]
            else:
                species_list = args.sth_species.split()
            
            # Validate species choices
            valid_sth_species = [s.value for s in STHSpecies]
            for species in species_list:
                if species not in valid_sth_species:
                    print(f"Error: Invalid STH species '{species}'. Valid choices: {valid_sth_species}", file=sys.stderr)
                    sys.exit(1)
            
            sth_species = [STHSpecies(s) for s in species_list]
        else:
            sth_species = [STHSpecies.ALL]
    
    if disease in [Disease.SCH, Disease.ALL]:
        if args.sch_species:
            # Handle comma-separated or space-separated species
            species_list = []
            if ',' in args.sch_species:
                species_list = [s.strip() for s in args.sch_species.split(',')]
            else:
                species_list = args.sch_species.split()
            
            # Validate species choices
            valid_sch_species = [s.value for s in SCHSpecies]
            for species in species_list:
                if species not in valid_sch_species:
                    print(f"Error: Invalid SCH species '{species}'. Valid choices: {valid_sch_species}", file=sys.stderr)
                    sys.exit(1)
            
            sch_species = [SCHSpecies(s) for s in species_list]
        else:
            sch_species = [SCHSpecies.ALL]
    
    # Execute the requested stage
    stage = Stage(args.stage)
    
    if stage == Stage.FITTING_PREP:
        success = True
        
        if sth_species:
            print(f"\nRunning fitting-prep for STH species: {[s.value for s in sth_species]}")
            success = success and run_fitting_prep_sth(sth_species, args.id)
        
        if sch_species:
            print(f"\nRunning fitting-prep for SCH species: {[s.value for s in sch_species]}")
            success = success and run_fitting_prep_sch(sch_species, args.id)
        
        if success:
            print("\nFitting preparation completed successfully!")
        else:
            print("\nFitting preparation failed!", file=sys.stderr)
            sys.exit(1)
    
    elif stage == Stage.FITTING:
        success = True
        
        # Validate required arguments for fitting
        if args.id is None:
            print("Error: --id is required for fitting stage", file=sys.stderr)
            sys.exit(1)
        
        if sth_species:
            print(f"\nRunning fitting for STH species: {[s.value for s in sth_species]}")
            success = success and run_fitting_sth(sth_species, args)
        
        if sch_species:
            print(f"\nRunning fitting for SCH species: {[s.value for s in sch_species]}")
            success = success and run_fitting_sch(sch_species, args)
        
        if success:
            print("\nFitting completed successfully!")
        else:
            print("\nFitting failed!", file=sys.stderr)
            sys.exit(1)
    
    elif stage == Stage.PROJECTIONS_PREP:
        success = True
        
        if sth_species:
            print(f"\nRunning projections-prep for STH species: {[s.value for s in sth_species]}")
            success = success and run_projections_prep_sth(sth_species, args)
        
        if sch_species:
            print(f"\nRunning projections-prep for SCH species: {[s.value for s in sch_species]}")
            success = success and run_projections_prep_sch(sch_species, args)
        
        if success:
            print("\nProjections preparation completed successfully!")
        else:
            print("\nProjections preparation failed!", file=sys.stderr)
            sys.exit(1)
    
    elif stage == Stage.NEARTERM_PROJECTIONS:
        success = True
        
        # Validate required arguments for nearterm-projections
        if args.id is None:
            print("Error: --id is required for nearterm-projections stage", file=sys.stderr)
            sys.exit(1)
        
        if sth_species:
            print(f"\nRunning nearterm-projections for STH species: {[s.value for s in sth_species]}")
            success = success and run_nearterm_projections_sth(sth_species, args)
        
        if sch_species:
            print(f"\nRunning nearterm-projections for SCH species: {[s.value for s in sch_species]}")
            success = success and run_nearterm_projections_sch(sch_species, args)
        
        if success:
            print("\nNear-term projections completed successfully!")
        else:
            print("\nNear-term projections failed!", file=sys.stderr)
            sys.exit(1)
    
    elif stage == Stage.ALL:
        print("Running complete pipeline: fitting-prep → fitting → projections-prep → nearterm-projections")
        success = True
        
        # Stage 1: Fitting-prep
        print(f"\n{'='*60}")
        print("STAGE 1: FITTING-PREP")
        print(f"{'='*60}")
        
        if sth_species:
            print(f"Running fitting-prep for STH species: {[s.value for s in sth_species]}")
            success = success and run_fitting_prep_sth(sth_species, args.id)
        
        if sch_species:
            print(f"Running fitting-prep for SCH species: {[s.value for s in sch_species]}")
            success = success and run_fitting_prep_sch(sch_species, args.id)
        
        if not success:
            print("\nFitting preparation failed!", file=sys.stderr)
            sys.exit(1)
        
        # Stage 2: Fitting
        print(f"\n{'='*60}")
        print("STAGE 2: FITTING")
        print(f"{'='*60}")
        
        # Validate required arguments for fitting
        if args.id is None:
            print("Error: --id is required for fitting stage", file=sys.stderr)
            sys.exit(1)
        
        if sth_species:
            print(f"Running fitting for STH species: {[s.value for s in sth_species]}")
            success = success and run_fitting_sth(sth_species, args)
        
        if sch_species:
            print(f"Running fitting for SCH species: {[s.value for s in sch_species]}")
            success = success and run_fitting_sch(sch_species, args)
        
        if not success:
            print("\nFitting failed!", file=sys.stderr)
            sys.exit(1)
        
        # Stage 3: Projections-prep
        print(f"\n{'='*60}")
        print("STAGE 3: PROJECTIONS-PREP")
        print(f"{'='*60}")
        
        if sth_species:
            print(f"Running projections-prep for STH species: {[s.value for s in sth_species]}")
            success = success and run_projections_prep_sth(sth_species, args)
        
        if sch_species:
            print(f"Running projections-prep for SCH species: {[s.value for s in sch_species]}")
            success = success and run_projections_prep_sch(sch_species, args)
        
        if not success:
            print("\nProjections preparation failed!", file=sys.stderr)
            sys.exit(1)
        
        # Stage 4: Near-term projections
        print(f"\n{'='*60}")
        print("STAGE 4: NEAR-TERM PROJECTIONS")
        print(f"{'='*60}")
        
        # Validate required arguments for nearterm-projections
        if args.id is None:
            print("Error: --id is required for nearterm-projections stage", file=sys.stderr)
            sys.exit(1)
        
        if sth_species:
            print(f"Running nearterm-projections for STH species: {[s.value for s in sth_species]}")
            success = success and run_nearterm_projections_sth(sth_species, args)
        
        if sch_species:
            print(f"Running nearterm-projections for SCH species: {[s.value for s in sch_species]}")
            success = success and run_nearterm_projections_sch(sch_species, args)
        
        if not success:
            print("\nNear-term projections failed!", file=sys.stderr)
            sys.exit(1)
        
        print("✅ Complete pipeline executed successfully!")
        
        print(f"\n{'='*60}")
        print("COMPLETE PIPELINE SUMMARY")
        print(f"{'='*60}")
        print("✅ fitting-prep: COMPLETED")
        print("✅ fitting: COMPLETED") 
        print("✅ projections-prep: COMPLETED")
        print("✅ nearterm-projections: COMPLETED")
        print(f"{'='*60}")
    
    elif stage == Stage.SKIP_FITTING_PREP:
        print("Skip fitting-prep workflow not yet implemented")
        sys.exit(1)


if __name__ == "__main__":
    main()