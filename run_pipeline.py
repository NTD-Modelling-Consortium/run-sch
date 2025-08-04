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


def run_fitting_prep_sth(species_list, for_projections=False):
    """Run fitting preparation for STH species."""
    base_path = "/ntdmc/sth-sch-amis-integration" if os.path.exists('/.dockerenv') else "."
    scripts_path = f"{base_path}/fitting-prep/scripts"
    
    # First, prepare histories for fitting
    if not for_projections:
        commands = []
        
        # For STH maps preparation, we need BOTH ascaris/hookworm AND trichuris histories
        # because prepare_maps_allspecies.R reads both files
        needs_ascaris_hookworm = (STHSpecies.ASCARIS in species_list or 
                                 STHSpecies.HOOKWORM in species_list or 
                                 STHSpecies.ALL in species_list)
        needs_trichuris = (STHSpecies.TRICHURIS in species_list or 
                          STHSpecies.ALL in species_list)
        
        if needs_ascaris_hookworm:
            commands.append(("Rscript prepare_histories.R", "Preparing histories for ascaris/hookworm"))
        
        if needs_trichuris or needs_ascaris_hookworm:  # Always run trichuris if any STH species needed
            commands.append(("Rscript prepare_histories_trichuris.R", "Preparing histories for trichuris"))
        
        # Always run maps preparation after histories
        if commands:
            commands.append(("Rscript prepare_maps_allspecies.R", "Preparing maps for all STH species"))
        
        for cmd, desc in commands:
            if not run_command(cmd, desc, cwd=scripts_path):
                return False
    else:
        # Prepare histories for projections
        commands = []
        
        if STHSpecies.ASCARIS in species_list or STHSpecies.HOOKWORM in species_list or STHSpecies.ALL in species_list:
            commands.append(("Rscript prepare_histories_projections.R", "Preparing projection histories for ascaris/hookworm"))
        
        if STHSpecies.TRICHURIS in species_list or STHSpecies.ALL in species_list:
            commands.append(("Rscript prepare_histories_trichuris_projections.R", "Preparing projection histories for trichuris"))
        
        for cmd, desc in commands:
            if not run_command(cmd, desc, cwd=scripts_path):
                return False
    
    return True


def run_fitting_prep_sch(species_list, for_projections=False, batch_id=None):
    """Run fitting preparation for SCH species."""
    base_path = "/ntdmc/sth-sch-amis-integration" if os.path.exists('/.dockerenv') else "."
    scripts_path = f"{base_path}/fitting-prep/scripts"
    
    if not for_projections:
        commands = []
        
        batch_arg = f" --id {batch_id}" if batch_id else ""
        
        if SCHSpecies.HAEMATOBIUM in species_list or SCHSpecies.ALL in species_list:
            commands.append((f"Rscript prepare_histories_and_maps_haematobium.R{batch_arg}", "Preparing histories and maps for haematobium"))
        
        if SCHSpecies.MANSONI_HIGH in species_list or SCHSpecies.MANSONI_LOW in species_list or SCHSpecies.ALL in species_list:
            commands.append((f"Rscript prepare_histories_and_maps_mansoni.R{batch_arg}", "Preparing histories and maps for mansoni"))
        
        for cmd, desc in commands:
            if not run_command(cmd, desc, cwd=scripts_path):
                return False
    else:
        # Prepare histories for projections
        if not run_command("Rscript prepare_histories_projections_sch.R", 
                         "Preparing projection histories for SCH species", cwd=scripts_path):
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
    
    parser.add_argument(
        "--for-projections",
        action="store_true",
        help="Prepare data for projections (vs fitting)"
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
            success = success and run_fitting_prep_sth(sth_species, args.for_projections)
        
        if sch_species:
            print(f"\nRunning fitting-prep for SCH species: {[s.value for s in sch_species]}")
            success = success and run_fitting_prep_sch(sch_species, args.for_projections, args.id)
        
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
        print("Projections preparation stage not yet implemented")
        sys.exit(1)
    
    elif stage == Stage.NEARTERM_PROJECTIONS:
        print("Near-term projections stage not yet implemented")
        sys.exit(1)
    
    elif stage == Stage.ALL:
        print("Running all stages not yet implemented")
        sys.exit(1)
    
    elif stage == Stage.SKIP_FITTING_PREP:
        print("Skip fitting-prep workflow not yet implemented")
        sys.exit(1)


if __name__ == "__main__":
    main()