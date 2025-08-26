# STH/SCH AMIS Pipeline

## 🚀 Quick Start

### Build Docker Image
```bash
DOCKER_BUILDKIT=1 docker build --ssh default=$SSH_AUTH_SOCK . -t sch-sth-amis-pipeline
```

### Scenario 1: Skip Fitting-Prep (Use Pre-built Artifacts)
**Most Common**: Use embedded fitting-prep artifacts, run fitting→projections-prep→nearterm-projections

```bash
# STH - Complete pipeline from fitting onwards
docker run --rm \
  -v $(pwd)/fitting/artefacts:/ntdmc/sth-sch-amis-integration/fitting/artefacts \
  -v $(pwd)/projections-prep/artefacts:/ntdmc/sth-sch-amis-integration/projections-prep/artefacts \
  -v $(pwd)/projections/artefacts:/ntdmc/sth-sch-amis-integration/projections/artefacts \
  sch-sth-amis-pipeline \
  --stage=skip-fitting-prep --disease=sth --sth-species=ascaris --id=10 --num-cores=8

# SCH - Complete pipeline from fitting onwards  
docker run --rm \
  -v $(pwd)/fitting/artefacts:/ntdmc/sth-sch-amis-integration/fitting/artefacts \
  -v $(pwd)/projections-prep/artefacts:/ntdmc/sth-sch-amis-integration/projections-prep/artefacts \
  -v $(pwd)/projections/artefacts:/ntdmc/sth-sch-amis-integration/projections/artefacts \
  sch-sth-amis-pipeline \
  --stage=skip-fitting-prep --disease=sch --sch-species=haematobium --id=10 --num-cores=8
```

### Scenario 2: Fitting Only (Distributed Batch Processing)
**Cloud Execution**: Run only fitting stage with pre-built artifacts

```bash
# STH Fitting - Production parameters
docker run --rm \
  -v $(pwd)/fitting/artefacts:/ntdmc/sth-sch-amis-integration/fitting/artefacts \
  sch-sth-amis-pipeline \
  --stage=fitting --disease=sth --sth-species=ascaris --id=10 \
  --amis-n-samples=1000 --amis-n-iters=50 --num-cores=8

# SCH Fitting - Production parameters  
docker run --rm \
  -v $(pwd)/fitting/artefacts:/ntdmc/sth-sch-amis-integration/fitting/artefacts \
  sch-sth-amis-pipeline \
  --stage=fitting --disease=sch --sch-species=haematobium --id=10 \
  --amis-n-samples=500 --amis-n-iters=50 --num-cores=8
```

### Scenario 3: Nearterm-Projections Only
**Post-Processing**: Run projections with existing fitting results

```bash
# STH Projections
docker run --rm \
  -v $(pwd)/projections-prep/artefacts:/ntdmc/sth-sch-amis-integration/projections-prep/artefacts \
  -v $(pwd)/projections/artefacts:/ntdmc/sth-sch-amis-integration/projections/artefacts \
  sch-sth-amis-pipeline \
  --stage=nearterm-projections --disease=sth --sth-species=ascaris --id=10 --num-cores=8

# SCH Projections
docker run --rm \
  -v $(pwd)/projections-prep/artefacts:/ntdmc/sth-sch-amis-integration/projections-prep/artefacts \
  -v $(pwd)/projections/artefacts:/ntdmc/sth-sch-amis-integration/projections/artefacts \
  sch-sth-amis-pipeline \
  --stage=nearterm-projections --disease=sch --sch-species=haematobium --id=10 --num-cores=8
```

## 🚀 Complete Stage Reference

### Stage 1: Fitting-Prep (Run Once Per Species) 
*Only needed if creating custom artifacts*

```bash
# STH - creates data for all batches
docker run --rm \
  -v $(pwd)/fitting-prep/artefacts:/ntdmc/sth-sch-amis-integration/fitting-prep/artefacts \
  sch-sth-amis-pipeline \
  --stage=fitting-prep --disease=sth

# SCH - creates data for all batches  
docker run --rm \
  -v $(pwd)/fitting-prep/artefacts:/ntdmc/sth-sch-amis-integration/fitting-prep/artefacts \
  sch-sth-amis-pipeline \
  --stage=fitting-prep --disease=sch --sch-species=haematobium
```

### Stage 2: Fitting (Run Per Batch)
```bash
# STH
docker run --rm \
  -v $(pwd)/fitting-prep/artefacts:/ntdmc/sth-sch-amis-integration/fitting-prep/artefacts \
  -v $(pwd)/fitting/artefacts:/ntdmc/sth-sch-amis-integration/fitting/artefacts \
  sch-sth-amis-pipeline \
  --stage=fitting --disease=sth --sth-species=ascaris --id=10 \
  --amis-n-samples=1000 --amis-n-iters=50 --num-cores=8

# SCH
docker run --rm \
  -v $(pwd)/fitting-prep/artefacts:/ntdmc/sth-sch-amis-integration/fitting-prep/artefacts \
  -v $(pwd)/fitting/artefacts:/ntdmc/sth-sch-amis-integration/fitting/artefacts \
  sch-sth-amis-pipeline \
  --stage=fitting --disease=sch --sch-species=haematobium --id=10 \
  --amis-n-samples=500 --amis-n-iters=50 --num-cores=8
```

### Stage 3: Projections-Prep (Run Per Batch)
```bash
# STH
docker run --rm \
  -v $(pwd)/fitting-prep/artefacts:/ntdmc/sth-sch-amis-integration/fitting-prep/artefacts \
  -v $(pwd)/fitting/artefacts:/ntdmc/sth-sch-amis-integration/fitting/artefacts \
  -v $(pwd)/projections-prep/artefacts:/ntdmc/sth-sch-amis-integration/projections-prep/artefacts \
  sch-sth-amis-pipeline \
  --stage=projections-prep --disease=sth --sth-species=ascaris --id=10

# SCH  
docker run --rm \
  -v $(pwd)/fitting-prep/artefacts:/ntdmc/sth-sch-amis-integration/fitting-prep/artefacts \
  -v $(pwd)/fitting/artefacts:/ntdmc/sth-sch-amis-integration/fitting/artefacts \
  -v $(pwd)/projections-prep/artefacts:/ntdmc/sth-sch-amis-integration/projections-prep/artefacts \
  sch-sth-amis-pipeline \
  --stage=projections-prep --disease=sch --sch-species=haematobium --id=10
```

### Stage 4: Near-term Projections (Run Per Batch)
```bash
# STH
docker run --rm \
  -v $(pwd)/fitting-prep/artefacts:/ntdmc/sth-sch-amis-integration/fitting-prep/artefacts \
  -v $(pwd)/projections-prep/artefacts:/ntdmc/sth-sch-amis-integration/projections-prep/artefacts \
  -v $(pwd)/projections/artefacts:/ntdmc/sth-sch-amis-integration/projections/artefacts \
  sch-sth-amis-pipeline \
  --stage=nearterm-projections --disease=sth --sth-species=ascaris --id=10 --num-cores=8

# SCH
docker run --rm \
  -v $(pwd)/fitting-prep/artefacts:/ntdmc/sth-sch-amis-integration/fitting-prep/artefacts \
  -v $(pwd)/projections-prep/artefacts:/ntdmc/sth-sch-amis-integration/projections-prep/artefacts \
  -v $(pwd)/projections/artefacts:/ntdmc/sth-sch-amis-integration/projections/artefacts \
  sch-sth-amis-pipeline \
  --stage=nearterm-projections --disease=sch --sch-species=haematobium --id=10 --num-cores=8
```

## 📁 Expected Output Files

### After Fitting Stage
```
fitting/artefacts/AMIS_output/{species}_amis_output{id}.Rdata
```
**Location**: `fitting/artefacts/AMIS_output/`
**Examples**: `ascaris_amis_output10.Rdata`, `haematobium_amis_output25.Rdata`

### After Projections-Prep Stage  
```
projections-prep/artefacts/InputPars_MTP_{species}/InputPars_MTP_{iu}.csv
```
**Location**: `projections-prep/artefacts/InputPars_MTP_{species}/`
**Examples**: `InputPars_MTP_ascaris/InputPars_MTP_12345.csv`

### After Near-term Projections Stage
```
projections/artefacts/projections/{species}/{country}/{country}{iu}/{Species}_{country}{iu}.p
```
**Location**: `projections/artefacts/projections/{species}/{country}/{country}{iu}/`
**Examples**: `projections/ascaris/TZA/TZA12345/Asc_TZA12345.p`

## 🔧 Key Parameters

| Parameter | STH Production | SCH Production | Purpose |
|-----------|---------------|----------------|---------|
| `--id` | `10,11,12...` | `10,11,12...` | Batch ID |
| `--amis-n-samples` | `1000` | `500` | AMIS samples per iteration |
| `--amis-n-iters` | `50` | `50` | Maximum AMIS iterations |
| `--num-cores` | `8-16` | `8-16` | CPU cores (match instance) |
| `--ess-threshold` | `200` | `200` | ESS threshold for file selection |

## 🔄 Two-Stage Fitting Workflow

**Important**: If fitting shows low ESS, you need a second run with higher sigma.

### 1. Check Fitting Results
After fitting, check the manifest:
```
📋 FITTING MANIFEST - Batch 10
🔬 ASCARIS
   📊 ESS: 1 IUs, min: 4.8 max: 4.8 mean: 4.8 below_threshold: 1
   ⚠️  1 IUs with ESS < 200 - consider --amis-sigma=0.025
```

### 2. Rerun Low-ESS Batches
```bash
# Rerun with higher sigma for better convergence
docker run --rm \
  -v $(pwd)/fitting-prep/artefacts:/ntdmc/sth-sch-amis-integration/fitting-prep/artefacts \
  -v $(pwd)/fitting/artefacts:/ntdmc/sth-sch-amis-integration/fitting/artefacts \
  sch-sth-amis-pipeline \
  --stage=fitting --disease=sth --sth-species=ascaris --id=10 --amis-sigma=0.025
```

This creates: `ascaris_amis_output10_sigma0.025.Rdata`

**Note**: Projections-prep will automatically use the sigma file for low-ESS IUs.

## 🚨 Common Issues

### 1. Projections-prep fails: "cannot open compressed file"
**Error**: `Error: cannot open compressed file 'ascaris_amis_output10_sigma0.0025.Rdata'`
**Fix**: Rerun fitting with `--amis-sigma=0.025` OR use `--ess-threshold=1` for testing

### 2. Missing files
**Error**: Various file not found errors
**Fix**: Ensure all volume mounts are specified and fitting-prep was run

## 🗂️ Supported Species

### STH
- `ascaris`
- `hookworm`
- `trichuris`

### SCH
- `haematobium`
- `mansoni_low_burden`
- `mansoni_high_burden`

---

# 📖 Detailed Reference

## Pipeline Overview

The pipeline consists of four sequential stages:
1. **fitting-prep**: Prepares data for fitting (run once per species)
2. **fitting**: Runs AMIS fitting algorithm (run per batch) 
3. **projections-prep**: Prepares data for projections (run per batch)
4. **nearterm-projections**: Runs projection simulations (run per batch)

## 📊 Pipeline Data Flow Diagram

```mermaid
graph TB
    %% Input Data Sources
    subgraph "Input Data (Google Cloud Storage)"
        MAPS_STH[("📁 Maps-STH<br/>• STHCleaned_1.csv<br/>• sartorious_2021.csv<br/>• LF_MDA_Africa_2024_IU_updated.csv")]
        MAPS_SCH[("📁 Maps-SCH<br/>• Schisto_IU_Cleaned_1.csv<br/>• S.haematobium_G*.xlsx<br/>• G*.mansoni_*.csv")]
        ESPEN[("📁 ESPEN_IU_2021<br/>• ESPEN_IU_2021.shp<br/>• ESPEN_IU_2021.dbf")]
    end

    %% Stage 1: Fitting-Prep
    subgraph "Stage 1: Fitting-Prep"
        FP_STH["🔧 STH Scripts<br/>• prepare_histories.R<br/>• prepare_histories_trichuris.R<br/>• prepare_maps_allspecies.R<br/>• prepare_histories_projections.R"]
        FP_SCH["🔧 SCH Scripts<br/>• prepare_histories_and_maps_haematobium.R<br/>• prepare_histories_and_maps_mansoni.R<br/>• prepare_histories_projections_sch.R"]
        
        FP_OUT[("📤 Artefacts<br/>• endgame_inputs/{species}/InputMDA_MTP_{id}.xlsx<br/>• endgame_inputs/{species}/InputMDA_MTP_projections_{iu}.xlsx<br/>• Maps/table_iu_idx_{species}.csv<br/>• Maps/iu_task_lookup_{species}.rds<br/>• Maps/proj_iu_task_lookup_{species}.rds")]
    end

    %% Stage 2: Fitting
    subgraph "Stage 2: Fitting"
        FIT_STH["🔧 STH Fitting<br/>• sth_fitting.R<br/>• amis_integration.R<br/>• {species}_prior.R"]
        FIT_SCH["🔧 SCH Fitting<br/>• sch_fitting.R<br/>• amis_integration.R<br/>• sch_prior.R"]
        
        FIT_OUT[("📤 Artefacts<br/>• AMIS_output/{species}_amis_output{id}.Rdata<br/>• AMIS_output/{species}_amis_output{id}_sigma{value}.Rdata<br/>• fitting_manifest_batch_{id}.json")]
    end

    %% Stage 3: Projections-Prep
    subgraph "Stage 3: Projections-Prep"
        PP_SCRIPT["🔧 Scripts<br/>• preprocess_for_projections.R<br/>• IUsWithInsufficientESS.R"]
        
        PP_OUT[("📤 Artefacts<br/>• InputPars_MTP_{species}/InputPars_MTP_{iu}.csv<br/>• post_AMIS_analysis/proc_output_{species}_{id}.csv")]
    end

    %% Stage 4: Nearterm-Projections
    subgraph "Stage 4: Nearterm-Projections"
        NP_STH["🔧 STH Projections<br/>• sth_projections_per_IU.py"]
        NP_SCH["🔧 SCH Projections<br/>• sch_projections_per_IU.py"]
        
        NP_OUT[("📤 Artefacts<br/>• projections/{species}/{country}/{country}{iu}/<br/>  - {Species}_{country}{iu}.p<br/>  - PrevDataset_{Species}_{country}{iu}.csv")]
    end

    %% Data Flow Connections
    MAPS_STH --> FP_STH
    MAPS_SCH --> FP_SCH
    ESPEN --> FP_SCH
    
    FP_STH --> FP_OUT
    FP_SCH --> FP_OUT
    
    FP_OUT --> FIT_STH
    FP_OUT --> FIT_SCH
    
    FIT_STH --> FIT_OUT
    FIT_SCH --> FIT_OUT
    
    FIT_OUT --> PP_SCRIPT
    FP_OUT -.->|"lookup tables"| PP_SCRIPT
    
    PP_SCRIPT --> PP_OUT
    
    PP_OUT --> NP_STH
    PP_OUT --> NP_SCH
    FP_OUT -.->|"coverage files"| NP_STH
    FP_OUT -.->|"coverage files"| NP_SCH
    
    NP_STH --> NP_OUT
    NP_SCH --> NP_OUT

    %% Styling
    classDef inputData fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    classDef script fill:#fff3e0,stroke:#e65100,stroke-width:2px
    classDef artefact fill:#f3e5f5,stroke:#4a148c,stroke-width:2px
    
    class MAPS_STH,MAPS_SCH,ESPEN inputData
    class FP_STH,FP_SCH,FIT_STH,FIT_SCH,PP_SCRIPT,NP_STH,NP_SCH script
    class FP_OUT,FIT_OUT,PP_OUT,NP_OUT artefact
```

## 📁 Detailed Input/Output Specifications

### Stage 1: Fitting-Prep

**Inputs:**
- **STH**: 
  - `Maps-STH/STHCleaned_1.csv` - STH prevalence data
  - `Maps-STH/sartorious_2021.csv` - Disability weights
  - `Maps-STH/LF_MDA_Africa_2024_IU_updated.csv` - MDA coverage
- **SCH**:
  - `Maps-SCH/Schisto_IU_Cleaned_1.csv` - SCH prevalence data
  - `Maps-SCH/S.haematobium_G*.xlsx` - Haematobium group data (G1-G6, 2005/2013/2023)
  - `Maps-SCH/G*.mansoni_*.csv` - Mansoni group data (G1-G6, 2005/2013/2023)
  - `ESPEN_IU_2021/*` - ESPEN shapefiles for IU boundaries

**Outputs:**
```
fitting-prep/artefacts/
├── endgame_inputs/
│   ├── STH/
│   │   ├── InputMDA_MTP_{id}.xlsx                    # Coverage files for fitting (per batch)
│   │   └── InputMDA_MTP_projections_{iu}.xlsx        # Coverage files for projections (per IU)
│   ├── sch-haematobium/
│   │   ├── InputMDA_MTP_{id}.xlsx
│   │   └── InputMDA_MTP_projections_{iu}.xlsx
│   └── sch-mansoni/
│       ├── InputMDA_MTP_{id}.xlsx
│       └── InputMDA_MTP_projections_{iu}.xlsx
└── Maps/
    ├── table_iu_idx_{species}.csv                    # IU-TaskID mapping table
    ├── iu_task_lookup_{species}.rds                  # Fitting batch lookup
    └── proj_iu_task_lookup_{species}.rds             # Projection batch lookup
```

### Stage 2: Fitting

**Inputs:**
- `fitting-prep/artefacts/endgame_inputs/{species}/InputMDA_MTP_{id}.xlsx`
- `fitting-prep/artefacts/Maps/iu_task_lookup_{species}.rds`
- Model priors: `{species}_prior.R`

**Outputs:**
```
fitting/artefacts/
├── AMIS_output/
│   ├── {species}_amis_output{id}.Rdata              # Default sigma (0.0025)
│   └── {species}_amis_output{id}_sigma{value}.Rdata # Custom sigma (e.g., 0.025)
└── fitting_manifest_batch_{id}.json                 # ESS analysis and metadata
```

### Stage 3: Projections-Prep

**Inputs:**
- `fitting/artefacts/AMIS_output/{species}_amis_output{id}*.Rdata`
- `fitting-prep/artefacts/Maps/table_iu_idx_{species}.csv`
- `fitting-prep/artefacts/Maps/proj_iu_task_lookup_{species}.rds`

**Outputs:**
```
projections-prep/artefacts/
├── InputPars_MTP_{species}/
│   └── InputPars_MTP_{iu}.csv                       # Sampled parameters per IU
└── post_AMIS_analysis/
    └── proc_output_{species}_{id}.csv               # Processed AMIS output
```

### Stage 4: Nearterm-Projections

**Inputs:**
- `projections-prep/artefacts/InputPars_MTP_{species}/InputPars_MTP_{iu}.csv`
- `fitting-prep/artefacts/endgame_inputs/{species}/InputMDA_MTP_projections_{iu}.xlsx`
- `fitting-prep/artefacts/Maps/table_iu_idx_{species}.csv`
- Model parameters: `{species}_params_projections.txt`

**Outputs:**
```
projections/artefacts/
└── projections/
    └── {species}/
        └── {country}/
            └── {country}{iu}/
                ├── {Species}_{country}{iu}.p        # Projection results (pickle)
                └── PrevDataset_{Species}_{country}{iu}.csv  # Prevalence dataset
```

## 🔀 Cross-Stage Dependencies

The pipeline has both linear and cross-stage dependencies:

```mermaid
graph LR
    subgraph "Direct Dependencies"
        FP[fitting-prep] --> F[fitting]
        F --> PP[projections-prep]
        PP --> NP[nearterm-projections]
    end
    
    subgraph "Cross-Stage Dependencies"
        FP2[fitting-prep] -.->|"coverage files<br/>lookup tables"| NP2[nearterm-projections]
    end
    
    style FP fill:#e1f5fe
    style F fill:#fff3e0
    style PP fill:#f3e5f5
    style NP fill:#e8f5e9
    style FP2 fill:#e1f5fe
    style NP2 fill:#e8f5e9
```

**Key Cross-Dependencies:**
- `nearterm-projections` requires files from **both** `fitting-prep` (coverage files, lookup tables) and `projections-prep` (parameter files)
- `projections-prep` uses lookup tables from `fitting-prep` to map TaskIDs to IUs

## Volume Mount Requirements

| Mount Path | Purpose | Required For |
|------------|---------|--------------|
| `/fitting-prep/artefacts` | Input data and prepared files | All stages |
| `/fitting/artefacts` | AMIS fitting outputs | fitting, projections-prep, nearterm-projections |
| `/projections-prep/artefacts` | Projection prep outputs | projections-prep, nearterm-projections |
| `/projections/artefacts` | Final projection results | nearterm-projections |

## Complete Parameter Reference

### Required Parameters
- `--stage`: `fitting-prep`, `fitting`, `projections-prep`, `nearterm-projections`, `all`
- `--disease`: `sth` or `sch`
- `--id`: Batch ID (required for fitting, projections-prep, nearterm-projections)

### Species Selection
- `--sth-species`: `ascaris`, `hookworm`, `trichuris`, `all`
- `--sch-species`: `haematobium`, `mansoni_low_burden`, `mansoni_high_burden`, `all`

### AMIS Parameters (Fitting Stage)
- `--amis-sigma`: Sigma parameter (default: 0.0025, use 0.025 for reruns)
- `--amis-n-samples`: Samples per iteration (default: 1000 STH, 500 SCH)
- `--amis-target-ess`: Target ESS (default: 500)
- `--amis-n-iters`: Max iterations (default: 50)
- `--num-cores`: CPU cores (default: auto-detect)

### Projection Parameters
- `--ess-threshold`: ESS threshold (default: 200)
- `--failed-ids`: Comma-separated failed batch IDs

## Development/Testing Parameters

For quick testing (not production):
```bash
--amis-n-samples=10 --amis-n-iters=10 --amis-target-ess=1 --ess-threshold=1 --num-cores=2
```

## File Naming Conventions

### AMIS Output Files
- Default sigma: `{species}_amis_output{id}.Rdata`
- Custom sigma: `{species}_amis_output{id}_sigma{value}.Rdata`

### Species Directory Structure
- **STH**: All species use `endgame_inputs/STH/`
- **SCH haematobium**: Uses `endgame_inputs/sch-haematobium/`
- **SCH mansoni**: Uses `endgame_inputs/sch-mansoni/`

## Pipeline Dependencies

```mermaid
graph LR
    FP[fitting-prep] --> F[fitting]
    F --> PP[projections-prep]
    PP --> NP[nearterm-projections]
    FP --> NP
```

**Cross-stage dependencies**:
- `nearterm-projections` needs files from both `fitting-prep` and `projections-prep`
- `projections-prep` intelligently selects AMIS files based on ESS thresholds

## ESS Analysis Script

Identify batches needing sigma reruns:
```bash
docker run --rm \
  -v $(pwd)/fitting-prep/artefacts:/ntdmc/sth-sch-amis-integration/fitting-prep/artefacts \
  -v $(pwd)/fitting/artefacts:/ntdmc/sth-sch-amis-integration/fitting/artefacts \
  --entrypoint Rscript \
  sch-sth-amis-pipeline \
  fitting/scripts/find_lowESS_ids.R --species=ascaris --ess-threshold=200
```

## Debug Mode

For troubleshooting:
```bash
docker run --name debug-run \
  -v $(pwd)/fitting-prep/artefacts:/ntdmc/sth-sch-amis-integration/fitting-prep/artefacts \
  sch-sth-amis-pipeline \
  --stage=fitting-prep --disease=sth --id=10

# Check logs
docker logs debug-run

# Clean up  
docker rm debug-run
```

## 🔧 Technical Notes

### SCH Parameter Files Workaround

**Issue**: The `run-amis-fitting` branch of `ntd-model-sch` is missing the `*_params_projections.txt` files required by SCH nearterm-projections stage.

**Solution**: During Docker build, we temporarily clone the `updateImportation` branch and copy the missing SCH projection parameter files:
- `haematobium_params_projections.txt`
- `mansoni_low_burden_params_projections.txt` 
- `mansoni_high_burden_params_projections.txt`

**Location in Dockerfile**:
```dockerfile
# WORKAROUND: Copy missing SCH projection parameter files from updateImportation branch
RUN mkdir -p /tmp/updateImportation
ADD --keep-git-dir git@github.com:NTD-Modelling-Consortium/ntd-model-sch.git#updateImportation /tmp/updateImportation
RUN cp -f /tmp/updateImportation/sch_simulation/data/SCH_params/*_params_projections.txt \
    ${STH_SCH_MODEL_DIR}/sch_simulation/data/SCH_params/ || echo "No projection parameter files found in updateImportation"
RUN rm -rf /tmp/updateImportation
```

**Future**: This workaround can be removed once the missing parameter files are added to the `run-amis-fitting` branch.
