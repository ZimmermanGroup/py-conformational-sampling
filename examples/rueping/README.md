# Rueping Catalyst Conformational Sampling

This workflow generates conformers for the Rueping catalyst system and computes reaction pathways for the cationic 2-aza-Cope rearrangement.

## System Description

- **Molecular formula:** C67H56NO4P (129 atoms)
- **Catalyst + substrate complex** with critical hydrogen bond between catalyst O-H and substrate N
- **Reaction:** Cationic 2-aza-Cope [3,3]-sigmatropic rearrangement

## Boron Bridge Strategy

Due to OpenBabel not maintaining hydrogen bonds during conformer generation, we use a **boron substitution strategy**:

1. Replace bridging H with B in `full_system.xyz` (B bonds to both O and N)
2. Generate conformers with B maintaining the bridge connectivity
3. Replace B → H and reperceive bonds via XYZ format
4. Optimize directly with xTB (skip force field - it damages H-bond geometry)

## Configuration Parameters

### Conformer Generation (`rueping.py`)

```python
initial_conformers=100     # Number of initial conformers to generate
num_cpus=64               # CPUs for parallel optimization
pre_xtb_rms_threshold=0.175  # RMSD threshold (Å) for duplicate filtering
```

**RMSD Threshold Recommendations:**
- **0.175 Å** (current): Only removes very similar structures → keeps MORE unique conformers
- **0.3-0.5 Å**: Moderate filtering, removes somewhat similar structures  
- **1.0-2.0 Å** (default in Config): Aggressive filtering → keeps FEWER conformers

**Smaller threshold = more conformers kept** (only very close duplicates removed)  
**Larger threshold = fewer conformers kept** (more aggressive deduplication)

The current `pre_xtb_rms_threshold=0.175` Å is conservative and will keep more conformational diversity, which is good for exploring the full reaction landscape with 100 initial conformers.

### GSM Calculations (`rueping_gsm.py`)

```python
num_cpus=64  # CPUs for parallel GSM runs (one per conformer)
```

**Driving Coordinates** (1-indexed):
- `ADD: C(100)-C(86)` - Form new C-C bond
- `BREAK: C(93)-C(94)` - Break C-C bond

## Running on HPC (SLURM)

### Quick Start

**Recommended: Use the automated submission script**
```bash
# This script handles pixi setup and job submission automatically
./submit_all.sh
```

The script will:
1. Set up the pixi environment on the head node (downloads dependencies)
2. Verify Python and xTB are available
3. Submit the conformer generation job
4. Submit the GSM job with dependency (runs after conformers complete)

**Alternative: Submit jobs manually**
```bash
# First ensure pixi environment is set up (run on head node)
cd /path/to/py-conformational-sampling
pixi install

# Then submit jobs
sbatch run_conformers.slurm
sbatch run_gsm.slurm  # Run after first job completes
```

### SLURM Configuration

The provided SLURM scripts are configured for the **zimA10 partition** with:
- 1 node, 64 CPUs per node (full node)
- Memory: All available (512 GB)
- Time limits: 12h (conformers), 24h (GSM)
- Uses `pixi run python` for environment isolation

**Files:**
- `submit_all.sh` - **Automated setup and submission** (recommended)
- `run_conformers.slurm` - Conformer generation job
- `run_gsm.slurm` - GSM pathway calculations job

### Important Notes

- **Compute nodes lack internet access:** The pixi environment must be set up on the head node before jobs run (handled automatically by `submit_all.sh`)
- **Full node allocation:** Scripts request all 64 CPUs, which provides access to all 512 GB memory
- **Automatic dependency:** GSM job waits for conformer job to complete successfully

### Expected Outputs

**Step 1 (Conformers):**
- `conformers_0_with_boron.xyz` - Initial conformers with B bridge
- `conformers_0_unoptimized.xyz` - After B→H replacement
- `conformers_1_xtb.xyz` - xTB optimized, sorted by energy (~50-100 conformers)
- Runtime: ~2-8 hours on 64 CPUs

**Step 2 (GSM):**
- `scratch/pystring_*/opt_converged_001.xyz` - Optimized reaction pathways
- Energy profiles for each conformer
- Runtime: ~4-12 hours on 64 CPUs (depends on conformer count)

**Step 3 (Analysis - automatic in run_gsm.slurm):**
- Barrier heights for each conformer
- Reaction energies  
- Full energy profile comparison
- Identification of lowest barrier pathway

### Monitor Jobs

```bash
# Check job status
squeue -u $USER

# View output (while running)
tail -f conformers_*.out
tail -f gsm_*.out

# Cancel jobs if needed
scancel <job_id>
```

## Results from Test Run (3 conformers)

From our initial test with 3 conformers:

| Conformer | Barrier (kcal/mol) | Reaction ΔE (kcal/mol) |
|-----------|-------------------|----------------------|
| 0         | **0.67**          | -0.89                |
| 1         | 1.54              | -1.42                |
| 2         | 2.37              | -0.99                |

All pathways show very low barriers (< 2.4 kcal/mol), indicating the reaction proceeds readily at room temperature.

## Key Findings

1. **Hydrogen bond preservation:** xTB maintains the H-bond geometry (166° → 165.9°)
2. **Force field harmful:** MMFF/UFF damage H-bond (166° → 143°), so we skip this step
3. **Low barriers:** Cationic 2-aza-Cope has very accessible barriers across conformers
4. **Exothermic:** All pathways are thermodynamically favorable

## Notes for HPC Users

- No AI model required - all scripts run standalone
- Uses xTB via ASE for quantum calculations
- Parallelization via Python multiprocessing
- Results stored in XYZ format for easy visualization
- Adjust `num_cpus` to match your allocation
