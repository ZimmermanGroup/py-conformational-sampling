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

## Running on HPC

### Step 1: Conformer Generation
```bash
# Generates conformers_1_xtb.xyz with ~100 unique conformers
python rueping.py
```

**Expected output files:**
- `conformers_0_with_boron.xyz` - Initial conformers with B bridge
- `conformers_0_unoptimized.xyz` - After B→H replacement
- `conformers_1_xtb.xyz` - xTB optimized, sorted by energy

**Runtime estimate:** ~2-8 hours on 64 CPUs (depends on system)

### Step 2: GSM Reaction Pathways
```bash
# Computes reaction pathways for all conformers
python rueping_gsm.py
```

**Expected output:**
- `scratch/pystring_*/opt_converged_001.xyz` - Optimized reaction pathways
- Energy profiles for each conformer

**Runtime estimate:** ~4-12 hours on 64 CPUs (depends on number of conformers)

### Step 3: Analyze Results
```bash
python analyze_gsm_results.py
```

**Output:**
- Barrier heights for each conformer
- Reaction energies
- Full energy profiles
- Identification of lowest barrier pathway

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
