#!/bin/bash
#SBATCH --job-name=rueping_gsm
#SBATCH --partition=zimA10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --mem=0                    # Request all memory on node
#SBATCH --output=gsm_%j.out
#SBATCH --error=gsm_%j.err
# #SBATCH --account=YOUR_ACCOUNT   # Uncomment and set if required

# Print job information
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"
echo ""

# Load required modules (adjust as needed for your cluster)
# module load python/3.11
# module load xtb
# Or activate conda/pixi environment:
# source /path/to/conda/etc/profile.d/conda.sh
# conda activate pycosa
# Or for pixi:
# eval "$(pixi shell-hook)"

# Ensure we're in the correct directory
cd $SLURM_SUBMIT_DIR

# Check that conformers_1_xtb.xyz exists
if [ ! -f conformers_1_xtb.xyz ]; then
    echo "ERROR: conformers_1_xtb.xyz not found!"
    echo "Please run run_conformers.slurm first"
    exit 1
fi

# Count number of conformers
NUM_CONFORMERS=$(grep -c "^129$" conformers_1_xtb.xyz)
echo "Found $NUM_CONFORMERS conformers from previous step"
echo "Running GSM calculations on 64 CPUs (parallel across conformers)"
echo "Expected runtime: 4-12 hours for ~50-100 conformers"
echo ""

# Run GSM calculations
python rueping_gsm.py

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "GSM calculations completed successfully at: $(date)"
    echo ""
    echo "Analyzing results..."
    python analyze_gsm_results.py
    echo ""
    echo "Results are in scratch/pystring_*/ directories"
    echo "See analysis output above for barrier heights and reaction energies"
else
    echo ""
    echo "ERROR: GSM calculations failed at: $(date)"
    exit 1
fi
