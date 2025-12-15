#!/bin/bash
#SBATCH --job-name=rueping_conformers
#SBATCH --partition=zimA10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --mem=0                    # Request all memory on node
#SBATCH --output=conformers_%j.out
#SBATCH --error=conformers_%j.err

# Print job information
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"
echo ""

# Ensure we're in the correct directory
cd $SLURM_SUBMIT_DIR

# Run conformer generation using pixi
# Note: pixi environment must be set up on head node before submitting
# (compute nodes don't have internet access for downloading dependencies)
echo "Starting conformer generation with 100 initial conformers on 64 CPUs"
echo "Expected runtime: 2-8 hours depending on system"
echo ""

pixi run python rueping.py

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "Conformer generation completed successfully at: $(date)"
    echo "Output files:"
    ls -lh conformers_*.xyz
    echo ""
    echo "Number of conformers generated:"
    grep -c "^129$" conformers_1_xtb.xyz
    echo ""
    echo "Next step: Submit run_gsm.slurm to compute reaction pathways"
else
    echo ""
    echo "ERROR: Conformer generation failed at: $(date)"
    exit 1
fi
