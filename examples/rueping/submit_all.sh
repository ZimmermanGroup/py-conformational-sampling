#!/bin/bash
# Helper script to submit both jobs with dependency
# This script should be run on the head node

set -e  # Exit on error

echo "========================================="
echo "Rueping Catalyst SLURM Job Submission"
echo "========================================="
echo ""

# Step 1: Setup pixi environment on head node
echo "Step 1: Setting up pixi environment (compute nodes lack internet access)"
echo "----------------------------------------"

# Check if pixi is available
if ! command -v pixi &> /dev/null; then
    echo "ERROR: pixi not found in PATH"
    echo "Please install pixi or ensure it's available on the head node"
    exit 1
fi

echo "Found pixi: $(which pixi)"
echo ""

# Navigate to repository root to run pixi install
REPO_ROOT=$(git rev-parse --show-toplevel 2>/dev/null || echo ".")
echo "Repository root: $REPO_ROOT"
cd "$REPO_ROOT"

echo "Running pixi install to download dependencies..."
pixi install

echo ""
echo "Verifying environment..."
pixi run python --version
pixi run xtb --version 2>/dev/null || echo "Warning: xTB not found in environment"

echo ""
echo "✓ Pixi environment ready"
echo ""

# Return to submission directory
cd - > /dev/null

# Step 2: Submit SLURM jobs
echo "Step 2: Submitting SLURM jobs"
echo "----------------------------------------"
echo ""

echo "Submitting conformer generation job..."
JOB1=$(sbatch --parsable run_conformers.slurm)

if [ -z "$JOB1" ]; then
    echo "ERROR: Failed to submit conformer generation job"
    exit 1
fi

echo "Conformer job submitted with ID: $JOB1"
echo ""
echo "Submitting GSM job to run after conformer job completes..."

JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 run_gsm.slurm)

if [ -z "$JOB2" ]; then
    echo "ERROR: Failed to submit GSM job"
    exit 1
fi

echo "GSM job submitted with ID: $JOB2 (will start after job $JOB1 completes)"
echo ""
echo "========================================="
echo "Submission Complete"
echo "========================================="
echo "Job status:"
echo "  Conformers: $JOB1 (run_conformers.slurm)"
echo "  GSM:        $JOB2 (run_gsm.slurm, dependency: afterok:$JOB1)"
echo ""
echo "Monitor with: squeue -u $USER"
echo "View logs:    tail -f conformers_$JOB1.out"
echo "              tail -f gsm_$JOB2.out"
echo "Cancel with:  scancel $JOB1 $JOB2"
echo "========================================="
