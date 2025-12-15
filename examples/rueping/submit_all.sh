#!/bin/bash
# Helper script to submit both jobs with dependency

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
echo "Job status:"
echo "  Conformers: $JOB1"
echo "  GSM:        $JOB2 (dependency: afterok:$JOB1)"
echo ""
echo "Monitor with: squeue -u $USER"
echo "Cancel with:  scancel $JOB1 $JOB2"
