import os

def num_cpus():
    try:
        return int(os.environ["SLURM_CPUS_PER_TASK"])
    except KeyError:
        return 2
