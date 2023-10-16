from dataclasses import dataclass, field
from pathlib import Path

from ase.calculators.calculator import Calculator

from conformational_sampling import utils

@dataclass
class Config:
    xtb_path: str = '/export/apps/CentOS7/xtb/xtb/bin/xtb'
    initial_conformers: int = 100
    # set if you want to pass multiple ancillary ligand conformers and generate more
    # conformers based on each of them
    combinatorial_ancillary_ligand_confs = False
    # initial_rms_threshold: float = 0.6 # NOT NEEDED IN OPENBABEL IMPLEMENTATION
    max_connectivity_changes: int = 2
    pre_xtb_rms_threshold: float = 2.0
    # max number of BFGS geometry optimization steps; low default for debugging speed
    max_dft_opt_steps: int = 2
    # number of cpus to use for each dft geometry optimization
    # the number of simultaneous conformers that can be optimized is num_cpus//dft_cpus_per_opt
    dft_cpus_per_opt: int = 1
    num_cpus: int = field(default_factory=utils.num_cpus)
    ase_calculator: Calculator = None
    restart_gsm: Path = None
