from dataclasses import dataclass, field
from pathlib import Path

from ase.calculators.calculator import Calculator

from conformational_sampling import utils

@dataclass
class Config:
    xtb_path: str = 'xtb'
    initial_conformers: int = 100
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
