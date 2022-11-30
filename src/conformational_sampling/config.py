from dataclasses import dataclass, field

from ase.calculators.calculator import Calculator
from xtb.ase.calculator import XTB as XTBCalc

from conformational_sampling import utils

@dataclass
class Config:
    xtb_path: str
    initial_conformers: int = 100
    # initial_rms_threshold: float = 0.6 # NOT NEEDED IN OPENBABEL IMPLEMENTATION
    max_connectivity_changes: int = 2
    pre_xtb_rms_threshold: float = 2.0
    num_cpus: int = field(default_factory=utils.num_cpus)
    ase_calculator: Calculator = XTBCalc()
