from dataclasses import dataclass, field
from conformational_sampling import utils

@dataclass
class Config:
    num_cpus: int = field(default_factory=utils.num_cpus)
    initial_conformers: int = 100
    initial_rms_threshold: float = 0.6
    max_connectivity_changes = 2
    pre_xtb_rms_threshold: float = 2.0