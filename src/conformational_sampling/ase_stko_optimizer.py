from dataclasses import dataclass

import ase
import stko
from ase.calculators.calculator import Calculator
from ase.optimize import BFGS

from conformational_sampling.utils import stk_mol_to_ase_atoms


@dataclass
class ASE(stko.optimizers.Optimizer):
    calculator: Calculator

    def optimize(self, stk_mol):
        ase_mol = stk_mol_to_ase_atoms(stk_mol)
        ase_mol.calc = self.calculator
        opt = BFGS(ase_mol)
        opt.run(fmax=0.2)
        stk_mol.with_position_matrix(opt.atoms.get_positions())
        return stk_mol