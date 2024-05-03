from dataclasses import dataclass
import os

import ase
import stk
import stko
from ase.calculators.calculator import CalculationFailed, Calculator
from ase.optimize import BFGS

from conformational_sampling.utils import stk_mol_to_ase_atoms


@dataclass
class ASE(stko.optimizers.Optimizer):
    calculator: Calculator

    def optimize(self, stk_mol: stk.Molecule):
        ase_mol = stk_mol_to_ase_atoms(stk_mol)
        ase_mol.calc = self.calculator
        opt = BFGS(ase_mol)
        try:
            opt.run(fmax=0.1)
        except CalculationFailed:
            return None
        return stk_mol.with_position_matrix(opt.atoms.get_positions())