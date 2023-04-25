import ase
import stko
from ase.calculators.calculator import Calculator
from ase.optimize import BFGS


from dataclasses import dataclass


@dataclass
class ASE(stko.optimizers.Optimizer):
    calculator: Calculator

    def optimize(self, stk_mol):
        ase_mol = ase.Atoms(
            positions=list(stk_mol.get_atomic_positions()),
            numbers=[atom.get_atomic_number() for atom in stk_mol.get_atoms()]
        )
        ase_mol.calc = self.calculator
        opt = BFGS(ase_mol)
        opt.run(fmax=0.2)
        stk_mol.with_position_matrix(opt.atoms.get_positions())
        return stk_mol