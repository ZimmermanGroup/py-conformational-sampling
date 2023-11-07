import numpy as np
import stk
from xtb.ase import calculator
from xtb.interface import Calculator, Param, Environment
from xtb.libxtb import VERBOSITY_FULL

from conformational_sampling.ase_stko_optimizer import ASE
from conformational_sampling.utils import stk_mol_to_ase_atoms


def test_xtb_solvent():
    stk_mol = stk.BuildingBlock('CCCC')
    # ase_mol.calc = calculator.XTB(solvent='gbsa')
    ase_mol = stk_mol_to_ase_atoms(stk_mol)
    env = Environment()
    env.set_output('error.log')
    numbers = np.array([8, 1, 1])
    positions = np.array([
    [ 0.00000000000000, 0.00000000000000,-0.73578586109551],
    [ 1.44183152868459, 0.00000000000000, 0.36789293054775],
    [-1.44183152868459, 0.00000000000000, 0.36789293054775]
    ])
    calc = Calculator(Param.GFN2xTB, numbers=numbers, positions=positions)
    env.set_verbosity(VERBOSITY_FULL)
    calc.set_verbosity(VERBOSITY_FULL)
    calc.singlepoint()
    env.release_output()
    
    # ase_mol.calc = calculator.XTB(solvent='acetonitrile', write='test.out')
    # energy = ase_mol.get_potential_energy()
    # print(energy)
    pass

# -371.55510438052335
test_xtb_solvent()