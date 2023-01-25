from ase.build import molecule
from ase.calculators.qchem import QChem
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.io import read, write
import ase
from ase.io.trajectory import Trajectory
import os
from pathlib import Path
from rdkit.Chem.rdmolfiles import MolToXYZBlock
from xtb.ase.calculator import XTB

import stk

from conformational_sampling.main import stk_list_to_xyz_file

os.environ['QCSCRATCH'] = 'qc_scratch'

stk_mol = stk.BuildingBlock(smiles='[Ne][P](F)(F)(F)')
ase_mol = ase.Atoms(
    positions=list(stk_mol.get_atomic_positions()),
    numbers=[atom.get_atomic_number() for atom in stk_mol.get_atoms()]
)
# calc = QChem(
#     label='ethane',
#     method='PBE',
#     basis='6-31G',
#     nt=4
# )

calc = XTB()

ase_mol.calc = calc
opt = BFGS(ase_mol, trajectory='test.traj')
opt.run()
trajectory = Trajectory('test.traj')
stk_trajectory = [stk_mol.with_position_matrix(atoms.get_positions()) for atoms in trajectory]
stk_list_to_xyz_file(stk_trajectory, 'test_trajectory.xyz')