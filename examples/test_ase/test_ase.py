from ase.build import molecule
from ase.calculators.qchem import QChem
from ase.calculators.emt import EMT
from ase.optimize import LBFGS
from ase.io import read, write
import ase
import os
from pathlib import Path
from rdkit.Chem.rdmolfiles import MolToXYZBlock

import stk

from conformational_sampling.main import stk_list_to_xyz_file

os.environ['QCSCRATCH'] = 'qc_scratch'

stk_mol = stk.BuildingBlock(smiles='[Ne][P](F)(F)(F)')
ase_mol = ase.Atoms(
    positions=list(stk_mol.get_atomic_positions()),
    numbers=[atom.get_atomic_number() for atom in stk_mol.get_atoms()]
)
# write('test.xyz', ase_mol)
# calc = QChem(
#     label='NePF3',
#     method='PBE',
#     basis='6-31G',
#     nt=4
# )

ase_mol = molecule('H2O')
calc = EMT()

ase_mol.calc = calc
opt = LBFGS(ase_mol)
opt.run()
print(opt.s) # list 
stk_trajectory = [stk_mol.with_position_matrix(pos_matrix) for pos_matrix in opt.s] # shapes of ndarrays don't match
stk_list_to_xyz_file(stk_trajectory, 'test_trajectory.xyz')