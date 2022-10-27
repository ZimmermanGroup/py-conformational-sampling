from ase.build import molecule
from ase.calculators.qchem import QChem
from ase.optimize import LBFGS
from ase.io import read, write
import ase
import os
from pathlib import Path
from rdkit.Chem.rdmolfiles import MolToXYZBlock

import stk

os.environ['QCSCRATCH'] = 'qc_scratch'
stk_mol = stk.BuildingBlock(smiles='[Ne][P](F)(F)(F)')
ase_mol = ase.Atoms(
    positions=list(stk_mol.get_atomic_positions()),
    numbers=[atom.get_atomic_number() for atom in stk_mol.get_atoms()]
)
write('test.xyz', ase_mol)
# rdkit_mol = stk_mol.to_rdkit_mol()
# MolToXYZBlock(rdkit_mol)
# mol = read('NePF3.xyz')
# mol = molecule('C2H6')
calc = QChem(
    label='NePF3',
    method='B3LYP',
    basis='6-31+G*',
    nt=4
)
ase_mol.calc = calc
opt = LBFGS(ase_mol)
opt.run()
print(opt.s) # list 
