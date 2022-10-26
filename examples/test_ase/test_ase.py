from ase.build import molecule
from ase.calculators.qchem import QChem
from ase.optimize import LBFGS
from ase.io import read, write
import os
from pathlib import Path

os.environ['QCSCRATCH'] = 'qc_scratch'
mol = read('NePF3.xyz')
# mol = molecule('C2H6')
calc = QChem(label='ethane',
             method='B3LYP',
             basis='6-31+G*')
mol.calc = calc
opt = LBFGS(mol)
opt.run()
