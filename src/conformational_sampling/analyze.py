#!/usr/bin/python

import os
from pathlib import Path
import re
import numpy as np

from openbabel import pybel as pb
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from conformational_sampling.utils import pybel_mol_to_rdkit_mol


path = Path('scratch')
# path = Path('/export/zimmerman/soumikd/py-conformational-sampling/scratch')
file_list = []

for name in os.listdir(path):
    if name.startswith("pystring"):
        #dir_list.append(names)
        in_path = path / f'{name}'
        if 'opt_converged_000.xyz' in os.listdir(in_path):
            file_list.append(str(in_path)+"/opt_converged_000.xyz")
        # new_path = Path(f'scratch/{names}/scratch')
        # for name in os.listdir(new_path):
        #     if name.startswith("pystring"):
        #         in_in_path = Path(f'scratch/{names}/scratch/{name}')
        #         if 'opt_converged_000.xyz' in os.listdir(in_in_path):
        #             file_list.append(str(in_in_path)+"/opt_converged_000.xyz")

#print(sorted(file_list))

infofile = open('TS_and_stereo.txt', 'w')
infofile.write("   Conformers      TS_Energy     Stereochemistry   Torsion angle (deg)\n")

# Conformer Stereochemistry: Dihedral C(OMe)-C(Ph)-C(Ph)-C(Me)
# -180 to 0 --> R and 0 to 180 --> S


for file in file_list:
    with open(file, 'r') as f:
        flines = f.read().splitlines()

        natoms = int(flines[0])
        nlines = natoms + 2

        line = 1
        en_list = []
        while line < len(flines):
            en_list.append(float(flines[line]))
            line += nlines
        
        if max(en_list) != 0.0:

            pybel_product = list(pb.readfile(format='xyz', filename=file))[-1]
            rdkit_product = pybel_mol_to_rdkit_mol(pybel_product)
            torsion_deg = rdMolTransforms.GetDihedralDeg(rdkit_product.GetConformer(), 56, 55, 79, 78)
            if torsion_deg >= 0:
                stereochem = 'S'
            else:
                stereochem = 'R'
            # 
            idx = re.search(r"pystring_(\d+)", file).groups()[0]
            
            infofile.write(f"  Conformer_{idx}     {max(en_list)}      {stereochem}           {torsion_deg}\n")


## Script for submitting the TS jobs