#!/usr/bin/python

import os
from pathlib import Path
import re

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
        if 'grown_string1_000.xyz' in os.listdir(in_path):
            file_list.append(str(in_path)+"/grown_string1_000.xyz")
        # new_path = Path(f'scratch/{names}/scratch')
        # for name in os.listdir(new_path):
        #     if name.startswith("pystring"):
        #         in_in_path = Path(f'scratch/{names}/scratch/{name}')
        #         if 'opt_converged_000.xyz' in os.listdir(in_in_path):
        #             file_list.append(str(in_in_path)+"/opt_converged_000.xyz")

# print(sorted(file_list))

infofile = open('TS_and_stereo.txt', 'w')
infofile.write("   Conformers      TS_Energy     Stereochemistry   Torsion angle (deg)\n")

ensemble_TS = open('ensemble_TS.xyz', 'w')

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
        
        # todo: improve mechanism for determining which node is the TS
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

            TS_index = en_list.index(max(en_list))*nlines

            ensemble_TS.write('\n'.join(flines[TS_index:TS_index+nlines]))
            ensemble_TS.write('\n')

ensemble_TS.close()
    
## Script for submitting the TS jobs
## subprocess.run(['sbatch', f'--array=0-{len(file_list)-1}', './tests/ts_job_array.py'])

## Script for making the TS job input files

new_path = Path('OptTS')

with open('ensemble_TS.xyz', 'r') as f:
    flines = f.read().splitlines()

    natoms = int(flines[0])
    nlines = natoms + 2

    os.chdir(new_path)

    for i, n in enumerate(range(0, len(flines), nlines)):
        with open('qstart2.inp', 'r') as f_1, open(f'q{i:03d}.TS.inp', 'w') as f_2:
            for lines in f_1:
                f_2.write(lines)
        with open(f'q{i:03d}.TS.inp', 'a+') as f_2:
            f_2.write('\n'.join(flines[n+2:n+nlines]))
            f_2.write('\n')
        with open(f'q{i:03d}.TS.inp', 'a+') as f_2, open('qend2.inp', 'r') as f_3:
            for lines in f_3:
                f_2.write(lines)



        


