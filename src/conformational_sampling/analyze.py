#!/usr/bin/python

from itertools import combinations
import os
from pathlib import Path
import re

from openbabel import pybel as pb
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from conformational_sampling.utils import pybel_mol_to_rdkit_mol


# Function to determine the existence and position of a unique TS node along the string
 
def ts_node(en_list):
    
    res_list = list(combinations(en_list, 2))

    max_diff = 0.0
    ts_node_energy = None
    ts_barrier = 0.0
    for item in res_list:
        tmp = item[1] - item[0]
        if tmp > max_diff:
            max_diff = tmp
            ts_node_energy = item[1]
            ts_barrier = max_diff
            
    return (max_diff, ts_node_energy, ts_barrier)


def analyze():
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

    # print(sorted(file_list))

    infofile = open('TS_and_stereo.txt', 'w')
    infofile.write("   Conformers      TS_Energy_Barrier     Stereochemistry   Torsion angle (deg)\n")

    ensemble_TS = open('ensemble_TS.xyz', 'w')

    # Conformer Stereochemistry: Dihedral C(OMe)-C(Ph)-C(Ph)-C(Me)
    # -180 to 0 --> R and 0 to 180 --> S




    for file in file_list:
        with open(file, 'r') as f:
            flines = f.read().splitlines()

            try:
                natoms = int(flines[0])
                nlines = natoms + 2

                line = 1
                en_list = []
                while line < len(flines):
                    en_list.append(float(flines[line]))
                    line += nlines
            
            except:
                natoms = int(flines[2])
                nlines = natoms + 2

                start_index = flines.index('energy')
                end_index = flines.index('max-force')

                en_list = []

                for i in range(start_index+1, end_index):
                    en_list.append(float(flines[i]))

                with open('file_xyz.xyz', 'w') as f_xyz:
                    f_xyz.write('\n'.join(flines[2:start_index-1]))    
            
            if ts_node(en_list)[0] != 0.0:

                try:
                    pybel_product = list(pb.readfile(format='xyz', filename='file_xyz.xyz'))[-1]
                    TS_index = en_list.index(ts_node(en_list)[1])*nlines + 2
                except:
                    pybel_product = list(pb.readfile(format='xyz', filename=file))[-1]
                    TS_index = en_list.index(ts_node(en_list)[1])*nlines

                rdkit_product = pybel_mol_to_rdkit_mol(pybel_product)
                torsion_deg = rdMolTransforms.GetDihedralDeg(rdkit_product.GetConformer(), 56, 55, 79, 78)
                if torsion_deg >= 0:
                    stereochem = 'S'
                else:
                    stereochem = 'R'
                # 
                idx = re.search(r"pystring_(\d+)", file).groups()[0]
                
                infofile.write(f"  Conformer_{idx}     {ts_node(en_list)[2]:.6f}      {stereochem}           {torsion_deg:.3f}\n")

                ensemble_TS.write('\n'.join(flines[TS_index:TS_index+nlines]))
                ensemble_TS.write('\n')

    ensemble_TS.close()
    try:
        os.remove('file_xyz.xyz')
    except:
        pass
        
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

if __name__ == '__main__':
    analyze()

