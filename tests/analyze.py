#!/usr/bin/python

import os
from pathlib import Path
import numpy as np

def dihedral(p):
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))  

path = Path('scratch')
file_list = []

for names in os.listdir(path):
    if names.startswith("pystring"):
        #dir_list.append(names)
        in_path = Path(f'scratch/{names}')
        if 'opt_converged_000.xyz' in os.listdir(in_path):
            file_list.append(str(in_path)+"/opt_converged_000.xyz")
        new_path = Path(f'scratch/{names}/scratch')
        # for name in os.listdir(new_path):
        #     if name.startswith("pystring"):
        #         in_in_path = Path(f'scratch/{names}/scratch/{name}')
        #         if 'opt_converged_000.xyz' in os.listdir(in_in_path):
        #             file_list.append(str(in_in_path)+"/opt_converged_000.xyz")

#print(sorted(file_list))

infofile = open('TS_and_stereo.txt', 'a+')
infofile.write("   Conformers      TS_Energy     Stereochemistry\n")

# Conformer Stereochemistry: Dihedral C(OMe)-C(Ph)-C(Ph)-C(Me)
# -180 to 0 --> R and 0 to 180 --> S


idx = 1
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

            infofile.write(f"  Conformer_{idx:02d}     {max(en_list)}           \n")
        
    idx += 1


## Script for submitting the TS jobs