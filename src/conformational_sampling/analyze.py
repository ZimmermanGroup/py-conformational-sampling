#!/usr/bin/python

from itertools import combinations
import os
from pathlib import Path
import re

from openbabel import pybel as pb
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdmolops

from conformational_sampling.utils import pybel_mol_to_rdkit_mol


# Function to determine the existence and position of a unique TS node along the string
 
def ts_node(en_list):
    
    res_list = list(combinations(en_list, 2))

    max_diff = 0.0
    ts_node_energy = en_list[0]
    ts_barrier = 0.0
    for item in res_list:
        tmp = item[1] - item[0]
        if tmp > max_diff:
            max_diff = tmp
            ts_node_energy = item[1]
            ts_barrier = max_diff
            
    return (max_diff, ts_node_energy, ts_barrier)


def truncate_string_at_bond_formation(string_nodes: list[Chem.rdchem.Mol], atom_1: int, atom_2: int, thresh=1.8):
    for idx, node in enumerate(string_nodes):
        dist_matrix = rdmolops.Get3DDistanceMatrix(node)
        if dist_matrix[atom_1][atom_2] < thresh:
            return string_nodes[:idx]
    return string_nodes
