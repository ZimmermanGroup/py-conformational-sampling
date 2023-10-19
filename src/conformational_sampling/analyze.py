#!/usr/bin/python

from dataclasses import dataclass
from itertools import combinations
import os
from pathlib import Path
import re

from openbabel import pybel as pb
from pygsm.utilities.units import KCAL_MOL_PER_AU
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdmolops
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock

from conformational_sampling.utils import pybel_mol_to_rdkit_mol

@dataclass
class System:
    reductive_elim_torsion: tuple
    pro_dis_torsion: tuple
    mol_path: Path


systems = {
    'ligand_l1': System(
        reductive_elim_torsion=(56, 55, 79, 78),
        pro_dis_torsion=(21, 11, 55, 65),
        mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l1_xtb')
    ),
    'ligand_achiral': System(
        reductive_elim_torsion=(36, 35, 59, 58),
        pro_dis_torsion=(21, 11, 35, 45),
        mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l1_symm_xtb_crest')
    ),
    'ligand_l3': System(
        reductive_elim_torsion=(78, 77, 101, 100),
        pro_dis_torsion=(17, 38, 77, 87),
        mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l3_xtb')
    ),
    'ligand_l4': System(
        reductive_elim_torsion=(74, 73, 97, 96),
        pro_dis_torsion=(11, 34, 73, 83),
        mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l4_xtb')
    ),
    'ligand_l6': System(
        reductive_elim_torsion=(82, 81, 105, 104),
        pro_dis_torsion=(11, 42, 81, 91),
        mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l6_xtb')
    ),
    'ligand_l8': System(
        reductive_elim_torsion=(74, 73, 97, 96),
        pro_dis_torsion=(47, 9, 73, 83),
        mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l8_xtb_crest')
    ),
    # 'ligand_l8_dft': System(
    #     reductive_elim_torsion=(74, 73, 97, 96),
    #     pro_dis_torsion=(47, 9, 73, 83),
    #     mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l8_xtb')
    # ),
}


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


def truncate_string_at_bond_formation(string_nodes: list[Chem.rdchem.Mol], atom_1: int, atom_2: int, thresh=1.7):
    for idx, node in enumerate(string_nodes):
        dist_matrix = rdmolops.Get3DDistanceMatrix(node)
        if dist_matrix[atom_1][atom_2] < thresh:
            return string_nodes[:idx]
    return string_nodes


# setup_mol()
@dataclass
class Conformer:
    system: System
    string_path: Path

    def __post_init__(self):
        string_nodes = list(pb.readfile('xyz', str(self.string_path)))
        self.string_nodes = [MolFromMolBlock(node.write('mol'), removeHs=False)
                             for node in string_nodes] # convert to rdkit

        # previous way of getting reactant energy
        # raw_dft_energy_path = self.string_path.parent / 'scratch' / '000' / 'E_0.txt'
        # raw_dft_energy_au = float(raw_dft_energy_path.read_text().split()[2])
        # self.string_energies = [(raw_dft_energy_au + float(MolToMolBlock(node).split()[0])) * KCAL_MOL_PER_AU
        #                         for node in self.string_nodes]

        # getting reactant energy
        raw_dft_energy_path = self.string_path.parent / 'de_output.txt'
        with open(raw_dft_energy_path) as file:
            while line := file.readline():
                if line.startswith(' Energy of the end points are'):
                    # extracting energy from output file formatting
                    raw_dft_energy_kcal_mol = float(line.split()[-2][:-1])
        self.string_energies = [raw_dft_energy_kcal_mol + float(MolToMolBlock(node).split()[0]) * KCAL_MOL_PER_AU
                                for node in self.string_nodes]

        self.truncated_string = truncate_string_at_bond_formation(
            self.string_nodes,
            *self.system.reductive_elim_torsion[1:3] # middle atoms of torsion
        )
        if not self.truncated_string:
            return
        max_diff, self.ts_energy, self.activation_energy = ts_node(self.string_energies[:len(self.truncated_string)])
        self.ts_node_num = self.string_energies.index(self.ts_energy)
        self.ts_rdkit_mol = self.string_nodes[self.ts_node_num]
        self.pdt_rdkit_mol = self.string_nodes[-1]

        # compute properties of the transition state
        self.forming_bond_torsion = rdMolTransforms.GetDihedralDeg(
            self.ts_rdkit_mol.GetConformer(),
            *self.system.reductive_elim_torsion
        )
        self.pro_dis_torsion = rdMolTransforms.GetDihedralDeg(
            self.ts_rdkit_mol.GetConformer(),
            *self.system.pro_dis_torsion
        )
        #compute properties of the product 
        self.formed_bond_torsion = rdMolTransforms.GetDihedralDeg(
            self.pdt_rdkit_mol.GetConformer(),
            *self.system.reductive_elim_torsion
        )

        self.tau_4_prime = tau_4_prime(self.ts_rdkit_mol, 0)

        self.pro_dis = 'proximal' if -90 <= self.pro_dis_torsion <= 90 else 'distal'
        # ts is exo if the torsion of the bond being formed is positive and the ts is proximal
        # if distal, the relationship is reversed
        self.endo_exo = 'exo' if (self.forming_bond_torsion >= 0) ^ (self.pro_dis == 'distal') else 'endo'
        self.syn_anti = 'syn' if -90 <= self.forming_bond_torsion <= 90 else 'anti'
        self.pdt_stereo = 'R' if self.formed_bond_torsion <= 0 else 'S'
        

def tau_4_prime(rdkit_mol, atom_id: int) -> float:
    'formula from https://en.wikipedia.org/wiki/Geometry_index'
    
    rdkit_atom = rdkit_mol.GetAtomWithIdx(atom_id)
    neighbors = rdkit_atom.GetNeighbors()
    if len(neighbors) != 4:
        return -1.0
    neighbors = [neighbor.GetIdx() for neighbor in neighbors]

    # list of tuples containing atom ids of each angle through the atom
    angles_atoms = [(comb[0], atom_id, comb[1]) for comb in combinations(neighbors, 2)]
    angles = [rdMolTransforms.GetAngleDeg(rdkit_mol.GetConformer(), *angle_atoms)
              for angle_atoms in angles_atoms]
    alpha, beta = sorted(angles, reverse=True)[:2]
    return -0.00399 * alpha - 0.01019 * beta + 2.55