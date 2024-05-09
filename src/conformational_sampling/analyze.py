#!/usr/bin/python

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import List, Optional

import numpy as np
import stk
from openbabel import pybel as pb
from pyGSM.utilities.units import EV_TO_AU, KCAL_MOL_PER_AU
from rdkit import Chem
from rdkit.Chem import rdmolops, rdMolTransforms
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock

from conformational_sampling.utils import rdkit_mol_to_stk_mol


@dataclass
class System:
    reductive_elim_torsion: tuple
    pro_dis_torsion: tuple
    mol_path: Path
    atrop_torsion: Optional[tuple] = None

systems = {
    # 'ligand_l1': System(
    #     reductive_elim_torsion=(56, 55, 79, 78),
    #     pro_dis_torsion=(21, 11, 55, 65),
    #     mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l1_xtb')
    # ),
    'ligand_achiral': System(
        reductive_elim_torsion=(36, 35, 59, 58),
        pro_dis_torsion=(21, 11, 35, 45),
        mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l1_symm_xtb_crest')
    ),
    'ligand_l3': System(
        reductive_elim_torsion=(78, 77, 101, 100),
        pro_dis_torsion=(17, 38, 77, 87),
        atrop_torsion=(3, 9, 17, 19),
        mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l3_xtb')
    ),
    'ligand_l4': System(
        reductive_elim_torsion=(74, 73, 97, 96),
        pro_dis_torsion=(11, 34, 73, 83),
        mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l4_xtb_new')
    ),
    'ligand_l6': System(
        reductive_elim_torsion=(82, 81, 105, 104),
        pro_dis_torsion=(11, 42, 81, 91),
        mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l6_xtb_new')
    ),
    'ligand_l8': System(
        reductive_elim_torsion=(74, 73, 97, 96),
        pro_dis_torsion=(47, 9, 73, 83),
        mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l8_xtb_new')
    ),
    # 'ligand_l8_dft': System(
    #     reductive_elim_torsion=(74, 73, 97, 96),
    #     pro_dis_torsion=(47, 9, 73, 83),
    #     mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l8_xtb')
    # ),
}

exclude_confs = {
    'ligand_l4': [315]
}


# Function to determine the existence and position of a unique TS node along the string
 
def ts_node(en_list):
    
    res_list = list(combinations(en_list, 2))

    max_diff = 0.0
    ts_node_energy = en_list[0]
    ts_barrier = 0.0
    min_reac_energy = en_list[0]
    for item in res_list:
        tmp = item[1] - item[0]
        if tmp > max_diff:
            max_diff = tmp
            ts_node_energy = item[1]
            min_reac_energy = item[0]
            ts_barrier = max_diff
            
    return (max_diff, ts_node_energy, ts_barrier, min_reac_energy)


def truncate_string_at_bond_formation(
    string_nodes: List[Chem.rdchem.Mol], atom_1: int, atom_2: int, thresh=1.7
):
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
    singlepoints: bool = False

    def __post_init__(self):
        ob_string_nodes = list(pb.readfile("xyz", str(self.string_path)))
        if self.singlepoints:
            self.string_nodes = []
            self.string_energies = []
            for node in ob_string_nodes:
                mol_block = node.write('mol')
                self.string_nodes.append(MolFromMolBlock(mol_block, removeHs=False))
                self.string_energies.append(float(mol_block.split()[0]) * EV_TO_AU * KCAL_MOL_PER_AU)
        else:
            self.string_nodes = [
                MolFromMolBlock(node.write("mol"), removeHs=False)
                for node in ob_string_nodes
            ]  # convert to rdkit

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
        max_diff, self.ts_energy, self.activation_energy, self.min_reac_energy = ts_node(self.string_energies[:len(self.truncated_string)])
        self.ts_node_num = self.string_energies.index(self.ts_energy)
        self.reac_node_num = self.string_energies.index(self.min_reac_energy)
        self.reac_rdkit_mol = self.string_nodes[self.reac_node_num]
        self.ts_rdkit_mol = self.string_nodes[self.ts_node_num]
        self.pdt_rdkit_mol = self.string_nodes[-1]

        # compute properties of the reactant
        phos = [atom for atom in self.reac_rdkit_mol.GetAtoms() if atom.GetSymbol() == 'P'][0]
        self.improper_torsion = rdMolTransforms.GetDihedralDeg(
            self.reac_rdkit_mol.GetConformer(),
            self.system.reductive_elim_torsion[1],
            0,
            phos.GetIdx(),
            self.system.reductive_elim_torsion[2],
        )
        self.ligands_angle = rdMolTransforms.GetAngleDeg(
            self.reac_rdkit_mol.GetConformer(),
            self.system.reductive_elim_torsion[1],
            0,
            self.system.reductive_elim_torsion[2],
        )
        
        self.out_of_plane_angle = out_of_plane_angle(
            self.reac_rdkit_mol,
            0,
            *self.system.reductive_elim_torsion[1:3],
            phos.GetIdx(),
        )
        
        self.tau_4_prime = tau_4_prime(self.reac_rdkit_mol, 0)
        
        # compute properties of the transition state
        self.forming_bond_torsion = rdMolTransforms.GetDihedralDeg(
            self.ts_rdkit_mol.GetConformer(),
            *self.system.reductive_elim_torsion
        )
        self.pro_dis_torsion = rdMolTransforms.GetDihedralDeg(
            self.ts_rdkit_mol.GetConformer(),
            *self.system.pro_dis_torsion
        )
        self.improper_torsion_ts = rdMolTransforms.GetDihedralDeg(
            self.ts_rdkit_mol.GetConformer(),
            self.system.reductive_elim_torsion[1],
            0,
            phos.GetIdx(),
            self.system.reductive_elim_torsion[2],
        )
        self.tau_4_prime_ts = tau_4_prime(self.ts_rdkit_mol, 0)
        
        #compute properties of the product 
        self.formed_bond_torsion = rdMolTransforms.GetDihedralDeg(
            self.pdt_rdkit_mol.GetConformer(),
            *self.system.reductive_elim_torsion
        )

        self.pro_dis = 'proximal' if -90 <= self.pro_dis_torsion <= 90 else 'distal'
        # ts is exo if the torsion of the bond being formed is positive and the ts is proximal
        # if distal, the relationship is reversed
        self.endo_exo = 'exo' if (self.forming_bond_torsion >= 0) ^ (self.pro_dis == 'distal') else 'endo'
        self.syn_anti = 'syn' if -90 <= self.forming_bond_torsion <= 90 else 'anti'
        self.pdt_stereo = 'R' if self.formed_bond_torsion <= 0 else 'S'
        

def tau_4_prime(rdkit_mol, atom_id: int = 0) -> float:
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


def out_of_plane_angle(
    rdkit_mol,
    central_atom: int,
    plane_atom_1: int,
    plane_atom_2: int,
    out_of_plane_atom: int,
):
    stk_mol = rdkit_mol_to_stk_mol(rdkit_mol)
    plane_normal = stk_mol.get_plane_normal((central_atom, plane_atom_1, plane_atom_2))
    pos_matrix = stk_mol.get_position_matrix()
    out_of_plane_vector = pos_matrix[out_of_plane_atom] - pos_matrix[central_atom]
    plane_normal = stk.get_acute_vector(out_of_plane_vector, plane_normal)
    return 90 - np.rad2deg(stk.vector_angle(plane_normal, out_of_plane_vector))