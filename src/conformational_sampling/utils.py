import os
import ase
from openbabel import pybel as pb
from rdkit.Chem.rdmolfiles import MolToXYZBlock, MolFromMolBlock

import stk
from rdkit import Chem

def num_cpus():
    try:
        return int(os.environ["SLURM_CPUS_PER_TASK"])
    except KeyError:
        return 2


def pybel_mol_to_stk_mol(pybel_mol):
    rdkit_mol = MolFromMolBlock(pybel_mol.write('mol'), removeHs=False)
    Chem.rdmolops.Kekulize(rdkit_mol)
    stk_mol = stk.BuildingBlock.init_from_rdkit_mol(rdkit_mol)
    return stk_mol


def stk_mol_to_pybel_mol(stk_mol, reperceive_bonds=False):
    if reperceive_bonds:
        return pb.readstring('xyz', MolToXYZBlock(stk_mol.to_rdkit_mol()))
    else:
        return pb.readstring('mol', MolToMolBlock(stk_mol.to_rdkit_mol()))


def stk_mol_to_ase_atoms(stk_mol: stk.Molecule) -> ase.Atoms:
    return ase.Atoms(
        positions=list(stk_mol.get_atomic_positions()),
        numbers=[atom.get_atomic_number() for atom in stk_mol.get_atoms()]
    )
    