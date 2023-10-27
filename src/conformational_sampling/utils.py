import os
import numpy as np
from scipy.constants import gas_constant
import ase
from openbabel import pybel as pb
from rdkit.Chem.rdmolfiles import MolToXYZBlock, MolFromMolBlock, MolToMolBlock

import stk
from rdkit import Chem
from pygsm.utilities.units import KJ_PER_KCAL

def num_cpus():
    try:
        return int(os.environ["SLURM_NTASKS_PER_NODE"])
    except KeyError:
        return 2


def pybel_mol_to_rdkit_mol(pybel_mol):
    rdkit_mol = MolFromMolBlock(pybel_mol.write('mol'), removeHs=False)
    try:
        Chem.rdmolops.Kekulize(rdkit_mol)
        return rdkit_mol
    except:
        return rdkit_mol


def rdkit_mol_to_stk_mol(rdkit_mol):
    return stk.BuildingBlock.init_from_rdkit_mol(rdkit_mol)


def pybel_mol_to_stk_mol(pybel_mol):
    return rdkit_mol_to_stk_mol(pybel_mol_to_rdkit_mol(pybel_mol))


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
    

def stk_metal(metal: str) -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles=f'[{metal}]',
        functional_groups=(
            stk.SingleAtom(stk.Pd(0))
            for i in range(6)
        ),
        position_matrix=[[0, 0, 0]],
    )


def partition_function(energies, temperature) -> float:
    'energies are assumed to be in kcal/mol; temperature is assumed to be in Kelvin'
    energies = np.array(energies)
    return sum(np.exp(-1 * energies * KJ_PER_KCAL * 1000 / (gas_constant  * temperature)))
    

def free_energy_diff(energies_1, energies_2, temperature: float) -> float:
    'free energy returned in kcal/mol; free energy 2 - free energy 1'
    partition_function_1 = partition_function(energies_1, temperature)
    partition_function_2 = partition_function(energies_2, temperature)
    return (gas_constant / (KJ_PER_KCAL * 1000) * temperature
            * (np.log(partition_function_1) - np.log(partition_function_2)))
    
