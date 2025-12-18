import os
# MUST be set before importing xtb to prevent OpenMP thread oversubscription
# when parallelizing over conformers. See README for details.
os.environ["OMP_NUM_THREADS"] = "1"

import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

FORMAT = (
    '[%(asctime)s %(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s'
)
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

from rdkit import Chem
from openbabel import pybel as pb
from xtb.ase import calculator

from conformational_sampling.main import (
    load_stk_mol, 
    gen_confs_openbabel, 
    stk_list_to_xyz_file,
    xtb_energy,
    ConformerOptimizationSequence,
    UNOPTIMIZED,
    METAL_OPTIMIZER,
    XTB,
)
from conformational_sampling.config import Config
from conformational_sampling.ase_stko_optimizer import ASE
from conformational_sampling.utils import stk_mol_to_pybel_mol, pybel_mol_to_stk_mol


def replace_boron_with_hydrogen(stk_mol):
    """Replace boron atom with hydrogen in an stk molecule.
    
    Changes the atomic number from B (5) to H (1) and removes the extra bond
    to nitrogen (keeping only O-H). Does NOT reperceive bonds to avoid 
    connectivity changes from atoms being too close together.
    """
    # Convert to pybel to modify
    pybel_mol = stk_mol_to_pybel_mol(stk_mol)
    obmol = pybel_mol.OBMol
    
    # Find boron atom
    boron_atom = None
    boron_idx = None
    for atom in pb.ob.OBMolAtomIter(obmol):
        if atom.GetAtomicNum() == 5:  # Boron
            boron_atom = atom
            boron_idx = atom.GetIdx()
            break
    
    if boron_atom is None:
        logging.warning('No boron atom found in molecule')
        return stk_mol
    
    # Find bonds to boron and identify O and N neighbors
    bonds_to_remove = []
    for bond in pb.ob.OBMolBondIter(obmol):
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if boron_idx in (begin_idx, end_idx):
            other_idx = end_idx if begin_idx == boron_idx else begin_idx
            other_atom = obmol.GetAtom(other_idx)
            # Keep bond to O (atomic num 8), remove bond to N (atomic num 7)
            if other_atom.GetAtomicNum() == 7:  # Nitrogen
                bonds_to_remove.append(bond)
                logging.debug(f'Will remove B-N bond (B idx {boron_idx}, N idx {other_idx})')
    
    # Remove the B-N bond(s)
    for bond in bonds_to_remove:
        obmol.DeleteBond(bond)
    
    # Change B to H
    boron_atom.SetAtomicNum(1)
    logging.debug(f'Changed B at index {boron_idx} to H')
    
    # Convert back to stk molecule
    return pybel_mol_to_stk_mol(pybel_mol)


def compute_rmsd(mol1, mol2):
    """Compute RMSD between two stk molecules with identical atom ordering.
    
    Uses RDKit's AlignMol with an explicit atom map (identity mapping)
    since conformers have identical atom indices - no substructure matching needed.
    """
    rdkit_mol1 = mol1.to_rdkit_mol()
    rdkit_mol2 = mol2.to_rdkit_mol()
    
    num_atoms = rdkit_mol1.GetNumAtoms()
    if num_atoms != rdkit_mol2.GetNumAtoms():
        raise ValueError(f"Molecules have different number of atoms: {num_atoms} vs {rdkit_mol2.GetNumAtoms()}")
    
    # Create identity atom map: [(0,0), (1,1), (2,2), ...]
    atom_map = [(i, i) for i in range(num_atoms)]
    
    # AlignMol returns RMSD after optimal alignment using the provided atom map
    rmsd = Chem.AllChem.AlignMol(rdkit_mol1, rdkit_mol2, atomMap=atom_map)
    return rmsd


def get_unique_conformer_ids(conformers, stage, rms_threshold=0.175):
    """Get indices of unique conformers based on RMSD comparison.
    
    Uses direct coordinate-based RMSD since all conformers have identical
    atom indices (no substructure matching needed).
    
    Args:
        conformers: List of ConformerOptimizationSequence objects
        stage: Stage name to use for comparison
        rms_threshold: RMSD threshold in Ångströms (default 0.175 Å)
                      - Smaller (0.1-0.2 Å): keeps MORE conformers
                      - Larger (0.5-2.0 Å): keeps FEWER conformers
    """
    stk_mols = {i: conformer.stages[stage]
                for i, conformer in enumerate(conformers)
                if stage in conformer.stages}
    
    unique_indices = []
    for i, stk_mol in stk_mols.items():
        is_unique = True
        for j in unique_indices:
            rmsd = compute_rmsd(stk_mol, stk_mols[j])
            if rmsd < rms_threshold:
                is_unique = False
                break
        if is_unique:
            unique_indices.append(i)
    return unique_indices


# Load the full molecular system with boron bridge
full_system_path = Path('full_system.xyz')
stk_mol = load_stk_mol(full_system_path)

# Configuration
config = Config(
    # larger value samples more thoroughly, set to roughly 20-100 for thorough
    # sampling depending on ligand flexibility and resource availability
    initial_conformers=100,
    num_cpus=64,  # set to number of cpus available for use by pycosa
    # RMSD threshold for removing duplicate conformers before xTB optimization
    # 0.175 Å removes only very similar structures (keeps MORE conformers)
    # Increase to 0.5-2.0 Å for more aggressive filtering (keeps FEWER conformers)
    pre_xtb_rms_threshold=0.175,
)

# Step 1: Generate conformers using OpenBabel with boron bridge
logging.debug('Step 1: Generating conformers with boron bridge')
stk_conformers_with_boron = gen_confs_openbabel(stk_mol, config)
logging.debug(f'Generated {len(stk_conformers_with_boron)} conformers')

# Save intermediate result
stk_list_to_xyz_file(stk_conformers_with_boron, 'conformers_0_with_boron.xyz')

# Step 2: Replace boron with hydrogen in all conformers
logging.debug('Step 2: Replacing boron with hydrogen in conformers')
stk_conformers = [replace_boron_with_hydrogen(conf) for conf in stk_conformers_with_boron]
logging.debug(f'Replaced boron with hydrogen in {len(stk_conformers)} conformers')

# Create conformer sequences
conformers = [ConformerOptimizationSequence(conformer) for conformer in stk_conformers]

# Save unoptimized with hydrogen
stk_list_to_xyz_file(stk_conformers, 'conformers_0_unoptimized.xyz')
logging.debug('Saved unoptimized conformers (with hydrogen)')

# Step 3: Remove duplicates before xTB (skip force field optimization)
logging.debug('Step 3: Removing duplicate conformers before xTB')
with ProcessPoolExecutor(max_workers=config.num_cpus) as executor:
    # Store unoptimized as the pre-xTB stage for duplicate detection
    for conformer in conformers:
        conformer.stages[METAL_OPTIMIZER] = conformer.stages[UNOPTIMIZED]
    
    unique_ids = get_unique_conformer_ids(conformers, METAL_OPTIMIZER, config.pre_xtb_rms_threshold)
    unique_conformers = [conformers[i] for i in unique_ids]
    logging.debug(f'Kept {len(unique_conformers)} unique conformers (removed {len(conformers) - len(unique_conformers)} duplicates)')
    
    unoptimized_unique = [conformer.stages[UNOPTIMIZED] for conformer in unique_conformers]
    
    # Step 4: Optimize directly with xTB
    logging.debug('Step 4: Optimizing directly with xTB (skipping force field)')
    (Path.cwd() / 'scratch').mkdir(exist_ok=True)
    xtb_optimized_mols = list(executor.map(ASE(calculator.XTB()).optimize, unoptimized_unique))
    xtb_optimized_mols = [mol for mol in xtb_optimized_mols if mol is not None]
    logging.debug(f'xTB optimized {len(xtb_optimized_mols)} conformers')
    
    # Compute xTB energies
    energies = list(executor.map(xtb_energy, xtb_optimized_mols))
    for i, conformer in enumerate(unique_conformers):
        if i < len(xtb_optimized_mols):
            conformer.stages[XTB] = xtb_optimized_mols[i]
            conformer.energies[XTB] = energies[i]
    
    # Sort by energy and save
    unique_conformers_with_xtb = [c for c in unique_conformers if XTB in c.stages]
    unique_conformers_with_xtb.sort(key=lambda c: c.energies[XTB])
    
    stk_list_to_xyz_file(
        [c.stages[XTB] for c in unique_conformers_with_xtb],
        'conformers_1_xtb.xyz'
    )
    
    # Log final energies
    logging.info('Final conformer energies (sorted):')
    for i, conformer in enumerate(unique_conformers_with_xtb):
        energy_hartree = conformer.energies[XTB]
        # Convert to kcal/mol relative to lowest
        if i == 0:
            min_energy = energy_hartree
        rel_energy_kcal = (energy_hartree - min_energy) * 627.509  # Hartree to kcal/mol
        logging.info(f'  Conformer {i+1}: {energy_hartree:.6f} Ha ({rel_energy_kcal:.2f} kcal/mol relative)')

logging.debug('Complete! Generated conformer files:')
logging.debug('  - conformers_0_with_boron.xyz (initial with B bridge)')
logging.debug('  - conformers_0_unoptimized.xyz (with H replacement)')
logging.debug('  - conformers_1_xtb.xyz (xTB optimized directly, sorted by energy)')
