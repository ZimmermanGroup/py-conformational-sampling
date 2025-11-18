import logging
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation

FORMAT = (
    '[%(asctime)s %(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s'
)
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

import stk
from conformational_sampling.main import load_stk_mol, gen_confs_openbabel, stk_list_to_xyz_file
from conformational_sampling.config import Config

# =============================================================================
# Phase 1: Load fragment structures and identify H-bond atoms
# =============================================================================

catalyst_path = Path('RUEPING-CATALYST.xyz')
substrate_path = Path('RUEPING-SUBSTRATE.xyz')

logging.debug('Phase 1: Loading catalyst structure')
catalyst = load_stk_mol(catalyst_path)
logging.debug(f'Catalyst loaded: {catalyst.get_num_atoms()} atoms')

logging.debug('Phase 1: Loading substrate structure')
substrate = load_stk_mol(substrate_path)
logging.debug(f'Substrate loaded: {substrate.get_num_atoms()} atoms')

# Verify total atom count matches full system (78 + 51 = 129)
total_atoms = catalyst.get_num_atoms() + substrate.get_num_atoms()
logging.debug(f'Total atoms in fragments: {total_atoms}')

# Find H-bond atoms using SMARTS patterns
def find_atoms_by_smarts(stk_mol, smarts_pattern, description):
    """Find atoms matching a SMARTS pattern in an stk molecule"""
    rdkit_mol = stk_mol.to_rdkit_mol()
    from rdkit import Chem
    pattern = Chem.MolFromSmarts(smarts_pattern)
    matches = rdkit_mol.GetSubstructMatches(pattern)
    logging.debug(f'{description}: SMARTS "{smarts_pattern}" found {len(matches)} match(es)')
    if len(matches) == 0:
        raise ValueError(f'No matches found for SMARTS pattern "{smarts_pattern}" in {description}')
    return matches[0]  # Return first match as tuple of atom indices

# Find H-bond donor in catalyst: O-H group
# SMARTS: [O][H] matches oxygen bonded to hydrogen (matches both atoms)
catalyst_oh_atoms = find_atoms_by_smarts(catalyst, '[O][H]', 'Catalyst O-H')
donor_o_idx = catalyst_oh_atoms[0]
donor_h_idx = catalyst_oh_atoms[1]
logging.debug(f'Catalyst H-bond donor: O={donor_o_idx}, H={donor_h_idx}')

# Find H-bond acceptor in substrate: nitrogen
# SMARTS: [N] matches any nitrogen
substrate_n_atoms = find_atoms_by_smarts(substrate, '[N]', 'Substrate N')
acceptor_n_idx = substrate_n_atoms[0]
logging.debug(f'Substrate H-bond acceptor N: atom index {acceptor_n_idx}')

logging.debug('Phase 1 complete: Fragments loaded and H-bond atoms identified')

# =============================================================================
# Phase 2: Generate conformers for each fragment independently
# =============================================================================

config = Config(
    initial_conformers=3,  # Small number for testing; increase for production
    num_cpus=3,
)

logging.debug('Phase 2: Starting conformer generation for catalyst')
catalyst_conformers = gen_confs_openbabel(catalyst, config)
logging.debug(f'Generated {len(catalyst_conformers)} catalyst conformers')
stk_list_to_xyz_file(catalyst_conformers, 'catalyst_conformers.xyz')

logging.debug('Phase 2: Starting conformer generation for substrate')
substrate_conformers = gen_confs_openbabel(substrate, config)
logging.debug(f'Generated {len(substrate_conformers)} substrate conformers')
stk_list_to_xyz_file(substrate_conformers, 'substrate_conformers.xyz')

logging.debug('Phase 2 complete: Conformers generated and saved')

# =============================================================================
# Phase 4: Create combinatorial assembly function
# =============================================================================

def find_neighbor_heavy_atom(rdkit_mol, atom_idx):
    """Find the first heavy (non-hydrogen) atom neighbor"""
    atom = rdkit_mol.GetAtomWithIdx(atom_idx)
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() != 'H':
            return neighbor.GetIdx()
    raise ValueError(f'No heavy atom neighbor found for atom {atom_idx}')

def combine_fragments_with_hbond(catalyst_stk, substrate_stk, 
                                  cat_o_idx, cat_h_idx, sub_n_idx,
                                  dihedral_angle_deg=0.0):
    """
    Combine catalyst and substrate conformers with a hydrogen bond.
    
    Uses stk to handle initial positioning by treating fragments as building blocks,
    then applies dihedral rotation using RDKit.
    
    The dihedral is defined by: [C-O-H...N] where C is a heavy neighbor of O
    
    Parameters:
    -----------
    catalyst_stk : stk.Molecule
        Catalyst conformer
    substrate_stk : stk.Molecule
        Substrate conformer
    cat_o_idx : int
        Index of oxygen atom in catalyst (H-bond donor)
    cat_h_idx : int
        Index of hydrogen atom in catalyst
    sub_n_idx : int
        Index of nitrogen atom in substrate (H-bond acceptor)
    dihedral_angle_deg : float
        Dihedral angle C-O-H...N in degrees
    
    Returns:
    --------
    stk.Molecule
        Combined structure with hydrogen bond
    """
    
    # Convert fragments to BuildingBlocks with functional groups for stk assembly
    # The functional groups define where the fragments will "connect" (near H-bond)
    
    catalyst_bb = stk.BuildingBlock.init_from_molecule(
        molecule=catalyst_stk,
        functional_groups=[
            # stk.SmartsFunctionalGroupFactory(
            #     smarts='[P][O][H]',
            #     bonders=(2,),  # H atom
            #     deleters=(),
            # )  # Functional group at OH
            stk.SingleAtom(stk.O(cat_o_idx))  # Functional group at O

        ]
    )
    
    substrate_bb = stk.BuildingBlock.init_from_molecule(
        molecule=substrate_stk,
        functional_groups=[
            # stk.SingleAtom(stk.N(sub_n_idx))  # Functional group at N
            stk.SmartsFunctionalGroupFactory(
                smarts='N(C)=C',
                bonders=(0,),  # N atom
                deleters=(),
                # placers=(1,2),  # C atoms bonded to N
            )  # Functional group at N
        ]
    )
    
    # Use stk's Linear polymer topology to position the two fragments
    # This will handle translation and rotation to align them
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(catalyst_bb, substrate_bb),
            repeating_unit='AB',
            num_repeating_units=1,
            optimizer=stk.Collapser(scale_steps=False)  # Bring fragments close together
        )
        # topology_graph=stk.small.NCore(
        #     core_building_block=catalyst_bb,
        #     arm_building_blocks=substrate_bb,
        #     repeating_unit='A',
        #     # optimizer=stk.Collapser(scale_steps=False)  # Bring fragments close together
        # )
    )
    
    # Now we have the fragments positioned near each other by stk
    # Apply dihedral rotation if needed
    if abs(dihedral_angle_deg) > 1e-6:
        from rdkit import Chem
        from rdkit.Chem import rdMolTransforms
        
        # Convert to RDKit for dihedral manipulation
        rdkit_mol = polymer.to_rdkit_mol()
        
        # Initialize ring info (required for RDKit operations)
        Chem.FastFindRings(rdkit_mol)
        
        # Find heavy atom neighbor of O for dihedral definition
        cat_c_idx = find_neighbor_heavy_atom(rdkit_mol, cat_o_idx)
        
        # In the combined molecule, substrate atoms are offset
        # We need to find where the N atom ended up
        # stk should preserve atom ordering, so substrate atoms start after catalyst atoms
        num_cat_atoms = catalyst_stk.get_num_atoms()
        sub_n_idx_in_combined = num_cat_atoms + sub_n_idx
        
        # Set the dihedral angle C-O-H-N
        rdMolTransforms.SetDihedralDeg(
            rdkit_mol.GetConformer(),
            cat_c_idx, cat_o_idx, cat_h_idx, sub_n_idx_in_combined,
            dihedral_angle_deg
        )
        
        # Convert back to stk
        combined = stk.BuildingBlock.init_from_rdkit_mol(rdkit_mol)
    else:
        combined = polymer
    
    return combined

logging.debug('Phase 4 complete: Assembly function created')

# =============================================================================
# Phase 5: Generate full combinatorial ensemble
# =============================================================================

# Test the assembly function with a single combination first
logging.debug('Phase 5: Testing assembly function with first conformers')

test_combined = combine_fragments_with_hbond(
    catalyst_conformers[0],
    substrate_conformers[0],
    donor_o_idx,
    donor_h_idx,
    acceptor_n_idx,
    dihedral_angle_deg=0.0
)

logging.debug(f'Test assembly: combined molecule has {test_combined.get_num_atoms()} atoms')
stk_list_to_xyz_file([test_combined], 'test_combined.xyz')

# Generate all combinatorial structures
dihedral_angles = [0, 60, 120, 180, 240, 300]  # degrees
logging.debug(f'Generating combinatorial structures with {len(dihedral_angles)} dihedral angles')

combinatorial_structures = []
for cat_idx, cat_conf in enumerate(catalyst_conformers):
    for sub_idx, sub_conf in enumerate(substrate_conformers):
        for angle in dihedral_angles:
            combined = combine_fragments_with_hbond(
                cat_conf,
                sub_conf,
                donor_o_idx,
                donor_h_idx,
                acceptor_n_idx,
                dihedral_angle_deg=angle
            )
            combinatorial_structures.append(combined)
            logging.debug(f'Generated: catalyst {cat_idx}, substrate {sub_idx}, angle {angle}°')

total_structures = len(combinatorial_structures)
expected = len(catalyst_conformers) * len(substrate_conformers) * len(dihedral_angles)
logging.debug(f'Generated {total_structures} combinatorial structures (expected {expected})')

# Save all combinatorial structures
stk_list_to_xyz_file(combinatorial_structures, 'combinatorial_conformers.xyz')
logging.debug('Phase 5 complete: Combinatorial ensemble generated and saved')

# TODO: Phase 6 - Integration with optimization pipeline
