"""
Check if OpenBabel properly recognizes the B atom bonding to both O and N atoms
in the full_system.xyz file, and add bonds manually if necessary.
"""
import logging
from pathlib import Path
from openbabel import pybel as pb
from openbabel import openbabel as ob

FORMAT = (
    '[%(asctime)s %(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s'
)
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

def check_and_fix_boron_bonding():
    """
    Load the full_system.xyz file, check if B is bonded to both O and N,
    and add the necessary bonds if missing.
    """
    # Load the molecule
    full_system_path = Path('full_system.xyz')
    logging.info(f'Loading molecule from {full_system_path}')
    
    pybel_mol = next(pb.readfile('xyz', str(full_system_path)))
    obmol = pybel_mol.OBMol
    
    # Find the boron atom
    boron_atom = None
    boron_idx = None
    for atom in ob.OBMolAtomIter(obmol):
        if atom.GetAtomicNum() == 5:  # B has atomic number 5
            boron_atom = atom
            boron_idx = atom.GetIdx()
            logging.info(f'Found Boron atom at index {boron_idx}')
            break
    
    if boron_atom is None:
        logging.error('No Boron atom found in the molecule!')
        return None
    
    # Get coordinates of boron
    b_x, b_y, b_z = boron_atom.GetX(), boron_atom.GetY(), boron_atom.GetZ()
    logging.info(f'Boron coordinates: ({b_x:.4f}, {b_y:.4f}, {b_z:.4f})')
    
    # Check current bonds to boron
    logging.info('Current bonds to Boron:')
    bonded_atoms = []
    for bond in ob.OBAtomBondIter(boron_atom):
        neighbor = bond.GetNbrAtom(boron_atom)
        bonded_atoms.append(neighbor)
        element = ob.GetSymbol(neighbor.GetAtomicNum())
        logging.info(f'  - Bonded to {element} at index {neighbor.GetIdx()}')
    
    # Find nearby O and N atoms (within reasonable bonding distance)
    # Typical B-O bond: ~1.4 Å, B-N bond: ~1.4-1.6 Å
    MAX_BOND_DISTANCE = 2.0  # Angstroms, generous threshold
    
    oxygen_candidates = []
    nitrogen_candidates = []
    
    for atom in ob.OBMolAtomIter(obmol):
        if atom.GetIdx() == boron_idx:
            continue
        
        # Calculate distance to boron
        dx = atom.GetX() - b_x
        dy = atom.GetY() - b_y
        dz = atom.GetZ() - b_z
        distance = (dx*dx + dy*dy + dz*dz) ** 0.5
        
        if distance < MAX_BOND_DISTANCE:
            if atom.GetAtomicNum() == 8:  # Oxygen
                oxygen_candidates.append((atom, distance))
            elif atom.GetAtomicNum() == 7:  # Nitrogen
                nitrogen_candidates.append((atom, distance))
    
    # Sort by distance
    oxygen_candidates.sort(key=lambda x: x[1])
    nitrogen_candidates.sort(key=lambda x: x[1])
    
    logging.info(f'\nNearby Oxygen atoms (within {MAX_BOND_DISTANCE} Å):')
    for atom, dist in oxygen_candidates:
        logging.info(f'  - O at index {atom.GetIdx()}, distance: {dist:.4f} Å')
    
    logging.info(f'\nNearby Nitrogen atoms (within {MAX_BOND_DISTANCE} Å):')
    for atom, dist in nitrogen_candidates:
        logging.info(f'  - N at index {atom.GetIdx()}, distance: {dist:.4f} Å')
    
    # Check if we need to add bonds
    bonded_to_oxygen = any(atom.GetAtomicNum() == 8 for atom in bonded_atoms)
    bonded_to_nitrogen = any(atom.GetAtomicNum() == 7 for atom in bonded_atoms)
    
    logging.info(f'\nBoron bonding status:')
    logging.info(f'  - Bonded to Oxygen: {bonded_to_oxygen}')
    logging.info(f'  - Bonded to Nitrogen: {bonded_to_nitrogen}')
    
    bonds_added = False
    
    # Add bond to closest oxygen if not already bonded
    if not bonded_to_oxygen and oxygen_candidates:
        closest_o, dist = oxygen_candidates[0]
        logging.info(f'\nAdding bond between B (idx {boron_idx}) and O (idx {closest_o.GetIdx()}), distance: {dist:.4f} Å')
        obmol.AddBond(boron_idx, closest_o.GetIdx(), 1)  # 1 = single bond
        bonds_added = True
    
    # Add bond to closest nitrogen if not already bonded
    if not bonded_to_nitrogen and nitrogen_candidates:
        closest_n, dist = nitrogen_candidates[0]
        logging.info(f'Adding bond between B (idx {boron_idx}) and N (idx {closest_n.GetIdx()}), distance: {dist:.4f} Å')
        obmol.AddBond(boron_idx, closest_n.GetIdx(), 1)  # 1 = single bond
        bonds_added = True
    
    if bonds_added:
        # Save the modified molecule
        output_path = 'full_system_with_bonds.xyz'
        logging.info(f'\nSaving modified molecule to {output_path}')
        pybel_mol.write('xyz', output_path, overwrite=True)
        
        # Also save as mol file to preserve bond information
        mol_output_path = 'full_system_with_bonds.mol'
        logging.info(f'Saving with bond information to {mol_output_path}')
        pybel_mol.write('mol', mol_output_path, overwrite=True)
        
        # Verify the bonds were added
        logging.info('\nVerifying bonds after modification:')
        for bond in ob.OBAtomBondIter(boron_atom):
            neighbor = bond.GetNbrAtom(boron_atom)
            element = ob.GetSymbol(neighbor.GetAtomicNum())
            logging.info(f'  - Bonded to {element} at index {neighbor.GetIdx()}')
        
        return output_path
    else:
        logging.info('\nNo bonds needed to be added - B is already properly bonded!')
        return None

if __name__ == '__main__':
    result = check_and_fix_boron_bonding()
    if result:
        print(f'\n✓ Modified molecule saved to: {result}')
    else:
        print('\n✓ Molecule already has correct bonding')
