"""
Verify that the B atom maintains bonds to both O and N in all generated conformers.
"""
import logging
from pathlib import Path
from openbabel import pybel as pb
from openbabel import openbabel as ob

FORMAT = (
    '[%(asctime)s %(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s'
)
logging.basicConfig(format=FORMAT, level=logging.INFO)

def verify_conformers():
    """
    Check all conformers in conformers_openbabel.xyz to verify B bonding.
    """
    conformers_path = Path('conformers_openbabel.xyz')
    
    if not conformers_path.exists():
        logging.error(f'{conformers_path} does not exist!')
        return
    
    logging.info(f'Loading conformers from {conformers_path}')
    conformers = list(pb.readfile('xyz', str(conformers_path)))
    logging.info(f'Found {len(conformers)} conformers\n')
    
    all_good = True
    
    for i, pybel_mol in enumerate(conformers, 1):
        logging.info(f'=== Conformer {i} ===')
        obmol = pybel_mol.OBMol
        
        # Find the boron atom
        boron_atom = None
        boron_idx = None
        for atom in ob.OBMolAtomIter(obmol):
            if atom.GetAtomicNum() == 5:  # B has atomic number 5
                boron_atom = atom
                boron_idx = atom.GetIdx()
                break
        
        if boron_atom is None:
            logging.error(f'  ✗ No Boron atom found!')
            all_good = False
            continue
        
        # Get boron coordinates
        b_x, b_y, b_z = boron_atom.GetX(), boron_atom.GetY(), boron_atom.GetZ()
        
        # Check bonds to boron
        bonded_atoms = []
        bonded_to_oxygen = False
        bonded_to_nitrogen = False
        
        for bond in ob.OBAtomBondIter(boron_atom):
            neighbor = bond.GetNbrAtom(boron_atom)
            bonded_atoms.append(neighbor)
            element = ob.GetSymbol(neighbor.GetAtomicNum())
            
            # Calculate distance
            dx = neighbor.GetX() - b_x
            dy = neighbor.GetY() - b_y
            dz = neighbor.GetZ() - b_z
            distance = (dx*dx + dy*dy + dz*dz) ** 0.5
            
            logging.info(f'  - B bonded to {element} at idx {neighbor.GetIdx()}, distance: {distance:.4f} Å')
            
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                bonded_to_oxygen = True
            elif neighbor.GetAtomicNum() == 7:  # Nitrogen
                bonded_to_nitrogen = True
        
        # Check if both O and N bonds are present
        if bonded_to_oxygen and bonded_to_nitrogen:
            logging.info('  ✓ B is properly bonded to both O and N')
        else:
            logging.warning(f'  ✗ B bonding issue: O={bonded_to_oxygen}, N={bonded_to_nitrogen}')
            all_good = False
        
        logging.info('')
    
    if all_good:
        print('\n✓✓✓ All conformers maintain proper B-O and B-N bonds! ✓✓✓')
    else:
        print('\n✗✗✗ Some conformers have bonding issues! ✗✗✗')
    
    return all_good

if __name__ == '__main__':
    verify_conformers()
