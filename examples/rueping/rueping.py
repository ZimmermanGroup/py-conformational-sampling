import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

import stko
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

FORMAT = (
    '[%(asctime)s %(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s'
)
logging.basicConfig(format=FORMAT, level=logging.DEBUG)


def replace_boron_with_hydrogen(stk_mol):
    """Replace boron atom with hydrogen in an stk molecule.
    
    Simply changes the atomic number from B (5) to H (1) and writes to XYZ format.
    When reading back from XYZ, OpenBabel will reperceive bonds based on distances,
    and hydrogen at that position will only form reasonable bonds.
    """
    # Convert to pybel
    pybel_mol = stk_mol_to_pybel_mol(stk_mol)
    obmol = pybel_mol.OBMol
    
    # Find boron atom and change it to hydrogen
    for atom in pb.ob.OBMolAtomIter(obmol):
        if atom.GetAtomicNum() == 5:  # Boron
            atom.SetAtomicNum(1)  # Change to Hydrogen
            logging.debug(f'Changed B at index {atom.GetIdx()} to H')
            break
    
    # Write to XYZ and read back to reperceive bonds
    xyz_string = pybel_mol.write('xyz')
    pybel_mol_new = pb.readstring('xyz', xyz_string)
    return pybel_mol_to_stk_mol(pybel_mol_new)


def get_unique_conformer_ids(conformers, stage, rms_threshold=0.175):
    """Get indices of unique conformers based on RMS comparison."""
    rdkit_mols = {i: Chem.RemoveHs(conformer.stages[stage].to_rdkit_mol())
                  for i, conformer in enumerate(conformers)
                  if stage in conformer.stages}
    unique_indices = []
    for i, rdkit_mol in rdkit_mols.items():
        for j in unique_indices:
            if Chem.AllChem.CalcRMS(rdkit_mol, rdkit_mols[j]) < rms_threshold:
                break
        else:
            unique_indices.append(i)
    return unique_indices


# Load the full molecular system with boron bridge
full_system_path = Path('full_system.xyz')
stk_mol = load_stk_mol(full_system_path)

# Configuration
config = Config(
    # larger value samples more thoroughly, set to roughly 20-100 for thorough
    # sampling depending on ligand flexibility and resource availability
    initial_conformers=3,
    num_cpus=3,  # set to number of cpus available for use by pycosa
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

# Step 3: Optimize with standard force field (skip MCHammer since no fragments, skip MetalOptimizer since no metal)
logging.debug('Step 3: Optimizing with standard force field (UFF)')
with ProcessPoolExecutor(max_workers=config.num_cpus) as executor:
    unoptimized_mols = [conformer.stages[UNOPTIMIZED] for conformer in conformers]
    
    # Use stko.UFF for force field optimization
    ff_optimized_mols = list(executor.map(stko.UFF().optimize, unoptimized_mols))
    logging.debug(f'Force field optimized {len(ff_optimized_mols)} conformers')
    
    for i, conformer in enumerate(conformers):
        conformer.stages[METAL_OPTIMIZER] = ff_optimized_mols[i]
    
    # Save force field optimized
    stk_list_to_xyz_file(ff_optimized_mols, 'conformers_1_force_field.xyz')
    
    # Step 4: Remove duplicates before xTB
    logging.debug('Step 4: Removing duplicate conformers')
    unique_ids = get_unique_conformer_ids(conformers, METAL_OPTIMIZER, config.pre_xtb_rms_threshold)
    unique_conformers = [conformers[i] for i in unique_ids]
    logging.debug(f'Kept {len(unique_conformers)} unique conformers (removed {len(conformers) - len(unique_conformers)} duplicates)')
    
    ff_optimized_unique = [conformer.stages[METAL_OPTIMIZER] for conformer in unique_conformers]
    
    # Step 5: Optimize with xTB
    logging.debug('Step 5: Optimizing with xTB')
    (Path.cwd() / 'scratch').mkdir(exist_ok=True)
    xtb_optimized_mols = list(executor.map(ASE(calculator.XTB()).optimize, ff_optimized_unique))
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
        'conformers_2_xtb.xyz'
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
logging.debug('  - conformers_1_force_field.xyz (UFF optimized)')
logging.debug('  - conformers_2_xtb.xyz (xTB optimized, sorted by energy)')
