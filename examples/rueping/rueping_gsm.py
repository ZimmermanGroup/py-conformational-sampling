"""
GSM (Growing String Method) script for Rueping catalyst system.

This runs SE-GSM (Single-Ended GSM) on the xTB-optimized conformers to find
reaction pathways on the xTB potential energy surface.

The driving coordinates model a cationic 2-aza-Cope rearrangement, where a C-C bond
forms while another C-C bond breaks, representing a key step in the reaction mechanism.
"""
import logging
from pathlib import Path

FORMAT = (
    '[%(asctime)s %(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s'
)
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

from conformational_sampling.config import Config
from conformational_sampling.gsm import stk_se_de_gsm_single_node_parallel
from conformational_sampling.main import load_stk_mol_list

# Load the xTB-optimized conformers from the previous step
conformer_path = Path('conformers_1_xtb.xyz')
conformer_mols = load_stk_mol_list(conformer_path)
logging.info(f'Loaded {len(conformer_mols)} conformers from {conformer_path}')

# Configuration
config = Config(
    num_cpus=3,  # set to number of cpus available for parallel GSM runs
)

# Define driving coordinates for the cationic 2-aza-Cope rearrangement
# NOTE: Atom indices are 1-indexed (as in Molden), not 0-indexed
# This matches the chemical convention and pyGSM's expected format
# 
# Cationic 2-aza-Cope rearrangement mechanism:
#   ADD:   C(100)-C(86)  [Form new C-C bond]
#   BREAK: C(93)-C(94)   [Break C-C bond]
#
# This [3,3]-sigmatropic rearrangement is a key step in the Rueping catalyst mechanism
driving_coordinates = [
    ('ADD', 100, 86),
    ('BREAK', 93, 94),
]

logging.info('Driving coordinates for cationic 2-aza-Cope rearrangement:')
logging.info('  ADD:   C(100)-C(86)  [new C-C bond formation]')
logging.info('  BREAK: C(93)-C(94)   [C-C bond cleavage]')

# Optional: For testing, limit to fewer conformers
# conformer_mols = conformer_mols[:config.num_cpus]

logging.info(f'Running SE-GSM on {len(conformer_mols)} conformers')
logging.info(f'Driving coordinates: {driving_coordinates}')

# Run SE-GSM followed by DE-GSM in parallel for all conformers
stk_se_de_gsm_single_node_parallel(
    stk_mols=conformer_mols,
    driving_coordinates=driving_coordinates,
    config=config,
)

logging.info('GSM calculations complete!')
logging.info('Results are in scratch/pystring_*/ directories')
