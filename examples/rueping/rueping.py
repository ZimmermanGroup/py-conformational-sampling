import logging
from pathlib import Path

FORMAT = (
    '[%(asctime)s %(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s'
)
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

from conformational_sampling.main import load_stk_mol, gen_confs_openbabel, stk_list_to_xyz_file
from conformational_sampling.config import Config

# Load the full molecular system
full_system_path = Path('full_system.xyz')
stk_mol = load_stk_mol(full_system_path)

# py-conformational-sampling configuration object
# (defaults are small for testing purposes)
config = Config(
    # larger value samples more thoroughly, set to roughly 20-100 for thorough
    # sampling depending on ligand flexibility and resource availability
    initial_conformers=3,
    num_cpus=3,  # set to number of cpus available for use by pycosa
)

# Generate conformers using OpenBabel
logging.debug('Start generating conformers')
stk_conformers = gen_confs_openbabel(stk_mol, config)
logging.debug(f'Generated {len(stk_conformers)} conformers')

# Save conformers to file
output_file = 'conformers_openbabel.xyz'
stk_list_to_xyz_file(stk_conformers, output_file)
logging.debug(f'Saved conformers to {output_file}')
