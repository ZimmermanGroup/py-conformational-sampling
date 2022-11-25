#!/export/zimmerman/joshkamm/apps/mambaforge/envs/conformational-sampling/bin/python
#SBATCH -p guest --job-name=conformational_sampling
#SBATCH -c16
#SBATCH -o output.txt

from pathlib import Path
import stk
from conformational_sampling.main import load_stk_mol, gen_ligand_library_entry
from conformational_sampling.config import Config

ligand_path = Path('ligand.sdf') # name of file containing ligand geometry
stk_ligand = load_stk_mol(ligand_path, fmt='sdf')

# specifies atoms of the ligand that bind to the metal, in this case as a smarts string
functional_group_factory = stk.SmartsFunctionalGroupFactory(
    smarts='P',
    bonders=(0,),
    deleters=(),
)
stk_ligand = stk.BuildingBlock.init_from_molecule(stk_ligand, functional_groups=[functional_group_factory])

# py-conformational-sampling configuration object
config = Config(initial_conformers=100, xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb')

# generates conformers, performs multiple step optimization and uniqueness filtering
gen_ligand_library_entry(stk_ligand, config)
