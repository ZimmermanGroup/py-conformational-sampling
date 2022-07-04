# %%
# jupyter only imports
%reload_ext autoreload
%autoreload 2

#!/export/zimmerman/joshkamm/apps/mambaforge/envs/conformational-sampling/bin/python
#SBATCH -p zimintel --job-name=conformational_sampling
#SBATCH -c4
#SBATCH -o output.txt

from pathlib import Path
import stk
from conformational_sampling.main import load_stk_mol, gen_ligand_library_entry

mol_path = Path('ligand.xyz')
stk_ligand = load_stk_mol(mol_path)

# functional groups define which atoms of the ligand bind to the metal
# functional_groups = stk.SmartsFunctionalGroupFactory(
#     smarts='P',
#     bonders=(0,),
#     deleters=(),
# ).get_functional_groups(stk_ligand)
# stk_ligand = stk_ligand.with_functional_groups(functional_groups)

functional_group_factory = stk.SmartsFunctionalGroupFactory(
    smarts='P',
    bonders=(0,),
    deleters=(),
)
stk_ligand = stk.BuildingBlock.init_from_molecule(stk_ligand, functional_groups=[functional_group_factory])

gen_ligand_library_entry(stk_ligand)
