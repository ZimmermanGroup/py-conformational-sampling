from pathlib import Path
import numpy as np

import pytest
import stk
import stko
import os

from conformational_sampling.ase_stko_optimizer import ASE
from conformational_sampling.config import Config
from conformational_sampling.gsm import stk_gsm, stk_se_de_gsm
from conformational_sampling.main import (ConformerEnsembleOptimizer, bind_ligands, bind_to_dimethyl_Pd, load_stk_mol, load_stk_mol_list, stk_list_to_xyz_file, suzuki_ligand_conf_gen)
from conformational_sampling.utils import stk_metal


# stk_ligand_5a = load_stk_mol(Path('example2_5a.xyz'))
# stk_ligand_6a = load_stk_mol(Path('example2_6a.xyz'))
# stk_ancillary_ligand = load_stk_mol(Path('example2_L1.xyz'))

# # specify functional groups of each ligand that bind to the metal, 
# # in this case as a smarts string
# functional_group_factory = stk.SmartsFunctionalGroupFactory(
#     smarts='cBr',
#     bonders=(0,),
#     deleters=(1,),
# )
# stk_ligand_5a = stk.BuildingBlock.init_from_molecule(stk_ligand_5a, functional_groups=[functional_group_factory])
# functional_group_factory = stk.SmartsFunctionalGroupFactory(
#     smarts='cB([O][H])[O][H]',
#     bonders=(0,),
#     deleters=tuple(range(1,6)),
# )
# stk_ligand_6a = stk.BuildingBlock.init_from_molecule(stk_ligand_6a, functional_groups=[functional_group_factory])
# functional_group_factory = stk.SmartsFunctionalGroupFactory(
#     smarts='P',
#     bonders=(0,),
#     deleters=(),
# )
# stk_ancillary_ligand = stk.BuildingBlock.init_from_molecule(stk_ancillary_ligand, functional_groups=[functional_group_factory])


config = Config(
    initial_conformers=2,
    num_cpus=2,
)
# suzuki_ligand_conf_gen(stk_ligand_5a, stk_ligand_6a, stk_ancillary_ligand, config)

conformer_path = Path('suzuki_conformers.xyz')
conformer_mols = load_stk_mol_list(conformer_path)

# for i in range(len(conformer_mols)):
for i in range(1):
    path = Path.cwd() / f'scratch/pystring_{i}'
    path.mkdir(parents=True, exist_ok=True)

os.chdir(path)

driving_coordinates = [('ADD',56,80),('BREAK',1,56),('BREAK',1,80)]

job_index=0
stk_se_de_gsm(
    stk_mol=conformer_mols[job_index],
    driving_coordinates=driving_coordinates,
    config=config,
)
# stk_gsm(
#     stk_mol=optimized_stk_mol,
#     driving_coordinates=driving_coordinates,
#     config=config,
# )
assert True

# if __name__ == "__main__":
#     suzuki()