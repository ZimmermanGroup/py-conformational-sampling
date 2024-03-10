import os
from pathlib import Path

import stk
from conformational_sampling.catalytic_reaction_complex import CatalyticReactionComplex
from conformational_sampling.config import Config
from conformational_sampling.gsm import stk_de_gsm, stk_se_gsm
from conformational_sampling.main import load_stk_mol, load_stk_mol_list
from conformational_sampling.utils import stk_metal

stk_ligand_5a = load_stk_mol(Path('example2_5a.xyz'))
stk_ligand_6a = load_stk_mol(Path('example2_6a.xyz'))
stk_ancillary_ligand = load_stk_mol(Path('example2_L1.xyz'))

# specify functional groups of each ligand that bind to the metal, 
# in this case as smarts strings
functional_group_factory = stk.SmartsFunctionalGroupFactory(
    smarts='cBr',
    bonders=(0,),
    deleters=(1,),
)
stk_ligand_5a = stk.BuildingBlock.init_from_molecule(stk_ligand_5a, functional_groups=[functional_group_factory])
functional_group_factory = stk.SmartsFunctionalGroupFactory(
    smarts='cB([O][H])[O][H]',
    bonders=(0,),
    deleters=tuple(range(1,6)),
)
stk_ligand_6a = stk.BuildingBlock.init_from_molecule(stk_ligand_6a, functional_groups=[functional_group_factory])
functional_group_factory = stk.SmartsFunctionalGroupFactory(
    smarts='P',
    bonders=(0,),
    deleters=(),
)
stk_ancillary_ligand = stk.BuildingBlock.init_from_molecule(stk_ancillary_ligand, functional_groups=[functional_group_factory])

config = Config(
    initial_conformers=2,
    num_cpus=2,
)
reactive_complex = CatalyticReactionComplex(
    metal=stk_metal('Pd'),
    ancillary_ligand=stk_ancillary_ligand,
    reactive_ligand_1=stk_ligand_5a,
    reactive_ligand_2=stk_ligand_6a,
    config=config,
)
reactive_complex.gen_conformers()

conformer_path = Path('conformers_3_xtb.xyz')
conformer_mols = load_stk_mol_list(conformer_path)

# for i in range(len(conformer_mols)):
for i in range(1):
    path = Path.cwd() / f'scratch/pystring_{i}'
    path.mkdir(parents=True, exist_ok=True)

os.chdir(path)

# specifying driving coordinates for pyGSM (1-indexed based on XYZ file)
# driving_coordinates = [('ADD',56,80),('BREAK',1,56),('BREAK',1,80)]
driving_coordinates = reactive_complex.gen_reductive_elim_drive_coords()

job_index=0
stk_se_gsm(
    stk_mol=conformer_mols[job_index],
    driving_coordinates=driving_coordinates,
    config=config,
)
stk_de_gsm(config=config)
assert True
