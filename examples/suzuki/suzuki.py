from pathlib import Path

import stk
from conformational_sampling.catalytic_reaction_complex import CatalyticReactionComplex
from conformational_sampling.config import Config
from conformational_sampling.gsm import stk_se_de_gsm_single_node_parallel
from conformational_sampling.main import load_stk_mol, load_stk_mol_list
from conformational_sampling.utils import stk_metal

# specify file paths and functional groups of each ligand that bind to the metal
ligand_5a = stk.BuildingBlock.init_from_molecule(
    molecule=load_stk_mol(Path('example2_5a.xyz')),
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='cBr',
            bonders=(0,),
            deleters=(1,),
        )
    ],
)
ligand_6a = stk.BuildingBlock.init_from_molecule(
    molecule=load_stk_mol(Path('example2_6a.xyz')),
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='cB([O][H])[O][H]',
            bonders=(0,),
            deleters=tuple(range(1, 6)),
        )
    ],
)
ancillary_ligand = stk.BuildingBlock.init_from_molecule(
    molecule=load_stk_mol(Path('example2_L1.xyz')),
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='P',
            bonders=(0,),
            deleters=(),
        )
    ],
)

config = Config(
    initial_conformers=3,
    num_cpus=3,
)
reactive_complex = CatalyticReactionComplex(
    metal=stk_metal('Pd'),
    ancillary_ligand=ancillary_ligand,
    reactive_ligand_1=ligand_5a,
    reactive_ligand_2=ligand_6a,
    config=config,
)
reactive_complex.gen_conformers()

conformer_path = Path('conformers_3_xtb.xyz')
conformer_mols = load_stk_mol_list(conformer_path)

driving_coordinates = reactive_complex.gen_reductive_elim_drive_coords()

stk_se_de_gsm_single_node_parallel(
    stk_mols=conformer_mols[:3],
    driving_coordinates=driving_coordinates,
    config=config,
)
