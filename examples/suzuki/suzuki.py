import logging
from pathlib import Path

FORMAT = (
    '[%(asctime)s %(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s'
)
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

import stk

from conformational_sampling.catalytic_reaction_complex import (
    ReductiveEliminationComplex,
)
from conformational_sampling.config import Config
from conformational_sampling.gsm import stk_se_de_gsm_single_node_parallel
from conformational_sampling.main import load_stk_mol, load_stk_mol_list
from conformational_sampling.utils import stk_metal

# specify file paths and functional groups of each ligand that bind to the metal
reactive_ligand_1 = stk.BuildingBlock.init_from_molecule(
    molecule=load_stk_mol(Path('example2_5a.xyz')),
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='cBr',
            bonders=(0,),
            deleters=(1,),
        )
    ],
)
reactive_ligand_2 = stk.BuildingBlock.init_from_molecule(
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

# py-conformational-sampling configuration object
# (defaults are small for testing purposes)
config = Config(
    # larger value samples more thoroughly, set to roughly 20-100 for thorough
    # sampling depending on ligand flexibility and resource availability
    initial_conformers=3,
    num_cpus=3,  # set to number of cpus available for use by pycosa
)

reactive_complex = ReductiveEliminationComplex(
    metal=stk_metal('Pd'),
    ancillary_ligand=ancillary_ligand,
    reactive_ligand_1=reactive_ligand_1,
    reactive_ligand_2=reactive_ligand_2,
    config=config,
)

# set to True to start visualization after conformer reaction paths have already
# been generated
start_visualization = False
if start_visualization:
    import panel as pn

    from conformational_sampling.visualization import ConformationalSamplingDashboard

    pn.serve(
        ConformationalSamplingDashboard(reactive_complex).app(), show=False, port=5006
    )
    raise SystemExit

# generates conformers including multi-phase optimization and uniqueness filtering
reactive_complex.gen_conformers()

conformer_path = Path('conformers_3_xtb.xyz')
conformer_mols = load_stk_mol_list(conformer_path)
# for testing, limit number of pyGSM reaction path optimizations to number of cpus
# remove or comment out this line for a full calculation
conformer_mols = conformer_mols[: config.num_cpus]

driving_coordinates = reactive_complex.gen_reductive_elim_drive_coords()

stk_se_de_gsm_single_node_parallel(
    stk_mols=conformer_mols,
    driving_coordinates=driving_coordinates,
    config=config,
)
