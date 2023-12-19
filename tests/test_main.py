from pathlib import Path
import numpy as np
import pytest
import stk
from conformational_sampling.ase_stko_optimizer import ASE
from conformational_sampling.config import Config
from conformational_sampling.main import ConformerEnsembleOptimizer, bind_to_dimethyl_Pd, gen_confs_openbabel, load_stk_mol
import stko
from xtb.ase.calculator import XTB
from conformational_sampling.utils import stk_mol_to_pybel_mol
from conformational_sampling.utils import pybel_mol_to_stk_mol

def test_imports():
    assert True

def test_stk_to_pybel_loop():
    'convert an stk molecule to a pybel molecule and then back to validate the conversion methods'
    butane_initial = stk.BuildingBlock('CCCC')
    butane_final = pybel_mol_to_stk_mol(stk_mol_to_pybel_mol(butane_initial))
    initial_positions = butane_initial.get_position_matrix()
    final_positions = butane_final.get_position_matrix()
    assert np.allclose(initial_positions, final_positions, atol=1e-3)


@pytest.mark.parametrize(
    "ligand_path", # test monodentate and bidentate examples
    [Path('examples/suzuki/example2_L1.xyz'), Path('examples/dppe/ligand.xyz')]
)
def test_metal_ligand_binding(ligand_path):
    stk_ligand = load_stk_mol(ligand_path)

    # specifies atoms of the ligand that bind to the metal, in this case as a smarts string
    functional_group_factory = stk.SmartsFunctionalGroupFactory(
        smarts='P',
        bonders=(0,),
        deleters=(),
    )
    stk_ligand = stk.BuildingBlock.init_from_molecule(
        stk_ligand,
        functional_groups=[functional_group_factory]
    )
    stk_mol = bind_to_dimethyl_Pd(stk_ligand)
    

@pytest.mark.skip(reason='implementation in progress')
def test_conformer_ensemble_optimizer():
    simple_ligand = stk.BuildingBlock('CPCC')
    functional_group_factory = stk.SmartsFunctionalGroupFactory(
        smarts='P',
        bonders=(0,),
        deleters=(),
    )
    simple_ligand = stk.BuildingBlock.init_from_molecule(
        simple_ligand,
        functional_groups=[functional_group_factory]
    )
    simple_complex = bind_to_dimethyl_Pd(simple_ligand)

    config = Config(initial_conformers=2)
    stk_confs = gen_confs_openbabel(simple_complex, config)
    assert len(stk_confs) == 2
    ConformerEnsembleOptimizer(stk_confs, config).optimize()
    