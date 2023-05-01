from pathlib import Path
import numpy as np

import pytest
import stk
import stko
from xtb.ase.calculator import XTB

from conformational_sampling.ase_stko_optimizer import ASE
from conformational_sampling.config import Config
from conformational_sampling.gsm import stk_gsm
from conformational_sampling.main import (ConformerEnsembleOptimizer, bind_ligands, bind_to_dimethyl_Pd, load_stk_mol, stk_list_to_xyz_file)
from conformational_sampling.utils import stk_metal


@pytest.mark.xfail(raises=RuntimeError)
def test_gsm():
    # create dppe bound complex as test stk molecule
    ligand_path = Path('examples/dppe/ligand.xyz') # name of file containing ligand geometry
    stk_ligand = load_stk_mol(ligand_path)

    # specifies atoms of the ligand that bind to the metal, in this case as a smarts string
    functional_group_factory = stk.SmartsFunctionalGroupFactory(
        smarts='P',
        bonders=(0,),
        deleters=(),
    )
    stk_ligand = stk.BuildingBlock.init_from_molecule(stk_ligand, functional_groups=[functional_group_factory])
    stk_mol = bind_to_dimethyl_Pd(stk_ligand)
    optimizer_sequence = stko.OptimizerSequence(stk.MCHammer(), stko.MetalOptimizer(), ASE(XTB()))
    stk_mol = optimizer_sequence.optimize(stk_mol)
    
    # run gsm to eliminate ethane
    # driving coordinates are 1-indexed
    driving_coordinates = [['BREAK',1,54],['BREAK',1,58],['ADD',54,58]]
    config = Config(
        xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb',
        ase_calculator=XTB(),
    )
    
    stk_gsm(
        stk_mol=stk_mol,
        driving_coordinates=driving_coordinates,
        config=config,
    )


def test_suzuki():
    stk_ancillary_ligand = load_stk_mol(Path('tests/test_data/example2_L1.xyz'))
    stk_ligand_5a = load_stk_mol(Path('tests/test_data/example2_5a.xyz'))
    stk_ligand_6a = load_stk_mol(Path('tests/test_data/example2_6a.xyz'))

    # specifies atoms of the ancillary ligand that bind to the metal, in this case as a smarts string
    functional_group_factory = stk.SmartsFunctionalGroupFactory(
        smarts='P',
        bonders=(0,),
        deleters=(),
    )
    stk_ancillary_ligand = stk.BuildingBlock.init_from_molecule(stk_ancillary_ligand, functional_groups=[functional_group_factory])
    # stk_mol = bind_to_dimethyl_Pd(stk_ancillary_ligand)

    # prepare functional groups for the reactive ligands
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
    stk_mol = bind_ligands(stk_metal('Pd'), stk_ancillary_ligand, stk_ligand_5a, stk_ligand_6a)
    stk_list_to_xyz_file([stk_mol], 'test_Pd_complex.xyz')
    
    config = Config(
        xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb',
        ase_calculator=XTB(),
    )
    # ConformerEnsembleOptimizer([stk_mol], config).optimize()

    # driving coordinates are 1-indexed
    # JOSH - look up if I already wrote an atom map from constructed molecule to building block
    # JOSH - try to generate the driving coordinates for the elimination reaction from stk or maybe just try to hardcode them for
    # the purpose of having something to demo
    func_group = list(stk_ligand_5a.get_functional_groups())[0]
    list(func_group.get_atom_ids())[0]
    for atom_info in 
    driving_coordinates = [['BREAK',1,54],['BREAK',1,58],['ADD',54,58]]
    stk_gsm(
        stk_mol=stk_mol,
        driving_coordinates=driving_coordinates,
        config=config,
    )

    # optimizer_sequence = stko.OptimizerSequence(stk.MCHammer(), stko.MetalOptimizer(), ASE(XTB()))
    # stk_mol = optimizer_sequence.optimize(stk_mol)
    
    # # update position matrix of ancillary ligand based on optimization while bound to dimethyl Pd
    # ancillary_ligand_infos = [atom_info for atom_info in stk_mol.get_atom_infos()
    #                           if atom_info.get_building_block_id() == 1]
    # ancillary_ligand_ids = [atom_info.get_atom().get_id() for atom_info in ancillary_ligand_infos]
    # ancillary_positions = np.row_stack(list(stk_mol.get_atomic_positions(ancillary_ligand_ids)))
    # stk_ancillary_ligand.with_position_matrix(ancillary_positions)
    
    # complex = bind_ligands(stk_metal('Pd'), stk_ancillary_ligand, stk_ligand_5a, stk_ligand_6a)
    # stk_list_to_xyz_file([complex], 'test_Pd_complex.xyz')
    assert True