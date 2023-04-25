from pathlib import Path

import pytest
import stk
import stko
from xtb.ase.calculator import XTB

from conformational_sampling.ase_stko_optimizer import ASE
from conformational_sampling.config import Config
from conformational_sampling.gsm import stk_gsm
from conformational_sampling.main import (bind_to_dimethyl_Pd, load_stk_mol)


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
