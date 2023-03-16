import numpy as np
import stk
from conformational_sampling.main import pybel_mol_to_stk_mol, stk_mol_to_pybel_mol

def test_imports():
    assert True

def test_stk_to_pybel_loop():
    'convert an stk molecule to a pybel molecule and then back to validate the conversion methods'
    butane_initial = stk.BuildingBlock('CCCC')
    butane_final = pybel_mol_to_stk_mol(stk_mol_to_pybel_mol(butane_initial))
    initial_positions = butane_initial.get_position_matrix()
    final_positions = butane_final.get_position_matrix()
    assert np.allclose(initial_positions, final_positions, atol=1e-3)
                                        