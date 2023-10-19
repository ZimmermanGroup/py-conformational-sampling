from pytest import approx
import numpy as np
import stk

from conformational_sampling.analyze import tau_4_prime


def test_tau_4_prime():
    # tetrahedral
    methane = stk.BuildingBlock('C')
    assert tau_4_prime(methane.to_rdkit_mol(), 0) == approx(1, abs=0.02)
    
    # square planar    
    XeF4 = stk.BuildingBlock('F[Xe](F)(F)F')
    XeF4 = XeF4.with_position_matrix(np.array([
        [1, 0, 0],
        [0.0, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0]
    ]))
    assert tau_4_prime(XeF4.to_rdkit_mol(), 1) == approx (0, abs=0.02)
    
    # seesaw - RDKit's ETKDG geometry is too bad to verify whether this actually works
    # SF4 = stk.BuildingBlock('FS(F)(F)F')
    # assert tau_4_prime(SF4.to_rdkit_mol(), 1) == approx (0.24, abs=0.02)
    