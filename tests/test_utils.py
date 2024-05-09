from pytest import approx
from conformational_sampling.utils import free_energy_diff

def test_free_energy_diff():
    assert free_energy_diff([5.0],[7.5],temperature=298) == approx(2.5)
    assert free_energy_diff([5,6], [5,6], temperature=100) == approx(0)
    # approx value from wolfram alpha
    assert free_energy_diff([5,10], [6,7], temperature=298) == approx(0.89972418865667)
