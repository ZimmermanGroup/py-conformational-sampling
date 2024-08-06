import runpy
import shutil
from pathlib import Path

import pytest
import stk
from conformational_sampling.catalytic_reaction_complex import (
    ReductiveEliminationComplex,
)
from conformational_sampling.main import load_stk_mol
from conformational_sampling.utils import stk_metal


@pytest.mark.slow
def test_dppe_example(tmp_path, monkeypatch):
    example_path = Path(__file__).parents[1] / 'examples' / 'dppe'
    example_tmp_path = shutil.copytree(example_path, tmp_path, dirs_exist_ok=True)
    monkeypatch.chdir(example_tmp_path)
    example_vars = runpy.run_path(example_path / 'dppe.py')
    assert True


@pytest.mark.slow
def test_suzuki_example(tmp_path, monkeypatch):
    example_path = Path(__file__).parents[1] / 'examples' / 'suzuki'
    example_tmp_path = shutil.copytree(example_path, tmp_path, dirs_exist_ok=True)
    monkeypatch.chdir(example_tmp_path)
    example_vars = runpy.run_path(example_path / 'suzuki.py')
    assert True
    # example of how to read variables from resulting environment
    # assert example_vars['hello'] == 'hello world'


def test_reductive_elim_drive_coords(monkeypatch):
    example_path = Path(__file__).parents[1] / 'examples' / 'suzuki'
    monkeypatch.chdir(example_path)
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
    stk_ligand_5a = stk.BuildingBlock.init_from_molecule(
        stk_ligand_5a, functional_groups=[functional_group_factory]
    )
    functional_group_factory = stk.SmartsFunctionalGroupFactory(
        smarts='cB([O][H])[O][H]',
        bonders=(0,),
        deleters=tuple(range(1, 6)),
    )
    stk_ligand_6a = stk.BuildingBlock.init_from_molecule(
        stk_ligand_6a, functional_groups=[functional_group_factory]
    )
    functional_group_factory = stk.SmartsFunctionalGroupFactory(
        smarts='P',
        bonders=(0,),
        deleters=(),
    )
    stk_ancillary_ligand = stk.BuildingBlock.init_from_molecule(
        stk_ancillary_ligand, functional_groups=[functional_group_factory]
    )

    reaction_complex = ReductiveEliminationComplex(
        metal=stk_metal('Pd'),
        ancillary_ligand=stk_ancillary_ligand,
        reactive_ligand_1=stk_ligand_5a,
        reactive_ligand_2=stk_ligand_6a,
    )
    assert set(reaction_complex.gen_reductive_elim_drive_coords()) == {
        ('ADD', 56, 80),
        ('BREAK', 1, 56),
        ('BREAK', 1, 80),
    }
