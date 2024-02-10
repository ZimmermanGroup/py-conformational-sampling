import os
from pathlib import Path
import runpy


def test_examples():
    assert True
    example_vars = runpy.run_path(Path(__file__).parents[1] / 'examples' / 'gsm' / 'simple_gsm.py')
    assert example_vars['hello'] == 'hello world'


def test_suzuki_example():
    example_path = Path(__file__).parents[1] / 'examples' / 'suzuki'
    os.chdir(example_path)
    example_vars = runpy.run_path(example_path / 'suzuki.py')
    