from pathlib import Path
import runpy


def test_examples():
    assert True
    example_vars = runpy.run_path(Path(__file__).parents[1] / 'examples' / 'gsm' / 'simple_gsm.py')
    assert example_vars['hello'] == 'hello world'


def test_suzuki_example():
    # CHANGE DIRECTORY INTO THE SUZUKI EXAMPLE BEFORE EXECUTING THE SCRIPT
    example_vars = runpy.run_path(Path(__file__).parents[1] / 'examples' / 'suzuki' / 'suzuki.py')
    