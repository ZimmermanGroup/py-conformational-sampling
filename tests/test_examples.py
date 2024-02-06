from pathlib import Path
import runpy


def test_examples():
    assert True
    example_vars = runpy.run_path(Path(__file__).parents[1] / 'examples' / 'gsm' / 'simple_gsm.py')
    assert example_vars['hello'] == 'hello world'