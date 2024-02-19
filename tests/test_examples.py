from pathlib import Path
import runpy
import shutil


def test_suzuki_example(tmp_path, monkeypatch):
    example_path = Path(__file__).parents[1] / 'examples' / 'suzuki'
    example_tmp_path = shutil.copytree(example_path, tmp_path, dirs_exist_ok=True)
    monkeypatch.chdir(example_tmp_path)
    example_vars = runpy.run_path(example_path / 'suzuki.py')
    assert True
    # example of how to read variables from resulting environment
    # assert example_vars['hello'] == 'hello world'
    