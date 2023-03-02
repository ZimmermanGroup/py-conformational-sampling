import pathlib
import setuptools

README = (pathlib.Path(__file__).parent / "README.md").read_text()

setuptools.setup(
    name="py-conformational-sampling",
    version="0.2.0",
    # description="",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/ZimmermanGroup/py-conformational-sampling",
    author="Joshua Kammeraad",
    # author_email="",
    license="MIT",
    package_dir={"": "src"},
    packages = setuptools.find_packages(where="src"),
    python_requires=">=3.8",
    install_requires=(
        'pytest',
        'numpy',
        'pandas',
        'param',
        'hvplot',
        'panel>=0.13.1',
        'panel-chemistry',
        'rdkit',
        'openbabel-wheel',
        'stk',
        'stko',
        'ase',
        'pyGSM @ git+https://github.com/ZimmermanGroup/pyGSM.git',
    ),
    extras_require={
        # "dev": ["sphinx", "sphinx_rtd_theme", "pytest", "coverage", "pytest-mock"],
        "dev": ["snakeviz", "dask-jobqueue"],
        "vscode": ["pylint", "jupyter"]
    },
)