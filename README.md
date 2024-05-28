py-conformational-sampling
==========================

py-conformational-sampling is an experimental python library for sampling conformers of a metal ligand complex and searching for a reaction path and corresponding transition state for each conformer.

![Conformational sampling](https://github.com/ZimmermanGroup/py-conformational-sampling/assets/5794446/ac17c431-3a02-4fb0-a923-d2ea1459b2d8)

The user provides an ancillary ligand and reactive ligands as structure files which this library reads and stores as [stk (supramolecular toolkit) objects](https://stk.readthedocs.io/en/stable/stk.molecular.molecules.building_block.html). The user also provides a description of binding atoms and a description of a reaction in the form of changes in bonding. This library generates an ensemble of conformers for each ligand, binds them to Pd, and performs an optimization and filtering funnel to refine the conformer ensemble. Then for each conformer, the user's reaction is carried out using [pyGSM](https://github.com/ZimmermanGroup/pyGSM) which looks for a reaction path and corresponding transition state.


# Draft installation instructions

Note: the code for this project is written in python which is generally platform independent. The author is unaware of any dependencies that are tied to a specific operating system. However, currently development and testing are primarily conducted within a CentOS 7 linux high performance computing environment using pip and venv or anaconda for management of python and other dependencies. Example commands below are for linux and will vary slightly by platform. Python versions 3.8 and 3.11 are tested in github actions workflows.

## Venv (built in virtual environment) installation

Note: [on Windows, the dependency openbabel may be difficult to install](https://github.com/openbabel/openbabel/issues/2408#issuecomment-1288847122).

* Use [git](https://git-scm.com/) to clone GitHub repository into a directory location of your choice.

```
git clone https://github.com/ZimmermanGroup/py-conformational-sampling.git
```

* From the command line, navigate to the directory where this git repository was cloned (py-conformational-sampling). The remaining sample commands are based on being run from this directory.

Create a [python](https://www.python.org/downloads/) virtual environment using the built in python module [venv](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment) as an isolated environment in which to install this library and its dependencies.

```
python -m venv .venv
```

* Activate the newly created environment. The environment must be activated whenever this library will be executed. Activation can be manual or automated through a job submission script or shell configuration file (e.g. .bashrc) as desired.

```
source ./.venv/bin/activate
```

* Use pip to install py-conformational-sampling and its dependencies.

```
pip install -e .
```

## Conda (alternative) installation

* Clone GitHub repository.
* Create a conda environment in which to install this library (links to more information on [anaconda installation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [anaconda environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)). The following command creates a new conda environment with python and a couple dependencies that the author was unable to install automatically using pip (see below).

```
conda create -c conda-forge --name conformational-sampling python openbabel xtb-python
```

Note: in the author's environment, manually compiling a version of xTB was found to run a few times faster than the precompiled binary from conda used in these instructions.

* Activate the conda environment.

```
conda activate conformational-sampling
```

* Navigate to the root directory of this project (py-conformational-sampling). The directory should contain this project's pyproject.toml which pip will look to for installation instructions.
* Use pip to install py-conformational-sampling and its dependencies.

```
pip install -e .
```

## For developers:

To install in editable mode with extra packages for development and testing within visual studio code, substitute the following during the pip installation:

```
pip install -e .[dev,vscode]
```

## Updating an existing installation:
If you have made changes to the example scripts or other files in this repository that you may want to keep, before updating you should first copy those files to a folder location outside this directory. If you are actively using this library but wish to try the latest version, having multiple versions of this repository in different difectory locations is an option. To update run:

```
# REMOVES ANY CHANGES YOU HAVE MADE TO FILES IN THIS REPOSITORY
git reset --hard

# updates this library to the latest version
git pull
pip install .
```


# Getting started

## Suzuki example from [Organic Letters publication](https://pubs.acs.org/doi/10.1021/acs.orglett.3c04047)

* Navigate to the examples/suzuki directory

```
cd examples/suzuki
```

* Ensure that the venv or conda environment containing py-conformational-sampling is activated.
* Execute the example:

```
python ./suzuki.py
```

## Interactive visualization example

Note: assumes that the visualization is running on a remote computing cluster. If serving the visualization locally, omit the ssh tunneling step

On the cluster containing the results of py-conformational-sampling calculations:
* Ensure that the venv or conda environment containing py-conformational-sampling is activated.
* Navigate to examples/dppe within py-conformational-sampling
* Set `start_visualization = True` in dppe.py, then run the script to start the server which loads the molecular data and hosts the visualization:

```
python dppe.py
```

* Establish an ssh SSH tunnel to the cluster in a separate terminal tab or window on the local machine with an available web browser to view the visualization:
```
ssh -NfL localhost:5006:localhost:5006 <user@remote.host>
```
* Access the visualization through a web browser at the following address:
```
http://localhost:5006
```

Note: the visualization is built using the Panel package. These instructions are based on the Panel documentation at https://panel.holoviz.org/how_to/server/index.html

## Examples of only conformer generation step (outdated)

### DPPE

* Navigate to the examples/dppe directory.

```
cd examples/dppe
```

* This example contains slurm.py which provides an example workflow that can serve as a template interface for working with the library. This python script can either be submitted as a slurm batch script or executed directly / interactively with python. To run interactively, acquire an interactive allocation of resources (e.g. interactive slurm job) and then execute:

```
python ./slurm.py
```

* Example output can be found in the examples/dppe_output directory. XYZ output files contain the conformer ensemble after each stage of optimization.

### Alonso's ligand

* Example is in the alonso_ligand directory. Execution is similar to DPPE example.

```
cd examples/alonso_ligand
python ./slurm.py
```

Output notes:

* For brevity of testing the workflow on the examples, `config.max_dft_opt_steps = 2` to test the integration and configuration of the DFT calculator in minutes. For fully optimized structures, increase to a larger number (e.g. 100) to allow most geometry optimizations to fully converge. This will likely take hours.
* Duplicate conformers are filtered out before and during the optimization pipeline to reduce the computational burden of higher level optimizations. Due to this, the number of output conformers may be only a small fraction of `config.initial_conformers` especially for smaller or more rigid molecules. Threshold for uniqueness can be configured with `config.pre_xtb_rms_threshold`.
* Conformers are ordered in the output to display relevant conformers first. Conformers with no more than 2 changes in bonding (3 conformers in the sample output) are output first. This is customizable with `config.max_connectivity_changes`. Conformers are then sorted by energy, with the lowest energy first.
* Current openbabel implementation may not maintain all types of stereochemistry such as atropisomerism.

## Using your own ligands

* Create a directory and copy files containing your ligand structures into the directory.
* Copy one of the example python files from the examples folder into the directory.
* If necessary or desired, customize the config and description of ligand binding

Notes:

* The venv or conda environment containing py-conformational-sampling should be activated when using the library
* Since openbabel is used to interpret the input file, any openbabel supported molecular file format can be used as input with slight modification (e.g. `stk_ligand = load_stk_mol(ligand_path, fmt='mol')` for a mol file)
* It is recommended that the environment variable OMP_NUM_THREADS is set to "1". If unset, it seems that xTB attempts to use multiple threads for each individual conformer, and this can be an order of magnitude less efficient than parallelizing over conformers. Py-conformational-sampling sets OMP_NUM_THREADS in its package's `__init__.py` file, but this will not take effect if you import xtb before conformational_sampling
* Binding atoms may be alternatively supplied by specifying the element and index of binding atoms based on zero-indexed ordering in the ligand structure file (example excerpt below)

```python
functional_groups = [stk.SingleAtom(stk.C(73)), stk.SingleAtom(stk.O(22))]
stk_ligand = stk.BuildingBlock.init_from_molecule(
    stk_ligand,
    functional_groups=functional_groups
)
```

