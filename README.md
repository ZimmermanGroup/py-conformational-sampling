# py-conformational-sampling
py-conformational-sampling is an experimental python library for sampling conformers of a ligand and binding it to a metal.

![Conformational sampling](https://user-images.githubusercontent.com/5794446/184696183-d74131bb-960a-4632-805c-12e6ae92f536.png)

## Draft installation instructions
py-conformational-sampling is developed and primarily tested in a CentOS 7 linux high performance computing environment using anaconda and pip for python and other dependency installation. While not specifically tested, the code and dependencies aren't intentionally tied to this specific platform and are in principle portable to other platforms.

* Clone GitHub repository.

* Create a conda environment in which to install this library (links to more information on [anaconda installation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [anaconda environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)). The following command creates a new conda environment with python and the dependencies that the author was unable to install using pip (see below).
```
conda create -c conda-forge --name conformational-sampling python openbabel xtb-python
```

* Activate the conda environment.
```
conda activate conformational-sampling
```

* Navigate to the root directory of py-conformational-sampling. The directory should contain this project's setup.py which pip will look to for installation instructions. 

* Use pip to install py-conformational-sampling and its dependencies.
```
pip install -e .
```

* Optional: to install with extra packages for development and testing within visual studio code, use:
```
pip install -e .[dev,vscode]
```

## Getting started
* Navigate to the examples/dppe directory.
```
cd examples/dppe
```
<!-- 
* Ensure xTB is installed and accessible. In the author's environment this is done using linux modules.
```
module load xtb
``` -->
* This example contains slurm.py which can either be submitted as a slurm batch script or executed directly / interactively with python. To run interactively, acquire an interactive allocation of resources (e.g. interactive slurm job) and then execute:
```
python ./slurm.py
```

* Example output can be found in the examples/dppe_output directory.
