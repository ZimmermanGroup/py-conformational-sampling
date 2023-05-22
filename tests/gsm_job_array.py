#!/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/.venv39/bin/python
#SBATCH -p zimgpu --job-name=py_gsm
#SBATCH -c1
#SBATCH --time=5-00:00:00
#SBATCH -o scratch/pystring_%a/output.txt
#SBATCH --array=0

import os
from pathlib import Path
from xtb.ase.calculator import XTB
from conformational_sampling.config import Config

from conformational_sampling.gsm import stk_gsm
from conformational_sampling.main import load_stk_mol_list

job_index = int(os.environ['SLURM_ARRAY_TASK_ID'])

conformer_path = Path('suzuki_conformers.xyz')
conformer_mols = load_stk_mol_list(conformer_path)

driving_coordinates = [('ADD',56,80),('BREAK',1,56),('BREAK',1,80)]


# create a directory for each specific job instance in the array
path = Path(f'scratch/pystring_{job_index}')
path.mkdir(exist_ok=True)
os.chdir(path)

# py-GSM configuration object
config = Config(
    xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb',
    ase_calculator=XTB(),
    #max_dft_opt_steps=2,
    num_cpus=1,
    #dft_cpus_per_opt=4,
)

# qchem ase calculator setup
# from ase.calculators.qchem import QChem
# os.environ['QCSCRATCH'] = os.environ['SLURM_LOCAL_SCRATCH']
# config.ase_calculator = QChem(
#     method='PBE',
#     # basis='6-31G',
#     # basis='STO-3G',
#     basis='LANL2DZ',
#     ecp='fit-LANL2DZ',
#     SCF_CONVERGENCE='5',
#     nt=config.dft_cpus_per_opt,
#     SCF_MAX_CYCLES='300',
#     SCF_ALGORITHM='DIIS',
# )

stk_gsm(
    stk_mol=conformer_mols[job_index],
    driving_coordinates=driving_coordinates,
    config=config,
)
