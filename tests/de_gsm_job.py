#!/export/zimmerman/soumikd/py-conformational-sampling/.venv/bin/python
#SBATCH -p zimintel --job-name=pygsm_l1_symm_degsm_dft
#SBATCH -c8
#SBATCH --time=10-00:00:00
#SBATCH -o de_dft_output.txt
#SBATCH -e de_dft_error.txt
#SBATCH --array=0

import os
from pathlib import Path
from xtb.ase.calculator import XTB
from conformational_sampling.config import Config

from conformational_sampling.gsm import stk_gsm, stk_gsm_command_line, stk_se_gsm, stk_de_gsm
from conformational_sampling.main import load_stk_mol_list

# job_index = int(os.environ['SLURM_ARRAY_TASK_ID'])

# conformer_path = Path('suzuki_conformers.xyz')
# conformer_mols = load_stk_mol_list(conformer_path)

# driving_coordinates = [('ADD',74,98),('BREAK',1,74),('BREAK',1,98)]


# create a directory for each specific job instance in the array
# path = Path.cwd() / f'scratch/pystring_{job_index}'
# path.mkdir(exist_ok=True)
# os.chdir(path)

# py-GSM configuration object for xTB
# config1 = Config(
#     xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb',
#     ase_calculator=XTB(),
#     # max_dft_opt_steps=2,
#     num_cpus=1,
#     dft_cpus_per_opt=1,
# )

# py-GSM configuration object for DFT
config2 = Config(
    # xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb',
    # ase_calculator=XTB(),
    max_dft_opt_steps=2,
    num_cpus=8,
    dft_cpus_per_opt=8,
)

# qchem ase calculator setup
from ase.calculators.qchem import QChem
os.environ['QCSCRATCH'] = os.environ['SLURM_LOCAL_SCRATCH']
config2.ase_calculator = QChem(
    method='B97-D',
    # DFT_D='EMPIRICAL_GRIMME3',
    # basis='6-31G',
    # basis='STO-3G',
    basis='GEN',
    ecp='fit-LANL2DZ',
    SCF_CONVERGENCE='6',
    nt=config2.dft_cpus_per_opt,
    SCF_MAX_CYCLES='500',
    SCF_ALGORITHM='RCA_DIIS',
    THRESH_RCA_SWITCH='4',
    basisfile='../../../basis',
    # ecpfile='../../../ecp',
)

# look for a partially completed string from which to restart pygsm
# if not (path / 'opt_converged_000.xyz').exists():
#     try:
#         config1.restart_gsm = sorted(path.glob('scratch/opt_iters_000_*.xyz'))[-1]
#     except IndexError: # no pygsm opt iters files exist, which should always occur in an initial run
#         pass

    # stk_gsm(
    #     stk_mol=conformer_mols[job_index],
    #     driving_coordinates=driving_coordinates,
    #     config=config,
    # )

    # stk_gsm_command_line(
    #     stk_mol=conformer_mols[job_index],
    #     driving_coordinates=driving_coordinates,
    #     config=config,
    # )

# stk_se_gsm(
#     stk_mol=conformer_mols[job_index],
#     driving_coordinates=driving_coordinates,
#     config=config1,
# )
    
stk_de_gsm(config=config2)