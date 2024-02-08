import os
from pathlib import Path
from itertools import repeat
import logging
import subprocess
FORMAT = "[%(asctime)s %(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s"
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

import stk
from xtb.ase.calculator import XTB
from concurrent.futures import ProcessPoolExecutor

from conformational_sampling.main import load_stk_mol, load_stk_mol_list, suzuki_ligand_conf_gen
from conformational_sampling.config import Config
from conformational_sampling.gsm import stk_gsm

##############################################
##  Conformer generation on Suzuki ligands  ##
##############################################

def conf_gsm():

    ligand_5a_path = Path('tests/test_data/example2_5a.xyz') # name of file containing ligand geometry
    ligand_6a_path = Path('tests/test_data/example2_6a.xyz')
    ancillary_ligand = Path('tests/test_data/example2_L1.xyz')
    stk_ligand_5a = load_stk_mol(ligand_5a_path)
    stk_ligand_6a = load_stk_mol(ligand_6a_path)
    stk_ancillary_ligand = load_stk_mol(ancillary_ligand)

    # specifies atoms of the ligand that bind to the metal, in this case as a smarts string
    functional_group_factory = stk.SmartsFunctionalGroupFactory(
        smarts='cBr',
        bonders=(0,),
        deleters=(1,),
    )
    stk_ligand_5a = stk.BuildingBlock.init_from_molecule(stk_ligand_5a, functional_groups=[functional_group_factory])
    functional_group_factory = stk.SmartsFunctionalGroupFactory(
        smarts='cB([O][H])[O][H]',
        bonders=(0,),
        deleters=tuple(range(1,6)),
    )
    stk_ligand_6a = stk.BuildingBlock.init_from_molecule(stk_ligand_6a, functional_groups=[functional_group_factory])
    functional_group_factory = stk.SmartsFunctionalGroupFactory(
        smarts='P',
        bonders=(0,),
        deleters=(),
    )
    stk_ancillary_ligand = stk.BuildingBlock.init_from_molecule(stk_ancillary_ligand, functional_groups=[functional_group_factory])

    # py-conformational-sampling configuration object
    config = Config(
        initial_conformers=50,
        xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb',
        #ase_calculator=XTB(),
        max_dft_opt_steps=10,
        num_cpus=28,
        dft_cpus_per_opt=4,
    )

    # qchem ase calculator setup
    from ase.calculators.qchem import QChem
    os.environ['QCSCRATCH'] = os.environ['SLURM_LOCAL_SCRATCH']
    config.ase_calculator = QChem(
        method='PBE',
        # basis='6-31G',
        # basis='STO-3G',
        basis='LANL2DZ',
        ecp='fit-LANL2DZ',
        SCF_CONVERGENCE='6',
        nt=config.dft_cpus_per_opt,
        SCF_MAX_CYCLES='500',
        SCF_ALGORITHM='RCA_DIIS',
        THRESH_RCA_SWITCH='4',
    )

    # generates conformers, performs multiple step optimization and uniqueness filtering
    # suzuki_ligand_conf_gen(stk_ligand_5a, stk_ligand_6a, stk_ancillary_ligand, config)

    ##############################################
    #####  py-GSM run on all the conformers  #####
    ##############################################

    # def run_gsm(idx, stk_mol, driving_coordinates, config):
    #     #print(idx)
    #     path = Path(f'scratch/pystring_{idx}')
    #     path.mkdir(exist_ok=True)
    #     os.chdir(path)
    #     stk_gsm(
    #        stk_mol=stk_mol,
    #        driving_coordinates=driving_coordinates,
    #        config=config)
    #     os.chdir(Path('../..'))

    conformer_path = Path('suzuki_conformers.xyz')
    conformer_mols = load_stk_mol_list(conformer_path)

    for i in range(len(conformer_mols)):
        path = Path.cwd() / f'scratch/pystring_{i}'
        path.mkdir(exist_ok=True)

    subprocess.run(['sbatch', f'--array=0-{len(conformer_mols)-1}', './tests/gsm_job_array.py'])
    # subprocess.run(['sbatch', f'--array=2', './tests/gsm_job_array.py'])

if __name__ == '__main__':
    conf_gsm()