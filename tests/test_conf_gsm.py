import os
from pathlib import Path
from itertools import repeat
import logging
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
    max_dft_opt_steps=2,
    num_cpus=20,
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
    SCF_CONVERGENCE='5',
    nt=config.dft_cpus_per_opt,
    SCF_MAX_CYCLES='300',
    SCF_ALGORITHM='DIIS',
)

# generates conformers, performs multiple step optimization and uniqueness filtering
suzuki_ligand_conf_gen(stk_ligand_5a, stk_ligand_6a, stk_ancillary_ligand, config)

##############################################
#####  py-GSM run on all the conformers  #####
##############################################

def run_gsm(idx, stk_mol, driving_coordinates, config):
    #print(idx)
    path = Path(f'scratch/pystring_{idx}')
    path.mkdir(exist_ok=True)
    os.chdir(path)
    stk_gsm(
       stk_mol=stk_mol,
       driving_coordinates=driving_coordinates,
       config=config)
    os.chdir(Path('../..'))


# py-GSM configuration object
config = Config(
    xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb',
    ase_calculator=XTB(),
    #max_dft_opt_steps=2,
    num_cpus=20,
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

conformer_path = Path('suzuki_conformers.xyz')
conformer_mols = load_stk_mol_list(conformer_path)

driving_coordinates = [('ADD',56,80),('BREAK',1,56),('BREAK',1,80)]

# Running py-GSM in parallel with xTB
with ProcessPoolExecutor(max_workers=config.num_cpus) as executor:
    #create separatre directiories and run pygsm for individual conformers in its
    # corresponding directory
    executor.map(run_gsm, range(len(conformer_mols)), conformer_mols, repeat(driving_coordinates), repeat(config))

# Running py-GSM in parallel with DFT
# with ProcessPoolExecutor(max_workers=config.num_cpus//config.dft_cpus_per_opt) as executor:
#    executor.map(run_gsm, range(len(conformer_mols)), conformer_mols, repeat(driving_coordinates), repeat(config))