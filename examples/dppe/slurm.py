#!/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/.venv/bin/python
#SBATCH -p zimintel --job-name=conformational_sampling
#SBATCH -c28
#SBATCH --time=1-0
#SBATCH -o output.txt

import os
from pathlib import Path
import logging
logging.basicConfig(level=logging.DEBUG)

import stk
from ase.calculators.qchem import QChem

from conformational_sampling.main import load_stk_mol, gen_ligand_library_entry
from conformational_sampling.config import Config

ligand_path = Path('ligand.xyz') # name of file containing ligand geometry
stk_ligand = load_stk_mol(ligand_path)

# specifies atoms of the ligand that bind to the metal, in this case as a smarts string
functional_group_factory = stk.SmartsFunctionalGroupFactory(
    smarts='P',
    bonders=(0,),
    deleters=(),
)
stk_ligand = stk.BuildingBlock.init_from_molecule(stk_ligand, functional_groups=[functional_group_factory])

# qchem ase calculator setup
os.environ['QCSCRATCH'] = os.environ['SLURM_LOCAL_SCRATCH']
ase_calculator = QChem(
    method='PBE',
    # basis='6-31G',
    # basis='STO-3G',
    basis='LANL2DZ',
    ecp='fit-LANL2DZ',
    SCF_CONVERGENCE='5',
    nt=4
)

# py-conformational-sampling configuration object
config = Config(
    initial_conformers=10,
    xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb',
    ase_calculator=ase_calculator,
)

# generates conformers, performs multiple step optimization and uniqueness filtering
gen_ligand_library_entry(stk_ligand, config)
