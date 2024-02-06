from pathlib import Path
import numpy as np

import pytest
import stk
import stko
from xtb.ase.calculator import XTB
import os

from conformational_sampling.ase_stko_optimizer import ASE
from conformational_sampling.config import Config
from conformational_sampling.gsm import stk_gsm
from conformational_sampling.main import (ConformerEnsembleOptimizer, bind_ligands, bind_to_dimethyl_Pd, load_stk_mol, stk_list_to_xyz_file)
from conformational_sampling.utils import stk_metal


def suzuki():
    stk_ancillary_ligand = load_stk_mol(Path('example2_L1.xyz'))
    stk_ligand_5a = load_stk_mol(Path('example2_5a.xyz'))
    stk_ligand_6a = load_stk_mol(Path('example2_6a.xyz'))

    # specifies atoms of the ancillary ligand that bind to the metal, in this case as a smarts string
    functional_group_factory = stk.SmartsFunctionalGroupFactory(
        smarts='P',
        bonders=(0,),
        deleters=(),
    )
    stk_ancillary_ligand = stk.BuildingBlock.init_from_molecule(stk_ancillary_ligand, functional_groups=[functional_group_factory])
    # stk_mol = bind_to_dimethyl_Pd(stk_ancillary_ligand)

    # prepare functional groups for the reactive ligands
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
    stk_mol = bind_ligands(stk_metal('Pd'), stk_ancillary_ligand, stk_ligand_5a, stk_ligand_6a)
    stk_list_to_xyz_file([stk_mol], 'test_Pd_complex.xyz')
    
    config = Config(
        xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb',
        #ase_calculator=XTB(),
        max_dft_opt_steps=30,
        # num_cpus=16,
        dft_cpus_per_opt=16,
    )

    from ase.calculators.qchem import QChem
    os.environ['QCSCRATCH'] = os.environ['SLURM_LOCAL_SCRATCH']
    config.ase_calculator = QChem(
        method='PBE',
        # basis='6-31G',
        #basis='STO-3G',
        basis='LANL2DZ',
        ecp='fit-LANL2DZ',
        SCF_CONVERGENCE='5',
        SCF_MAX_CYCLES='300',
        SCF_ALGORITHM='DIIS',
        nt=config.dft_cpus_per_opt,
    )

    conformer_ensemble_optimizer = ConformerEnsembleOptimizer([stk_mol], config)
    optimized_stk_mol = conformer_ensemble_optimizer.optimize()[0]

    driving_coordinates = [('ADD',56,80),('BREAK',1,56),('BREAK',1,80)]
    stk_gsm(
        stk_mol=optimized_stk_mol,
        driving_coordinates=driving_coordinates,
        config=config,
    )
    assert True

if __name__ == "__main__":
    suzuki()