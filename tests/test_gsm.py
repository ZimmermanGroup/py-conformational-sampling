from dataclasses import dataclass
from pathlib import Path
import subprocess as sp
import pytest
import ase
from ase.calculators.calculator import Calculator
from ase.optimize import BFGS
from xtb.ase.calculator import XTB

import stk
import stko
from conformational_sampling.config import Config
from conformational_sampling.gsm import stk_gsm

from conformational_sampling.main import bind_to_dimethyl_Pd, load_stk_mol, stk_list_to_xyz_file, xtb_optimize

@dataclass
class ASE(stko.optimizers.Optimizer):
    calculator: Calculator
        
    def optimize(self, stk_mol):
        ase_mol = ase.Atoms(
            positions=list(stk_mol.get_atomic_positions()),
            numbers=[atom.get_atomic_number() for atom in stk_mol.get_atoms()]
        )
        ase_mol.calc = self.calculator
        opt = BFGS(ase_mol)
        opt.run(fmax=0.2)
        stk_mol.with_position_matrix(opt.atoms.get_positions())
        return stk_mol
            

@pytest.mark.xfail(raises=RuntimeError)
def test_gsm():
    # create dppe bound complex as test stk molecule
    ligand_path = Path('examples/dppe/ligand.xyz') # name of file containing ligand geometry
    stk_ligand = load_stk_mol(ligand_path)

    # specifies atoms of the ligand that bind to the metal, in this case as a smarts string
    functional_group_factory = stk.SmartsFunctionalGroupFactory(
        smarts='P',
        bonders=(0,),
        deleters=(),
    )
    stk_ligand = stk.BuildingBlock.init_from_molecule(stk_ligand, functional_groups=[functional_group_factory])
    stk_mol = bind_to_dimethyl_Pd(stk_ligand)
    stk_list_to_xyz_file([stk_mol], Path('test.xyz'))
    optimizer_sequence = stko.OptimizerSequence(stk.MCHammer(), stko.MetalOptimizer(), ASE(XTB()))
    stk_mol = optimizer_sequence.optimize(stk_mol)
    
    # result = sp.run(
    #             'which xtb',
    #             stdout=sp.PIPE,
    #             shell=True
    #         )
    # xtb_path = result.stdout.decode('UTF-8')[:-1]
    # print(xtb_path)
    # stk_mol = xtb_optimize(0, stk_mol, xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb')
    # stk_mol = xtb_optimize(0, stk_mol, xtb_path=xtb_path)
    
    stk_list_to_xyz_file([stk_mol], Path('test.xyz'))

    # run gsm to eliminate ethane
    # driving coordinates are 1-indexed
    driving_coordinates = [['BREAK',1,54],['BREAK',1,58],['ADD',54,58]]
    config = Config(
        # xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb',
        # ase_calculator=XTB(method='GFN-FF'),
        ase_calculator=XTB(),
    )
    # with pytest.raises(RuntimeError, ma) as excinfo:
    stk_gsm(
        stk_mol=stk_mol,
        driving_coordinates=driving_coordinates,
        config=config,
    )
    # assert str(excinfo) == "TS node shouldn't be the first or last node"        
