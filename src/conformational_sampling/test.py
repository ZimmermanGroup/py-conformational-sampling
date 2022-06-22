# %%

# jupyter only imports
%reload_ext autoreload
%autoreload 2

import os
# from IPython.display import display

from pathlib import Path
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock, MolToXYZBlock
from concurrent.futures import ProcessPoolExecutor

import openbabel as ob
from openbabel import pybel as pb
from rdkit import Chem
import nglview

import stk
import stko
from conformational_sampling.mixed_square_planar import MixedSquarePlanar

def load_mol(molecule_path, fmt='xyz'):
    'loads a molecule from a file via pybel into an rdkit Mol object'
    pybel_mol = next(pb.readfile(fmt, str(molecule_path)))
    rdkit_mol = MolFromMolBlock(pybel_mol.write('mol'), removeHs=False)
    return rdkit_mol

def rdkit_conf_to_building_block(rdkit_mol, conf_id=-1):
    Chem.rdmolops.Kekulize(rdkit_mol)
    stk_mol = stk.BuildingBlock.init_from_rdkit_mol(
        rdkit_mol,
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts='P',
                bonders=(0,),
                deleters=(),
            )
        ]
    )
    return stk_mol.with_position_matrix(rdkit_mol.GetConformer(conf_id).GetPositions())

def bind_to_dimethyl_Pd(ligand):
    metal = stk.BuildingBlock(
        smiles='[Pd]',
        functional_groups=(
            stk.SingleAtom(stk.Pd(0))
            for i in range(6)
        ),
        position_matrix=[[0, 0, 0]],
    )

    methyl = stk.BuildingBlock(
        smiles='[CH3]',
        functional_groups=[
            stk.SingleAtom(stk.C(0))
        ]
    )

    return stk.ConstructedMolecule(
        topology_graph=MixedSquarePlanar(
            metals=metal,
            ligands={
                ligand: (0, ),
                methyl: (1, 2),
            },
        ),
    )

def stk_list_to_xyz_file(stk_mol_list, file_path):
    with open(file_path, 'w') as file:
        for stk_mol in stk_mol_list:
            file.write(MolToXYZBlock(stk_mol.to_rdkit_mol()))

def execute_xtb(idx, complex):
    return stko.XTB(
        '/export/apps/CentOS7/xtb/xtb/bin/xtb',
        output_dir=f'xtb_test_{idx}',
        calculate_hessian=False,
        max_runs=1,
        charge=0
    ).optimize(complex)

def gen_ligand_library_entry(mol):
    conf_ids = Chem.AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=40, pruneRmsThresh=0.6, numThreads=2)
    stk_ligands = [rdkit_conf_to_building_block(mol, conf_id) for conf_id in conf_ids]
    unoptimized_complexes = [bind_to_dimethyl_Pd(ligand) for ligand in stk_ligands]
    stk_list_to_xyz_file(unoptimized_complexes, 'conformers_0_unoptimized.xyz')
    mc_hammer_complexes = [stk.MCHammer().optimize(complex) for complex in unoptimized_complexes]
    stk_list_to_xyz_file(unoptimized_complexes, 'conformers_1_mc_hammer.xyz')
    metal_optimizer_complexes = [stko.MetalOptimizer().optimize(complex)
                                 for complex in mc_hammer_complexes]
    stk_list_to_xyz_file(mc_hammer_complexes, 'conformers_2_metal_optimizer.xyz')
    # with ProcessPoolExecutor(max_workers=2) as executor:
    #     xtb_complexes = list(executor.map(execute_xtb, range(len(metal_optimizer_complexes)),
    #                                       metal_optimizer_complexes))
    # stk_list_to_xyz_file(xtb_complexes, 'conformers_3_xtb.xyz')
    
# os.chdir(Path('/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/examples/vitek_dmpp/'))
os.chdir(Path('/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/examples/alonso_ligand/'))
mol_path = Path('ligand.xyz')
vitek_dmpp_ligand = load_mol(mol_path)
# vitek_dmpp_ligand = Path('/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/examples/bis_hydrazone/bis_hydrazone_ligand.xyz')
gen_ligand_library_entry(vitek_dmpp_ligand)
