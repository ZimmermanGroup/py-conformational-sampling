# %%

# jupyter only imports
%reload_ext autoreload
%autoreload 2
from IPython.display import display

from pathlib import Path
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock, MolToXYZBlock
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

def bind_to_dimethyl_Pd(ligand, conf_id=-1):
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

    Chem.rdmolops.Kekulize(ligand)
    ligand_stk = stk.BuildingBlock.init_from_rdkit_mol(
        ligand,
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts='P',
                bonders=(0,),
                deleters=(),
            )
        ]
    )
    ligand_stk = ligand_stk.with_position_matrix(ligand.GetConformer(conf_id).GetPositions())

    return stk.ConstructedMolecule(
        topology_graph=MixedSquarePlanar(
            metals=metal,
            ligands={
                ligand_stk: (0, ),
                methyl: (1, 2),
            },
        ),
    )

def stk_list_to_xyz_file(stk_mol_list, file_path):
    with open(file_path, 'w') as file:
        for stk_mol in stk_mol_list:
            file.write(MolToXYZBlock(stk_mol.to_rdkit_mol()))

def gen_ligand_library_entry(mol):
    conf_ids = Chem.AllChem.EmbedMultipleConfs(mol, numConfs=100, randomSeed=40, pruneRmsThresh=0.6)
    complexes = [bind_to_dimethyl_Pd(mol, conf_id=conf_id) for conf_id in conf_ids]
    stk_list_to_xyz_file(complexes, 'test_output.xyz')
    # mc_hammer_complexes = [stk.MCHammer().optimize(complex) for complex in complexes]
    
    for complex in complexes:
        display_mol = complex.to_rdkit_mol()
        Chem.SanitizeMol(display_mol)
        display(nglview.show_rdkit(display_mol))
    
mol_path = Path('/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/examples/vitek_dmpp/ligand.xyz')
vitek_dmpp_ligand = load_mol(mol_path)
# vitek_dmpp_ligand = Path('/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/examples/bis_hydrazone/bis_hydrazone_ligand.xyz')
gen_ligand_library_entry(vitek_dmpp_ligand)
# for i in conf_ids:
#     display(nglview.show_rdkit(vitek_dmpp_ligand, conf_id=i))




# %%

complex = bind_to_dimethyl_Pd(vitek_dmpp_ligand, 3)
complex = stk.MCHammer().optimize(complex)
mol = complex.to_rdkit_mol()
Chem.SanitizeMol(mol)
display(nglview.show_rdkit(mol))
complex = stko.MetalOptimizer().optimize(complex)
mol = complex.to_rdkit_mol()
Chem.SanitizeMol(mol)
display(nglview.show_rdkit(mol))

complexes = [complex] * 2
from concurrent.futures import ProcessPoolExecutor
def execute_xtb(idx, complex):
    return stko.XTB(
        '/export/apps/CentOS7/xtb/xtb/bin/xtb',
        output_dir=f'xtb_test_{idx}',
        calculate_hessian=False,
        max_runs=1,
        charge=0
    ).optimize(complex)
with ProcessPoolExecutor(max_workers=2) as executor:
    optimized_complexes = list(executor.map(execute_xtb, range(len(complexes)), complexes))
mol = optimized_complexes[0].to_rdkit_mol()
Chem.SanitizeMol(mol)
display(nglview.show_rdkit(mol))

# %%
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
with SLURMCluster(
    queue='zimintel',
    cores=1,
    memory='1GB',
    processes=1
) as cluster, Client(cluster) as client:
    cluster.scale(jobs=2)
    # futures = client.map(execute_xtb, range(len(complexes)), complexes)
    # mol = futures[0].result(timeout=10).to_rdkit_mol()
# Chem.SanitizeMol(mol)
# display(nglview.show_rdkit(mol))
