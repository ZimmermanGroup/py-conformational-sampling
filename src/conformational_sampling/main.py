import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from openbabel import pybel as pb
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock, MolToXYZBlock

import stk
import stko

from conformational_sampling.metal_complexes import OneLargeTwoSmallMonodentateTrigonalPlanar, TwoMonoOneBidentateSquarePlanar

XTB_PATH = '/export/apps/CentOS7/xtb/xtb/bin/xtb'

def num_cpus():
    try:
        return int(os.environ["SLURM_CPUS_PER_TASK"])
    except KeyError:
        return 2

def load_stk_mol(molecule_path, fmt='xyz'):
    'loads a molecule from a file via pybel into an stk Molecule object'
    pybel_mol = next(pb.readfile(fmt, str(molecule_path)))
    rdkit_mol = MolFromMolBlock(pybel_mol.write('mol'), removeHs=False)
    Chem.rdmolops.Kekulize(rdkit_mol)
    stk_mol = stk.BuildingBlock.init_from_rdkit_mol(rdkit_mol)
    return stk_mol

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
    
    if ligand.get_num_functional_groups() == 1: #monodentate
        return stk.ConstructedMolecule(
            topology_graph=OneLargeTwoSmallMonodentateTrigonalPlanar(
                metals=metal,
                ligands={
                    ligand: (0, ),
                    methyl: (1, 2),
                },
            ),
        )
    elif ligand.get_num_functional_groups() == 2: #bidentate
        return stk.ConstructedMolecule(
            topology_graph=TwoMonoOneBidentateSquarePlanar(
                metals=metal,
                ligands={
                    ligand: (0, ),
                    methyl: (1, 2),
                },
            ),
        )


def stk_list_to_xyz_file(stk_mol_list, file_path, energies=None):
    with open(file_path, 'w') as file:
        for i, stk_mol in enumerate(stk_mol_list):
            rdkit_mol = stk_mol.to_rdkit_mol()
            if energies is not None:
                rdkit_mol.SetProp('_Name', str(energies[i]))
            file.write(MolToXYZBlock(rdkit_mol))

def xtb_optimize(idx, complex):
    return stko.XTB(
        XTB_PATH,
        output_dir=f'scratch/xtb_optimize_{idx}',
        calculate_hessian=False,
        max_runs=1,
        charge=0
    ).optimize(complex)

def xtb_energy(idx, complex):
    return stko.XTBEnergy(
        XTB_PATH,
        output_dir=f'scratch/xtb_energy_{idx}',
    ).get_energy(complex)

def gen_ligand_library_entry(stk_ligand, numConfs=100):
    rdkit_mol = stk_ligand.to_rdkit_mol()
    conf_ids = Chem.AllChem.EmbedMultipleConfs(rdkit_mol, numConfs=numConfs, randomSeed=40, pruneRmsThresh=0.6, numThreads=num_cpus())
    stk_conformers = [stk_ligand.with_position_matrix(rdkit_mol.GetConformer(conf_id).GetPositions())
                   for conf_id in conf_ids]
    stk_list_to_xyz_file(stk_conformers, 'conformers_ligand_only.xyz')
    unoptimized_complexes = [bind_to_dimethyl_Pd(ligand) for ligand in stk_conformers]
    stk_list_to_xyz_file(unoptimized_complexes, 'conformers_0_unoptimized.xyz')
    mc_hammer_complexes = [stk.MCHammer().optimize(complex) for complex in unoptimized_complexes]
    stk_list_to_xyz_file(mc_hammer_complexes, 'conformers_1_mc_hammer.xyz')
    metal_optimizer_complexes = [stko.MetalOptimizer().optimize(complex)
                                 for complex in mc_hammer_complexes]
    stk_list_to_xyz_file(metal_optimizer_complexes, 'conformers_2_metal_optimizer.xyz')
    (Path.cwd() / 'scratch').mkdir(exist_ok=True)
    with ProcessPoolExecutor(max_workers=num_cpus()) as executor:
        xtb_complexes = list(executor.map(xtb_optimize, range(len(metal_optimizer_complexes)),
                                          metal_optimizer_complexes))
        xtb_energies = list(executor.map(xtb_energy, range(len(metal_optimizer_complexes)),
                                          metal_optimizer_complexes))
    stk_list_to_xyz_file(xtb_complexes, 'conformers_3_xtb.xyz', energies=xtb_energies)
    