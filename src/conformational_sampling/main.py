from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

from openbabel import pybel as pb
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock, MolToXYZBlock

import stk
import stko

from conformational_sampling.utils import num_cpus
from conformational_sampling.metal_complexes import (OneLargeTwoSmallMonodentateTrigonalPlanar,
                                                     TwoMonoOneBidentateSquarePlanar)

UNOPTIMIZED = 0
MC_HAMMER = 1
METAL_OPTIMIZER = 2
XTB = 3

NAMES = {UNOPTIMIZED: 'unoptimized', MC_HAMMER: 'mc_hammer', METAL_OPTIMIZER: 'metal_optimizer', XTB: 'xtb'}

class ConformerOptimizationSequence:
    def __init__(self, unoptimized) -> None:
        self.stages = {UNOPTIMIZED: unoptimized}
        self.energies = {}
        self.num_connectivity_changes = None

class ConformerEnsembleOptimizer:
    def __init__(self, unoptimized_conformers, config) -> None:
        self.conformers = [ConformerOptimizationSequence(conformer) for conformer in unoptimized_conformers]
        self.config = config
    
    def order_conformers(self):
        metal_optimized_conformers = []
        xtb_conformers = []
        final_conformers = []
        for conformer in self.conformers:
            if (conformer.num_connectivity_changes is not None
                and conformer.num_connectivity_changes <= self.config.max_connectivity_changes):
                final_conformers.append(conformer)
            elif XTB in conformer.stages:
                xtb_conformers.append(conformer)
            else: # only made it to metal optimizer stage
                metal_optimized_conformers.append(conformer)
        self.conformers = sorted(final_conformers, key=lambda conformer: conformer.energies[XTB])
        self.conformers += sorted(xtb_conformers, key=lambda conformer: conformer.energies[XTB])
        self.conformers += metal_optimized_conformers
    
    def get_unique_conformer_ids(self, stage):
        rdkit_mols = {i: Chem.RemoveHs(conformer.stages[stage].to_rdkit_mol())
                      for i, conformer in enumerate(self.conformers)
                      if stage in conformer.stages}
        unique_indices = []
        for i, rdkit_mol in rdkit_mols.items():
            for j in unique_indices:
                if Chem.AllChem.CalcRMS(rdkit_mol, rdkit_mols[j]) < self.config.pre_xtb_rms_threshold:
                    break
            else: # should only get here when conformer is sufficiently unique
                unique_indices.append(i)
        
        return unique_indices

    def optimize(self):
        with ProcessPoolExecutor(max_workers=num_cpus()) as executor:
            unoptimized_complexes = [conformer.stages[UNOPTIMIZED] for conformer in self.conformers]
            mc_hammer_complexes = list(executor.map(stk.MCHammer().optimize, unoptimized_complexes))
            for i, conformer in enumerate(self.conformers):
                conformer.stages[MC_HAMMER] = mc_hammer_complexes[i]
            metal_optimizer_complexes = list(executor.map(stko.MetalOptimizer().optimize, mc_hammer_complexes))
            for i, conformer in enumerate(self.conformers):
                conformer.stages[METAL_OPTIMIZER] = metal_optimizer_complexes[i]
            self.write()
        
            # remove duplicate molecules before running xTB
            unique_ids = self.get_unique_conformer_ids(METAL_OPTIMIZER)
            unique_conformers = [self.conformers[i] for i in unique_ids]
            metal_optimizer_complexes = [conformer.stages[METAL_OPTIMIZER] for conformer in unique_conformers]

            # run xTB on conformers in parallel
            (Path.cwd() / 'scratch').mkdir(exist_ok=True)
            xtb_complexes = list(executor.map(xtb_optimize, range(len(metal_optimizer_complexes)),
                                            metal_optimizer_complexes, repeat(self.config.xtb_path)))
            # compute energies
            energies = list(executor.map(xtb_energy, range(len(xtb_complexes)),
                                         xtb_complexes, repeat(self.config.xtb_path)))
            for i, conformer in enumerate(unique_conformers):
                conformer.stages[XTB] = xtb_complexes[i]
                conformer.energies[XTB] = energies[i]
            
            # compute number of connectivity changes to put most relevant conformers first in output
            xtb_complexes = list(executor.map(reperceive_bonds, xtb_complexes))
            for i, conformer in enumerate(unique_conformers):
                conformer.num_connectivity_changes = num_connectivity_differences(
                    unoptimized_complexes[0],
                    xtb_complexes[i]
                )
            
            # order conformers with the most relevant first and write to output file
            self.order_conformers()
            self.write()
    
    def write(self):
        for i, name in NAMES.items():
            complexes = [conformer.stages[i] for conformer in self.conformers if i in conformer.stages]
            energies = [conformer.energies[i] if i in conformer.energies else None
                        for conformer in self.conformers]
            with open(f'conformers_{i}_{name}.xyz', 'w') as file:
                for i, (stk_mol, energy) in enumerate(zip(complexes, energies)):
                    rdkit_mol = stk_mol.to_rdkit_mol()
                    if energy is not None:
                        rdkit_mol.SetProp('_Name', str(energy))
                    file.write(MolToXYZBlock(rdkit_mol))


def pybel_mol_to_stk_mol(pybel_mol):
    rdkit_mol = MolFromMolBlock(pybel_mol.write('mol'), removeHs=False)
    Chem.rdmolops.Kekulize(rdkit_mol)
    stk_mol = stk.BuildingBlock.init_from_rdkit_mol(rdkit_mol)
    return stk_mol

def stk_mol_to_pybel_mol(stk_mol, reperceive_bonds=False):
    if reperceive_bonds:
        return pb.readstring('xyz', MolToXYZBlock(stk_mol.to_rdkit_mol()))
    else:
        return pb.readstring('mol', MolToMolBlock(stk_mol.to_rdkit_mol()))
    

def load_stk_mol(molecule_path, fmt='xyz'):
    'loads a molecule from a file via pybel into an stk Molecule object'
    pybel_mol = next(pb.readfile(fmt, str(molecule_path)))
    return pybel_mol_to_stk_mol(pybel_mol)

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
        MetalComplexClass = OneLargeTwoSmallMonodentateTrigonalPlanar
    elif ligand.get_num_functional_groups() == 2: #bidentate
        MetalComplexClass = TwoMonoOneBidentateSquarePlanar
        
    return stk.ConstructedMolecule(
        topology_graph=MetalComplexClass(
            metals=metal,
            ligands={
                ligand: (0, ),
                methyl: (1, 2),
            },
        ),
    )


def stk_list_to_xyz_file(stk_mol_list, file_path):
    with open(file_path, 'w') as file:
        for i, stk_mol in enumerate(stk_mol_list):
            rdkit_mol = stk_mol.to_rdkit_mol()
            file.write(MolToXYZBlock(rdkit_mol))

def xtb_optimize(idx, complex, xtb_path):
    return stko.XTB(
        xtb_path,
        output_dir=f'scratch/xtb_optimize_{idx}',
        calculate_hessian=False,
        max_runs=1,
        charge=0
    ).optimize(complex)

def xtb_energy(idx, complex, xtb_path):
    return stko.XTBEnergy(
        xtb_path,
        output_dir=f'scratch/xtb_energy_{idx}',
    ).get_energy(complex)
    
def reperceive_bonds(stk_mol):
    # output to xyz and read in with pybel to reperceive bonding
    pybel_mol = pb.readstring('xyz', MolToXYZBlock(stk_mol.to_rdkit_mol()))
    return pybel_mol_to_stk_mol(pybel_mol)

def num_connectivity_differences(stk_mol_1, stk_mol_2):
    def bond_tuple(bond):
        return(tuple(sorted((bond.get_atom1().get_id(),  bond.get_atom2().get_id()))))
        
    bond_tuples_1 = {bond_tuple(bond) for bond in stk_mol_1.get_bonds()}
    bond_tuples_2 = {bond_tuple(bond) for bond in stk_mol_2.get_bonds()}
    # return the number of bonds that appear in one of the molecules but not both
    return len(bond_tuples_1 ^ bond_tuples_2)

def gen_confs_openbabel(stk_mol, config):
    pybel_mol = stk_mol_to_pybel_mol(stk_mol)
    cs = pb.ob.OBConformerSearch()
    # Setup arguments: OBMol, numConformers, numChildren, mutability, convergence
    cs.Setup(pybel_mol.OBMol, config.initial_conformers, 5, 5, 5) 
    cs.Search()
    cs.GetConformers(pybel_mol.OBMol)
    
    stk_conformers = []
    for i in range(pybel_mol.OBMol.NumConformers()):
        pybel_mol.OBMol.SetConformer(i)
        stk_conformers.append(pybel_mol_to_stk_mol(pybel_mol))
    # preserve any additional stk information from original ligand in returned stk conformers
    return [stk_mol.with_position_matrix(stk_conformer.get_position_matrix())
            for stk_conformer in stk_conformers]
    
def gen_ligand_library_entry(stk_ligand, config):
    stk_conformers = gen_confs_openbabel(stk_ligand, config)
    stk_list_to_xyz_file(stk_conformers, 'conformers_ligand_only.xyz')
    unoptimized_complexes = [bind_to_dimethyl_Pd(ligand) for ligand in stk_conformers]
    ConformerEnsembleOptimizer(unoptimized_complexes, config).optimize()
