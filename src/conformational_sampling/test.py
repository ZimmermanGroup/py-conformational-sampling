# %%
%reload_ext autoreload
%autoreload 2
from IPython.display import display

from pathlib import Path
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock
import openbabel as ob
from openbabel import pybel as pb
import nglview

from conformational_sampling.mixed_square_planar import MixedSquarePlanar

vitek_dmpp_ligand = Path('/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/examples/vitek_dmpp/ligand.xyz')
# vitek_dmpp_ligand = Path('/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/examples/bis_hydrazone/bis_hydrazone_ligand.xyz')
vitek_dmpp_ligand = next(pb.readfile('xyz', str(vitek_dmpp_ligand)))
vitek_dmpp_ligand = MolFromMolBlock(vitek_dmpp_ligand.write('mol'), removeHs=False)
nglview.show_rdkit(vitek_dmpp_ligand)


# %%
from rdkit import Chem
import stk
import stko

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

    return stk.ConstructedMolecule(
        topology_graph=MixedSquarePlanar(
            metals=metal,
            ligands={
                ligand_stk: (0, ),
                methyl: (1, 2),
            },
        ),
    )

complex = bind_to_dimethyl_Pd(vitek_dmpp_ligand)
mol = complex.to_rdkit_mol()
Chem.SanitizeMol(mol)
# display(nglview.show_rdkit(mol))
complex = stk.MCHammer().optimize(complex)
mol = complex.to_rdkit_mol()
Chem.SanitizeMol(mol)
# display(nglview.show_rdkit(mol))
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
