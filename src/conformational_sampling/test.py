# %%
from stk.molecular.topology_graphs import Edge
from stk.molecular.topology_graphs.metal_complex import MetalComplex
from stk.molecular.topology_graphs.metal_complex.vertices import BiDentateLigandVertex, MonoDentateLigandVertex, MetalVertex

class CustomSquarePlanar(MetalComplex):

    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )
    _ligand_vertex_prototypes = (
        BiDentateLigandVertex(1, (2.5, 2.5, 0)),
        MonoDentateLigandVertex(2, (-2.5, 0, 0)),
        MonoDentateLigandVertex(3, (0, -2.5, 0)),
    )
    
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(2.5, 0, 0),
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(0, 2.5, 0),
        ),

        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
            position=(-1, 0, 0),
        ),
        Edge(
            id=3,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
            position=(0, -1, 0),
        ),
    )

# %%
%reload_ext autoreload
%autoreload 2
from IPython.display import display

from pathlib import Path
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock
import openbabel as ob
from openbabel import pybel as pb
import nglview

vitek_dmpp_ligand = Path('/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/examples/vitek_dmpp/ligand.xyz')
target = Path('/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/examples/vitek_dmpp/target.xyz')
vitek_dmpp_ligand = next(pb.readfile('xyz', str(vitek_dmpp_ligand)))
target = next(pb.readfile('xyz', str(target)))
vitek_dmpp_ligand = MolFromMolBlock(vitek_dmpp_ligand.write('mol'), removeHs=False)
target = MolFromMolBlock(target.write('mol'), removeHs=False)
print(MolToMolBlock(target))
nglview.show_rdkit(vitek_dmpp_ligand)


# %%
from rdkit import Chem
import stk
import stko

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

Chem.rdmolops.Kekulize(vitek_dmpp_ligand)
vitek_dmpp_ligand_stk = stk.BuildingBlock.init_from_rdkit_mol(
    vitek_dmpp_ligand,
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='P',
            bonders=(0,),
            deleters=(),
        )
    ]
)

complex = stk.ConstructedMolecule(
    topology_graph=CustomSquarePlanar(
        metals=metal,
        ligands={
            vitek_dmpp_ligand_stk: (0, ),
            methyl: (1, 2),
        },
    ),
)

mol = complex.to_rdkit_mol()
Chem.SanitizeMol(mol)
display(nglview.show_rdkit(mol))
complex = stk.MCHammer().optimize(complex)
mol = complex.to_rdkit_mol()
Chem.SanitizeMol(mol)
display(nglview.show_rdkit(mol))
complex = stko.MetalOptimizer().optimize(complex)
mol = complex.to_rdkit_mol()
Chem.SanitizeMol(mol)
display(nglview.show_rdkit(mol))
complex = stko.XTB('/export/apps/CentOS7/xtb/xtb/bin/xtb').optimize(complex)
mol = complex.to_rdkit_mol()
Chem.SanitizeMol(mol)
display(nglview.show_rdkit(mol))
