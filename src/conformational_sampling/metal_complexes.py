from dataclasses import dataclass
from math import cos, pi, sin

import stk
from stk.molecular.topology_graphs import Edge
from stk.molecular.topology_graphs.metal_complex import MetalComplex
from stk.molecular.topology_graphs.metal_complex.vertices import (
    BiDentateLigandVertex,
    MetalVertex,
    MonoDentateLigandVertex,
)

from conformational_sampling.main import bind_ligands
from conformational_sampling.utils import stk_metal


class TwoMonoOneBidentateSquarePlanar(MetalComplex):

    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )
    _ligand_vertex_prototypes = (
        # BiDentateLigandVertex(1, (2.5, 2.5, 0)),
        BiDentateLigandVertex(1, (5, 5, 0)),
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

class OneLargeTwoSmallMonodentateTrigonalPlanar(MetalComplex):

    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )
    _ligand_vertex_prototypes = (
        MonoDentateLigandVertex(1, (5, 0, 0)),
        MonoDentateLigandVertex(2, (2.5*cos(2*pi/3), 2.5*sin(2*pi/3), 0)),
        MonoDentateLigandVertex(3, (2.5*cos(4*pi/3), 2.5*sin(4*pi/3), 0)),
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
            vertex2=_ligand_vertex_prototypes[1],
            position=(cos(2*pi/3), sin(2*pi/3), 0),
        ),

        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
            position=(cos(4*pi/3), sin(4*pi/3), 0),
        ),
    )


@dataclass
class CatalyticReactionComplex:
    metal: stk.Molecule
    ancillary_ligand: stk.Molecule
    reactive_ligand_1: stk.Molecule
    reactive_ligand_2: stk.Molecule
    
    
    def __post_init__(self) -> None:
        'create a complex with ane ancillary ligand and two reactive ligands, like bind_ligands'
        self.complex = bind_ligands(
            self.metal,
            self.ancillary_ligand,
            self.reactive_ligand_1,
            self.reactive_ligand_2,
        )
    
    
    def gen_reductive_elim_drive_coords(self):
        'generate reductive elimination driving coordinates for this complex'
        atom_infos = list(self.complex.get_atom_infos())
        driving_coordinates = []

        def is_metal_reactive_ligand_bond(
            atom_info_1: stk.AtomInfo, atom_info_2: stk.AtomInfo
        ) -> bool:
            building_block_1 = atom_info_1.get_building_block()
            building_block_2 = atom_info_2.get_building_block()
            return building_block_1 is self.metal and building_block_2 in (
                self.reactive_ligand_1,
                self.reactive_ligand_2,
            )

        # determine driving coordinates to break
        for bond_info in self.complex.get_bond_infos():
            if bond_info.get_building_block() is None:
                bond = bond_info.get_bond()
                atom_id_1, atom_id_2 = bond.get_atom1().get_id(), bond.get_atom2().get_id()
                atom_info_1 = next(self.complex.get_atom_infos(atom_id_1))
                atom_info_2 = next(self.complex.get_atom_infos(atom_id_2))
                if (is_metal_reactive_ligand_bond(atom_info_1, atom_info_2)
                    or is_metal_reactive_ligand_bond(atom_info_2, atom_info_1)):
                    # this bond is broken in reductive elimination
                    driving_coordinates += ('BREAK', atom_id_1, atom_id_2)
                
        # CHECK THAT THERE ARE TWO BREAKING DRIVING COORDS

        # go from zero-indexed to one-indexed
        return [(type, i + 1, j + 1) for (type, i, j) in driving_coordinates]
