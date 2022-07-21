from math import cos, pi, sin
from stk.molecular.topology_graphs import Edge
from stk.molecular.topology_graphs.metal_complex import MetalComplex
from stk.molecular.topology_graphs.metal_complex.vertices import BiDentateLigandVertex, MonoDentateLigandVertex, MetalVertex

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
