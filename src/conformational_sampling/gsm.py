import sys
# sys.path.append('/export/zimmerman/joshkamm/Lilly/pyGSM/pygsm')

import ase.io
import numpy as np
from ase.calculators.morse import MorsePotential

from pygsm.coordinate_systems import DelocalizedInternalCoordinates, PrimitiveInternalCoordinates, Topology
from pygsm.level_of_theories.ase import ASELoT
from pygsm.optimizers import eigenvector_follow
from pygsm.potential_energy_surfaces import PES
from pygsm.utilities import elements, manage_xyz, nifty
from pygsm.wrappers import Molecule
from pygsm.growing_string_methods import SE_GSM
from pygsm.wrappers.main import plot as gsm_plot
from pygsm.wrappers.main import get_driving_coord_prim, Distance
import stk

from conformational_sampling.config import Config

def stk_gsm(stk_mol: stk.Molecule, driving_coordinates, config: Config):
    nifty.printcool(" Building the LOT")
    # lot = ASELoT.from_options(MorsePotential(), geom=geom)
    ELEMENT_TABLE = elements.ElementData()
    atoms = [ELEMENT_TABLE.from_atomic_number(atom.get_atomic_number())
             for atom in stk_mol.get_atoms()]
    xyz = stk_mol.get_position_matrix()
    atom_symbols = np.array(atom.symbol for atom in atoms)
    geom = np.column_stack([atom_symbols, xyz]).tolist()
    lot = ASELoT.from_options(config.ase_calculator, geom)

    nifty.printcool(" Building the PES")
    pes = PES.from_options(
        lot=lot,
        ad_idx=0,
        multiplicity=1,
    )

    nifty.printcool("Building the topology")
    top = Topology.build_topology(
        xyz,
        atoms,
    )

    driving_coord_prims = []
    for dc in driving_coordinates:
        prim = get_driving_coord_prim(dc)
        if prim is not None:
            driving_coord_prims.append(prim)

    for prim in driving_coord_prims:
        if type(prim) == Distance:
            bond = (prim.atoms[0], prim.atoms[1])
            if bond in top.edges:
                pass
            elif (bond[1], bond[0]) in top.edges():
                pass
            else:
                print(" Adding bond {} to top1".format(bond))
                top.add_edge(bond[0], bond[1])

    nifty.printcool("Building Primitive Internal Coordinates")
    p1 = PrimitiveInternalCoordinates.from_options(
        xyz=xyz,
        atoms=atoms,
        addtr=False,  # Add TRIC
        topology=top,
    )

    nifty.printcool("Building Delocalized Internal Coordinates")
    coord_obj1 = DelocalizedInternalCoordinates.from_options(
        xyz=xyz,
        atoms=atoms,
        addtr=False,  # Add TRIC
        primitives=p1,
    )

    nifty.printcool("Building Molecule")
    reactant = Molecule.from_options(
        geom=geom,
        PES=pes,
        coord_obj=coord_obj1,
        Form_Hessian=True,
    )

    nifty.printcool("Creating optimizer")
    optimizer = eigenvector_follow.from_options(Linesearch='backtrack', OPTTHRESH=0.0005, DMAX=0.5, abs_max_step=0.5,
                                                conv_Ediff=0.5)

    nifty.printcool("initial energy is {:5.4f} kcal/mol".format(reactant.energy))
    # geoms, energies = optimizer.optimize(
    #     molecule=reactant,
    #     refE=reactant.energy,
    #     opt_steps=5,
    #     verbose=True,
    # )

    # nifty.printcool("Final energy is {:5.4f}".format(reactant.energy))
    # manage_xyz.write_xyz('minimized.xyz', geoms[-1], energies[-1], scale=1.)
    
    # trying out GSM
    se_gsm = SE_GSM.from_options(
        reactant=reactant,
        nnodes=7,
        optimizer=optimizer,
        xyz_writer=manage_xyz.write_std_multixyz,
        driving_coords=driving_coordinates,        
    )
    
    se_gsm.go_gsm()
    gsm_plot(se_gsm.energies, x=range(len(se_gsm.energies)), title=0)
