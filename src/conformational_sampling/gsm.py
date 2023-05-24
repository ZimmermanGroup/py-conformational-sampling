from pathlib import Path
import sys
import pygsm

from conformational_sampling.main import load_stk_mol_list
# workaround for issue with pygsm installation
sys.path.append(str(Path(pygsm.__file__).parent))
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

def stk_mol_to_gsm_objects(stk_mol: stk.Molecule):
    ELEMENT_TABLE = elements.ElementData()
    # atoms is a list of pygsm element objects
    atoms = [ELEMENT_TABLE.from_atomic_number(atom.get_atomic_number())
            for atom in stk_mol.get_atoms()]
    # xyz is a numpy array of the position matrix
    xyz = stk_mol.get_position_matrix()
    atom_symbols = np.array(list(atom.symbol for atom in atoms))

    # geom is an aggregate ndarray with the structure of the body of an XYZ file
    geom = np.column_stack([atom_symbols, xyz]).tolist()
    return atoms, xyz, geom


def stk_mol_list_to_gsm_objects(stk_mol_list):
    atoms, xyz, geom = stk_mol_to_gsm_objects(stk_mol_list[0])
    geoms = [stk_mol_to_gsm_objects(stk_mol)[2] for stk_mol in stk_mol_list]
    return atoms, xyz, geoms


def stk_gsm(stk_mol: stk.Molecule, driving_coordinates, config: Config):
    nifty.printcool(" Building the LOT")
    
    if config.restart_gsm:
        stk_string = load_stk_mol_list(config.restart_gsm)
        atoms, xyz, geoms = stk_mol_list_to_gsm_objects(stk_string)
        geom = geoms[0]
    else:
        atoms, xyz, geom = stk_mol_to_gsm_objects(stk_mol)
    
    lot = ASELoT.from_options(config.ase_calculator, geom=geom)

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
    if config.restart_gsm:
        product = Molecule.from_options(
            geom=geoms[-1],
            PES=pes,
            coord_obj=coord_obj1,
            Form_Hessian=True,
        )

    nifty.printcool("Creating optimizer")
    optimizer = eigenvector_follow.from_options(Linesearch='backtrack', OPTTHRESH=0.005, DMAX=0.5, abs_max_step=0.5,
                                                conv_Ediff=0.5)

    nifty.printcool("initial energy is {:5.4f} kcal/mol".format(reactant.energy))

    se_gsm = SE_GSM.from_options(
        reactant=reactant,
        product=product if config.restart_gsm else None,
        nnodes=len(geoms) if config.restart_gsm else 20,
        optimizer=optimizer,
        xyz_writer=manage_xyz.write_std_multixyz,
        driving_coords=driving_coordinates,
        DQMAG_MAX=0.6,
    )
    
    # run pyGSM, setting up restart if necessary
    if config.restart_gsm:
        se_gsm.setup_from_geometries(geoms, reparametrize=True, restart_energies=False)
    se_gsm.go_gsm()
    gsm_plot(se_gsm.energies, x=range(len(se_gsm.energies)), title=0)
