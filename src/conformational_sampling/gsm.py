from pathlib import Path
import sys
import pyGSM

from conformational_sampling.main import load_stk_mol_list
# workaround for issue with pyGSM installation
sys.path.append(str(Path(pyGSM.__file__).parent))
import ase.io
import numpy as np
from ase.calculators.morse import MorsePotential

from pyGSM.coordinate_systems import DelocalizedInternalCoordinates, PrimitiveInternalCoordinates, Distance, Topology
from pyGSM.level_of_theories.ase import ASELoT
from pyGSM.optimizers import eigenvector_follow
from pyGSM.potential_energy_surfaces import PES
from pyGSM.utilities import elements, manage_xyz, nifty
from pyGSM.molecule import Molecule
from pyGSM.growing_string_methods import SE_GSM, DE_GSM
from pyGSM.utilities.cli_utils import plot as gsm_plot
from pyGSM.utilities.cli_utils import get_driving_coord_prim
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

    # testing whether adding any edges near to a transition metal to the topology will help with
    # internal coordinate stability
    # HARDCODED TO SPECIFIC TEST MOLECULE
    metal_connections = ((0,9), (0,56), (0,73), (0,97), (0,83), (10,73), (10,97), (56,73), (56,97))
    for metal_connection in metal_connections:
        if metal_connection in top.edges:
            pass
        elif tuple(reversed(metal_connection)) in top.edges():
            pass
        else:
            print(" Adding bond {} to top1".format(metal_connection))
            top.add_edge(*metal_connection)

    nifty.printcool("Building Primitive Internal Coordinates")
    p1 = PrimitiveInternalCoordinates.from_options(
        xyz=xyz,
        atoms=atoms,
        addtr=True,  # Add TRIC
        topology=top,
    )

    nifty.printcool("Building Delocalized Internal Coordinates")
    coord_obj1 = DelocalizedInternalCoordinates.from_options(
        xyz=xyz,
        atoms=atoms,
        addtr=True,  # Add TRIC
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
    optimizer = eigenvector_follow.from_options(Linesearch='backtrack', OPTTHRESH=0.0005, DMAX=0.5, abs_max_step=0.5,
                                                conv_Ediff=0.1)

    nifty.printcool("initial energy is {:5.4f} kcal/mol".format(reactant.energy))

    nifty.printcool("REACTANT GEOMETRY NOT FIXED!!! OPTIMIZING")
    optimizer.optimize(
            molecule=reactant,
            refE=reactant.energy,
            opt_steps=50,
            # path=path
        )

    se_gsm = SE_GSM.from_options(
        reactant=reactant,
        product=product if config.restart_gsm else None,
        nnodes=len(geoms) if config.restart_gsm else 20,
        optimizer=optimizer,
        xyz_writer=manage_xyz.write_std_multixyz,
        driving_coords=driving_coordinates,
        DQMAG_MAX=0.5, #default value is 0.8
        ADD_NODE_TOL=0.01, #default value is 0.1
        CONV_TOL = 0.0005,
    )
    
    # run pyGSM, setting up restart if necessary
    if config.restart_gsm:
        se_gsm.setup_from_geometries(geoms, reparametrize=True, restart_energies=False)
    se_gsm.go_gsm()
    gsm_plot(se_gsm.energies, x=range(len(se_gsm.energies)), title=0)
    
def stk_gsm_command_line(stk_mol: stk.Molecule, driving_coordinates, config: Config):
    # guessing a command line simulation
    import sys
    stk_mol.write('initial0001.xyz') #### to xyz
    with open('isomers0001.txt', 'w') as f:
        f.write('ADD 56 80\n BREAK 1 56\n BREAK 1 80\n')

    sys.argv = ["gsm",
        "-coordinate_type", "TRIC",
        "-xyzfile", "initial0001.xyz",
        "-mode", "SE_GSM",
        "-package", "ase",
        "--ase-class", "ase.calculators.qchem.QChem",
        "--ase-kwargs", '{"method":"PBE", "basis":"LANL2DZ", "ecp":"fit-LANL2DZ", "SCF_CONVERGENCE": "5", "nt": 8, "SCF_MAX_CYCLES": "200", "SCF_ALGORITHM":"DIIS"}',
        "-DQMAG_MAX", "0.6",
        "-num_nodes", "15",
        "-isomers", "isomers0001.txt",
    ]
    # main()
    raise NotImplementedError(
        'the way to access the pyGSM command line functionality from python has changed'
    )
    
    
def stk_se_de_gsm(stk_mol: stk.Molecule, driving_coordinates, config: Config):
    
    nifty.printcool(" Building the LOT")
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
        addtr=True,  # Add TRIC
        topology=top,
    )

    nifty.printcool("Building Delocalized Internal Coordinates")
    coord_obj1 = DelocalizedInternalCoordinates.from_options(
        xyz=xyz,
        atoms=atoms,
        addtr=True,  # Add TRIC
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
                                                conv_Ediff=0.1)

    nifty.printcool("initial energy is {:5.4f} kcal/mol".format(reactant.energy))

    nifty.printcool("REACTANT GEOMETRY NOT FIXED!!! OPTIMIZING")
    optimizer.optimize(
            molecule=reactant,
            refE=reactant.energy,
            opt_steps=50,
            # path=path
        )
    
    se_gsm = SE_GSM.from_options(
        reactant=reactant,
        nnodes=20,
        optimizer=optimizer,
        xyz_writer=manage_xyz.write_std_multixyz,
        driving_coords=driving_coordinates,
        DQMAG_MAX=0.5, #default value is 0.8
        ADD_NODE_TOL=0.01, #default value is 0.1
        CONV_TOL = 0.0005,
    )

    se_gsm.set_V0()

    se_gsm.nodes[0].gradrms = 0.
    se_gsm.nodes[0].V0 = se_gsm.nodes[0].energy
    print(" Initial energy is %1.4f" % se_gsm.nodes[0].energy)
    se_gsm.add_GSM_nodeR()
    se_gsm.grow_string(max_iters=50, max_opt_steps=10)
    if se_gsm.tscontinue:
        se_gsm.pastts = se_gsm.past_ts()
        print("pastts {}".format(se_gsm.pastts))
        try:
            if se_gsm.pastts == 1: #normal over the hill
                se_gsm.add_GSM_nodeR(1)
                se_gsm.add_last_node(2)
            elif se_gsm.pastts == 2 or se_gsm.pastts==3: #when cgrad is positive
                se_gsm.add_last_node(1)
                if se_gsm.nodes[se_gsm.nR-1].gradrms > 5.*se_gsm.options['CONV_TOL']:
                    se_gsm.add_last_node(1)
            elif se_gsm.pastts == 3: #product detected by bonding
                se_gsm.add_last_node(1)
        except:
            print("Failed to add last node, continuing.")
            # probably need to make sure last node is optimized

    se_gsm.nnodes = se_gsm.nR
    se_gsm.nodes = se_gsm.nodes[:se_gsm.nR]
    energies = se_gsm.energies

    if se_gsm.TSnode == se_gsm.nR-1:
        print(" The highest energy node is the last")
        print(" not continuing with TS optimization.")
        se_gsm.tscontinue = False

    print(" Number of nodes is ", se_gsm.nnodes)
    print(" Warning last node still not optimized fully")
    se_gsm.xyz_writer('grown_string_{:03}.xyz'.format(se_gsm.ID), se_gsm.geometries, se_gsm.energies, se_gsm.gradrmss, se_gsm.dEs)
    print(" SSM growth phase over")
    se_gsm.done_growing = True

    product = se_gsm.nodes[se_gsm.nnodes-1]

    nifty.printcool("OPTIMIZING PRODUCT GEOMETRY")
    optimizer.optimize(
            molecule=product,
            refE=reactant.energy,
            opt_steps=50,
            # path=path
        )
    
    de_gsm = DE_GSM.from_options(
        reactant=reactant,
        product=product,
        nnodes=15,
        optimizer=optimizer,
        xyz_writer=manage_xyz.write_std_multixyz,
        ID=1,
    )

    de_gsm.go_gsm()

    gsm_plot(de_gsm.energies, x=range(len(de_gsm.energies)), title=1)

