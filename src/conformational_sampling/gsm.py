import os
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

import pyGSM

from conformational_sampling.analyze import ts_node
from conformational_sampling.main import load_stk_mol_list

# workaround for issue with pyGSM installation
sys.path.append(str(Path(pyGSM.__file__).parent))
import numpy as np
from openbabel import pybel as pb
import stk
from pyGSM.coordinate_systems import (
    DelocalizedInternalCoordinates,
    Distance,
    PrimitiveInternalCoordinates,
    Topology,
)
from pyGSM.growing_string_methods import DE_GSM, SE_GSM
from pyGSM.level_of_theories.ase import ASELoT
from pyGSM.molecule import Molecule
from pyGSM.optimizers import eigenvector_follow
from pyGSM.potential_energy_surfaces import PES
from pyGSM.utilities import elements, manage_xyz, nifty
from pyGSM.utilities.cli_utils import get_driving_coord_prim
from pyGSM.utilities.cli_utils import plot as gsm_plot
from xtb.ase import calculator

from conformational_sampling.config import Config

# from conformational_sampling.analyze import ts_node

OPT_STEPS = 50 # 10 for debugging, 50 for production
TS_OPT_FILENAME = 'ts_opt.xyz'

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


def stk_se_de_gsm(
    path: Path, stk_mol: stk.Molecule, driving_coordinates, config: Config
):
    """Run pyGSM (changes directory, so run in its own process)"""
    path.mkdir(parents=True, exist_ok=True)
    os.chdir(path)
    stk_se_gsm(
        stk_mol=stk_mol,
        driving_coordinates=driving_coordinates,
        config=config,
    )
    stk_de_gsm(config=config)


def stk_se_de_gsm_single_node_parallel(stk_mols, driving_coordinates, config: Config):
    paths = [Path.cwd() / f'scratch/pystring_{i}' for i in range(len(stk_mols))]
    with ProcessPoolExecutor(max_workers=config.num_cpus) as executor:
        executor.map(
            stk_se_de_gsm, paths, stk_mols, repeat(driving_coordinates), repeat(config)
        )
    ts_opt_paths = [path / TS_OPT_FILENAME for path in paths]
    ts_opts = [
        list(pb.readfile('xyz', str(ts_opt_path))) for ts_opt_path in ts_opt_paths
    ]
    transition_states = [ts_opt[-1] for ts_opt in ts_opts]
    # JOSH - get string energies like in Conformer.__post_init__ and include those in the output file
    # should probably sort by them too
    # maybe I can just leave them as
    # manage_xyz.write_std_multixyz('transition_states.xyz')


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
        DQMAG_MAX=0.5,  # default value is 0.8
        ADD_NODE_TOL=0.01,  # default value is 0.1
        CONV_TOL=0.0005,
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
    
    
def stk_se_gsm(stk_mol: stk.Molecule, driving_coordinates, config: Config):
    
    nifty.printcool(" Building the LOT")
    atoms, xyz, geom = stk_mol_to_gsm_objects(stk_mol)
    ase_calculator = config.ase_calculator
    if ase_calculator is None:
        ase_calculator = calculator.XTB()
    lot = ASELoT.from_options(ase_calculator, geom=geom)
    
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
        opt_steps=OPT_STEPS,
        # path=path
    )

    se_gsm = SE_GSM.from_options(
        reactant=reactant,
        nnodes=20,
        optimizer=optimizer,
        xyz_writer=manage_xyz.write_std_multixyz,
        driving_coords=driving_coordinates,
        DQMAG_MAX=0.5,  # default value is 0.8
        ADD_NODE_TOL=0.01,  # default value is 0.1
        CONV_TOL=0.0005,
    )

    se_gsm.set_V0()

    se_gsm.nodes[0].gradrms = 0.0
    se_gsm.nodes[0].V0 = se_gsm.nodes[0].energy
    print(" Initial energy is %1.4f" % se_gsm.nodes[0].energy)
    se_gsm.add_GSM_nodeR()
    se_gsm.grow_string(max_iters=50, max_opt_steps=10)
    if se_gsm.tscontinue:
        se_gsm.pastts = se_gsm.past_ts()
        print("pastts {}".format(se_gsm.pastts))
        try:
            if se_gsm.pastts == 1:  # normal over the hill
                se_gsm.add_GSM_nodeR(1)
                se_gsm.add_last_node(2)
            elif se_gsm.pastts == 2 or se_gsm.pastts == 3:  # when cgrad is positive
                se_gsm.add_last_node(1)
                if se_gsm.nodes[se_gsm.nR-1].gradrms > 5.*se_gsm.options['CONV_TOL']:
                    se_gsm.add_last_node(1)
            elif se_gsm.pastts == 3:  # product detected by bonding
                se_gsm.add_last_node(1)
        except:
            print("Failed to add last node, continuing.")
            # probably need to make sure last node is optimized

    se_gsm.nnodes = se_gsm.nR
    se_gsm.nodes = se_gsm.nodes[:se_gsm.nR]
    energies = se_gsm.energies

    if se_gsm.TSnode == se_gsm.nR - 1:
        print(" The highest energy node is the last")
        print(" not continuing with TS optimization.")
        se_gsm.tscontinue = False

    print(" Number of nodes is ", se_gsm.nnodes)
    print(" Warning last node still not optimized fully")
    se_gsm.xyz_writer(
        "grown_string_{:03}.xyz".format(se_gsm.ID),
        se_gsm.geometries,
        se_gsm.energies,
        se_gsm.gradrmss,
        se_gsm.dEs,
    )
    print(" SSM growth phase over")
    se_gsm.done_growing = True

    # product = se_gsm.nodes[se_gsm.nnodes-1]

    # nifty.printcool("OPTIMIZING PRODUCT GEOMETRY")
    # optimizer.optimize(
    #         molecule=product,
    #         refE=reactant.energy,
    #         opt_steps=50,
    #         # path=path
    #     )


def stk_de_gsm(config: Config):
    # geoms = manage_xyz.read_xyzs('opt_converged_001.xyz')
    geoms = manage_xyz.read_xyzs("grown_string_000.xyz")
    ase_calculator = config.ase_calculator
    if ase_calculator is None:
        ase_calculator = calculator.XTB()
    lot = ASELoT.from_options(ase_calculator, geom=geoms[0])

    pes = PES.from_options(lot=lot, ad_idx=0, multiplicity=1)

    atom_symbols = manage_xyz.get_atoms(geoms[0])

    ELEMENT_TABLE = elements.ElementData()
    atoms = [ELEMENT_TABLE.from_symbol(atom) for atom in atom_symbols]
    xyz1 = manage_xyz.xyz_to_np(geoms[0])
    xyz2 = manage_xyz.xyz_to_np(geoms[-1])

    top1 = Topology.build_topology(
        xyz1,
        atoms,
    )

    # find union bonds
    xyz2 = manage_xyz.xyz_to_np(geoms[-1])
    top2 = Topology.build_topology(
        xyz2,
        atoms,
    )

    # Add bonds to top1 that are present in top2
    # It's not clear if we should form the topology so the bonds
    # are the same since this might affect the Primitives of the xyz1 (slightly)
    # Later we stil need to form the union of bonds, angles and torsions
    # However, I think this is important, the way its formulated, for identifiyin
    # the number of fragments and blocks, which is used in hybrid TRIC.
    for bond in top2.edges():
        if bond in top1.edges:
            pass
        elif (bond[1], bond[0]) in top1.edges():
            pass
        else:
            print(" Adding bond {} to top1".format(bond))
            if bond[0] > bond[1]:
                top1.add_edge(bond[0], bond[1])
            else:
                top1.add_edge(bond[1], bond[0])

    addtr = True
    connect = addcart = False
    p1 = PrimitiveInternalCoordinates.from_options(
        xyz=xyz1,
        atoms=atoms,
        connect=connect,
        addtr=addtr,
        addcart=addcart,
        topology=top1,
    )

    p2 = PrimitiveInternalCoordinates.from_options(
        xyz=xyz2,
        atoms=atoms,
        addtr=addtr,
        addcart=addcart,
        connect=connect,
        topology=top1,  # Use the topology of 1 because we fixed it above
    )

    p1.add_union_primitives(p2)

    coord_obj1 = DelocalizedInternalCoordinates.from_options(
        xyz=xyz1,
        atoms=atoms,
        addtr=addtr,
        addcart=addcart,
        connect=connect,
        primitives=p1,
    )

    # coord_obj2 = DelocalizedInternalCoordinates.from_options(
    #     xyz=xyz2,
    #     atoms=atoms,
    #     addtr=addtr,
    #     addcart=addcart,
    #     connect=connect,
    #     primitives=p2,
    # )

    reactant = Molecule.from_options(
        geom=geoms[0],
        PES=pes,
        coord_obj=coord_obj1,
        Form_Hessian=True,
    )

    product = Molecule.copy_from_options(
        reactant,
        xyz=xyz2,
        new_node_id=len(geoms) - 1,
        copy_wavefunction=False,
    )

    optimizer = eigenvector_follow.from_options(
        Linesearch="backtrack",
        OPTTHRESH=0.0005,
        DMAX=0.5,
        abs_max_step=0.5,
        conv_Ediff=0.1,
    )

    nifty.printcool("OPTIMIZING REACTANT GEOMETRY")
    optimizer.optimize(
        molecule=reactant,
        refE=reactant.energy,
        opt_steps=OPT_STEPS,
    )

    nifty.printcool("OPTIMIZING PRODUCT GEOMETRY")
    optimizer.optimize(
        molecule=product,
        refE=reactant.energy,
        opt_steps=OPT_STEPS,
    )

    # For xTB

    de_gsm_max_nodes = 20
    path = Path.cwd() / "scratch/001"
    path.mkdir(exist_ok=True)

    for item in range(de_gsm_max_nodes):
        path = Path.cwd() / f"scratch/001/{item}"
        path.mkdir(exist_ok=True)

    de_gsm = DE_GSM.from_options(
        reactant=reactant,
        product=product,
        nnodes=15,
        optimizer=optimizer,
        xyz_writer=manage_xyz.write_std_multixyz,
        ID=1,
    )

    # For DFT

    # de_gsm_max_nodes = 20
    # path = Path.cwd() / 'scratch/002'
    # path.mkdir(exist_ok=True)

    # for item in range(de_gsm_max_nodes):
    #     path = Path.cwd() / f'scratch/002/{item}'
    #     path.mkdir(exist_ok=True)

    # de_gsm = DE_GSM.from_options(
    #     reactant=reactant,
    #     product=product,
    #     nnodes=15,
    #     optimizer=optimizer,
    #     xyz_writer=manage_xyz.write_std_multixyz,
    #     ID=2,
    # )

    de_gsm.go_gsm()

    # TS-Optimization following DE-GSM run

    ts_node_energy = ts_node(de_gsm.energies)[1]
    ts_node_index = de_gsm.energies.index(ts_node_energy)
    ts_node_geom = de_gsm.nodes[ts_node_index]

    nifty.printcool('Optimizing TS node')
    geoms, energies = optimizer.optimize(
        molecule=ts_node_geom,
        refE=de_gsm.energies[0],
        opt_steps=OPT_STEPS,
        opt_type='TS',
        ictan=de_gsm.ictan[ts_node_index],
    )
    manage_xyz.write_std_multixyz(
        TS_OPT_FILENAME, geoms, energies, gradrms=None, dEs=None
    )
    de_gsm.nodes[ts_node_index] = (
        ts_node_geom  # has ts_node_geom actually been updated during optimize?
    )

    gsm_plot(de_gsm.energies, x=range(len(de_gsm.energies)), title=1)