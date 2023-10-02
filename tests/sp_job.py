#!/export/zimmerman/soumikd/py-conformational-sampling/.venv/bin/python
#SBATCH -p zimintel --job-name=pygsm_singlePoint
#SBATCH -c14
#SBATCH --time=5-00:00:00
#SBATCH -o scratch/pystring_%a/sp_output.txt
#SBATCH -e scratch/pystring_%a/sp_error.txt
#SBATCH --array=0

import os
import pygsm
import sys
from pathlib import Path
sys.path.append(str(Path(pygsm.__file__).parent))

from xtb.ase.calculator import XTB
from conformational_sampling.config import Config
from pygsm.coordinate_systems import DelocalizedInternalCoordinates, PrimitiveInternalCoordinates, Topology
from pygsm.level_of_theories.ase import ASELoT
from pygsm.potential_energy_surfaces import PES
from pygsm.utilities import elements, manage_xyz, nifty
from pygsm.wrappers import Molecule
from conformational_sampling.analyze import ts_node, truncate_string_at_bond_formation

from rdkit.Chem import rdmolfiles

job_index = int(os.environ['SLURM_ARRAY_TASK_ID'])

path = Path.cwd() / f'scratch/pystring_{job_index}'
path.mkdir(exist_ok=True)
os.chdir(path)

# singlePoint configuration object for DFT
config3 = Config(
    # xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb',
    # ase_calculator=XTB(),
    # max_dft_opt_steps=2,
    num_cpus=14,
    dft_cpus_per_opt=14,
)

# qchem ase calculator setup
from ase.calculators.qchem import QChem
os.environ['QCSCRATCH'] = os.environ['SLURM_LOCAL_SCRATCH']
config3.ase_calculator = QChem(
    # method='wB97X-D',
    method='PBE',
    # basis='6-31G',
    # basis='STO-3G',
    # basis='gen',
    basis='LANL2DZ',
    # ecp='gen',
    ecp='fit-LANL2DZ',
    SCF_CONVERGENCE='6',
    nt=config3.dft_cpus_per_opt,
    SCF_MAX_CYCLES='750',
    SCF_ALGORITHM='RCA_DIIS',
    THRESH_RCA_SWITCH='4',
    # basisfile= '../../basis',
    # ecpfile= '../../ecp',
)

def sp_from_xyz():

    geoms = manage_xyz.read_xyzs('opt_converged_001.xyz')
    energies = []

    # class System:
    #     reductive_elim_torsion: tuple
    #     pro_dis_torsion: tuple

    # systems = {
    # 'ligand_l1': System(reductive_elim_torsion=(56, 55, 79, 78), pro_dis_torsion=(21, 11, 55, 65)),
    # 'ligand_l8': System(reductive_elim_torsion=(74, 73, 97, 96), pro_dis_torsion=(47, 9, 73, 83)),
    # 'ligand_achiral': System(reductive_elim_torsion=(36, 35, 59, 58), pro_dis_torsion=(21, 11, 35, 45)),
    # }

    # system = systems['ligand_achiral']

    # geoms_rdkit = [rdmolfiles.MolFromXYZBlock(geom) for geom in geoms]
    # truncated_geoms_rdkit = truncate_string_at_bond_formation(
    #         geoms_rdkit,
    #         *system.reductive_elim_torsion[1:3]
    #     )
    # for mol in truncated_geoms_rdkit:
    #     rdmolfiles.MolToXYZFile(mol, 'opt_converged_001_mod.xyz')


    with open('opt_converged_001.xyz') as f:
        lines = f.readlines()
        natoms = int(lines[0])
        total_lines = len(lines)
        num_geoms = total_lines/(natoms+2)

        se = 1
        for i in range(int(num_geoms)):
            energies.append(float(lines[se]))
            se = se + natoms + 2
        
    ts_node_energy = ts_node(energies)[1]
    ts_node_geom = geoms[energies.index(ts_node_energy)]
    reactant_geom = geoms[0]

    lot = ASELoT.from_options(config3.ase_calculator, geom=reactant_geom)

    pes = PES.from_options(lot=lot, ad_idx=0, multiplicity=1)
        
    atom_symbols = manage_xyz.get_atoms(reactant_geom)

    ELEMENT_TABLE = elements.ElementData()
    atoms = [ELEMENT_TABLE.from_symbol(atom) for atom in atom_symbols]
    xyz1 = manage_xyz.xyz_to_np(reactant_geom)
    xyz2 = manage_xyz.xyz_to_np(ts_node_geom)

    top1 = Topology.build_topology(
        xyz1,
        atoms,
    )

    top2 = Topology.build_topology(
            xyz2,
            atoms,
        )

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

    reactant = Molecule.from_options(
            geom=reactant_geom,
            PES=pes,
            coord_obj=coord_obj1,
            Form_Hessian=True,
        )

    transition_state =  Molecule.from_options(
            geom=ts_node_geom,
            PES=pes,
            coord_obj=coord_obj1,
            Form_Hessian=True,
        )

    nifty.printcool("Reactant singlePoint energy is {:5.4f} kcal/mol".format(reactant.energy))

    nifty.printcool("TS singlePoint energy is {:5.4f} kcal/mol".format(transition_state.energy))

def sp_from_string_node():
     pass

if __name__ == '__main__':
    sp_from_xyz()
