#!/export/zimmerman/joshkamm/apps/mambaforge/envs/conformational-sampling/bin/python
#SBATCH -p guest --job-name=conformational_sampling
#SBATCH -c16
#SBATCH -o output.txt

from pathlib import Path
import stk
from conformational_sampling.main import load_stk_mol, gen_ligand_library_entry, pybel_mol_to_stk_mol, stk_list_to_xyz_file
from conformational_sampling.config import Config

ligand_path = Path('ligand.xyz') # name of file containing ligand geometry
stk_ligand = load_stk_mol(ligand_path)

# specifies atoms of the ligand that bind to the metal, in this case as a smarts string
functional_group_factory = stk.SmartsFunctionalGroupFactory(
    smarts='P',
    bonders=(0,),
    deleters=(),
)
stk_ligand = stk.BuildingBlock.init_from_molecule(stk_ligand, functional_groups=[functional_group_factory])

# py-conformational-sampling configuration object
config = Config(initial_conformers=100, xtb_path='/export/apps/CentOS7/xtb/xtb/bin/xtb')

# generates conformers, performs multiple step optimization and uniqueness filtering
gen_ligand_library_entry(stk_ligand, config)

# EXPERIMENTING WITH RDKIT CONFORMER GENERATION
# rdkit_mol = stk_ligand.to_rdkit_mol()
# from rdkit import Chem
# from rdkit.Chem import AllChem
# print(Chem.FindMolChiralCenters(rdkit_mol, includeUnassigned=True, useLegacyImplementation=False))
# # AllChem.AssignStereochemistryFrom3D(rdkit_mol, confId=0, replaceExistingTags=True)
# si = Chem.FindPotentialStereo(rdkit_mol)
# for element in si:
#     print(f'  Type: {element.type}, Which: {element.centeredOn}, Specified: {element.specified}, Descriptor: {element.descriptor} ')
# pass

# EXPERIMENTING WITH OPEN BABEL CONFORMER GENERATION
# from openbabel import pybel as pb
# pybel_mol = next(pb.readfile('xyz', str(ligand_path)))
# # print(len(pybel_mol.conformers))
# cs = pb.ob.OBConformerSearch()
# cs.Setup(pybel_mol.OBMol, 10, 5, 5, 5) # numConformers, numChildren, mutability, convergence
# cs.Search()
# cs.GetConformers(pybel_mol.OBMol)
# print(pybel_mol.OBMol.NumConformers())
# stk_conformers = []
# for i in range(pybel_mol.OBMol.NumConformers()):
#     pybel_mol.OBMol.SetConformer(i)
#     stk_conformers.append(pybel_mol_to_stk_mol(pybel_mol))
# stk_list_to_xyz_file(stk_conformers, 'test_pybel_generation.xyz')