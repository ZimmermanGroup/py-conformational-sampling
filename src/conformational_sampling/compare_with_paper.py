from pathlib import Path
import pickle

from rdkit.Chem.rdMolAlign import AlignMol
from openbabel import pybel as pb

from conformational_sampling.utils import pybel_mol_to_rdkit_mol, pybel_mol_to_stk_mol

xyz_path = Path('/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/tests/test_data/L8_RE_Ra_AntiEndo_ts.xyz')
paper_pybel_mol = next(pb.readfile('xyz', str(xyz_path)))
paper_rdkit_mol = pybel_mol_to_rdkit_mol(paper_pybel_mol)

pickle_path = Path.home() / 'mols.pkl'
with pickle_path.open('rb') as file:
    mols = pickle.load(file)

for mol_name, mol_confs in mols.items():
    if mol_name == 'ligand_l8':
        for conf_idx, conformer in mol_confs.items():
            
            rmsd = AlignMol(conformer.ts_rdkit_mol, paper_rdkit_mol)
            print(rmsd)
            pass