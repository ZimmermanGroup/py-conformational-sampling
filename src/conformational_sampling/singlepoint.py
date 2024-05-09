from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import re

from xtb.ase import calculator
from rdkit.Chem.rdmolfiles import MolToXYZBlock
from conformational_sampling.analyze import systems


from conformational_sampling.main import load_stk_mol_list
from conformational_sampling.utils import num_cpus, stk_mol_to_ase_atoms


def run_singlepoints(string_path):
    stk_string = load_stk_mol_list(string_path)
    energies = []
    for stk_mol in stk_string:
        ase_mol = stk_mol_to_ase_atoms(stk_mol)
        ase_mol.calc = calculator.XTB(solvent='acetonitrile')
        energy = ase_mol.get_potential_energy()
        print(f'{string_path = }, {energy = }')
        energies.append(energy)
    
    # get the conformer index for this string
    # conf_idx = int(re.search(r"pystring_(\d+)", str(string_path)).groups()[0])
    # out_path = mol_path / f'scratch/pystring_{conf_idx}/xtb_singlepoints_001.xyz'
    out_path = string_path.parent / 'xtb_singlepoints_001.xyz'
    with open(out_path, 'w') as file:
        for i, (stk_mol, energy) in enumerate(zip(stk_string, energies)):
            rdkit_mol = stk_mol.to_rdkit_mol()
            if energy is not None:
                rdkit_mol.SetProp('_Name', str(energy))
            file.write(MolToXYZBlock(rdkit_mol))


def singlepoints_conformer_ensemble():
    mol_paths = [system.mol_path for system in systems.values()]
    for mol_path in mol_paths:
        # mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l1_symm_xtb_crest')
        # mol_path = Path('/export/zimmerman/joshkamm/Lilly/test_xtb/example_l1_symm_xtb_crest')
        string_paths = tuple(mol_path.glob('scratch/pystring_*/opt_converged_001.xyz'))
        print(string_paths[:4])
        with ProcessPoolExecutor(max_workers=num_cpus()) as executor:
            executor.map(run_singlepoints, string_paths[:4])

if __name__ == '__main__':
    singlepoints_conformer_ensemble()

