from pathlib import Path

from xtb.ase import calculator

from conformational_sampling.main import load_stk_mol_list
from conformational_sampling.utils import stk_mol_to_ase_atoms


def run_singlepoints():
    mol_path=Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l1_symm_xtb_crest')
    string_paths = tuple(mol_path.glob('scratch/pystring_*/opt_converged_001.xyz'))
    for string_path in string_paths[:2]:
        stk_string = load_stk_mol_list(string_path)
        for stk_mol in stk_string:
            ase_mol = stk_mol_to_ase_atoms(stk_mol)
            ase_mol.calc = calculator.XTB(solvent='acetonitrile')
            energy = ase_mol.get_potential_energy()
            print(f'{string_path = }, {energy = }')

if __name__ == '__main__':
    run_singlepoints()

