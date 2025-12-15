"""
Analyze GSM results for the cationic 2-aza-Cope rearrangement.

This script extracts energy profiles from the GSM calculations
and identifies the barrier heights and reaction energies.
"""
import numpy as np
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')

def read_energies_from_opt_converged(gsm_dir):
    """Read energies from opt_converged_001.xyz file."""
    opt_file = gsm_dir / 'opt_converged_001.xyz'
    if not opt_file.exists():
        return None
    
    with open(opt_file, 'r') as f:
        lines = f.readlines()
    
    # Extract energies from XYZ file comment lines
    energies = []
    natoms = int(lines[0].strip())
    step = natoms + 2
    
    for j in range(0, len(lines), step):
        if j + 1 < len(lines):
            comment = lines[j + 1].strip()
            # Energy is just a number on the comment line
            try:
                energy = float(comment)
                energies.append(energy)
            except ValueError:
                # Try to extract from formatted string
                parts = comment.split()
                if len(parts) > 0:
                    try:
                        energies.append(float(parts[0]))
                    except (ValueError, IndexError):
                        pass
    
    return np.array(energies) if energies else None


def analyze_conformer(conformer_idx):
    """Analyze GSM results for a single conformer."""
    gsm_dir = Path(f'scratch/pystring_{conformer_idx}')
    
    if not gsm_dir.exists():
        logging.warning(f'Directory {gsm_dir} does not exist')
        return None
    
    energies = read_energies_from_opt_converged(gsm_dir)
    
    if energies is None or len(energies) == 0:
        logging.warning(f'No energies found for conformer {conformer_idx}')
        return None
    
    # Convert to relative energies in kcal/mol
    rel_energies = (energies - energies[0]) * 627.5095
    
    # Find barrier and product
    max_energy = np.max(rel_energies)
    max_idx = np.argmax(rel_energies)
    final_energy = rel_energies[-1]
    
    return {
        'conformer': conformer_idx,
        'n_nodes': len(energies),
        'start_energy': energies[0],
        'barrier_height': max_energy,
        'ts_node': max_idx,
        'reaction_energy': final_energy,
        'final_energy': energies[-1],
        'profile': rel_energies
    }


if __name__ == '__main__':
    logging.info('\n' + '='*70)
    logging.info('GSM Results Analysis: Cationic 2-aza-Cope Rearrangement')
    logging.info('='*70)
    
    results = []
    for i in range(3):
        result = analyze_conformer(i)
        if result:
            results.append(result)
    
    if not results:
        logging.error('No results found!')
        exit(1)
    
    # Print summary
    logging.info(f'\nAnalyzed {len(results)} conformer(s)\n')
    
    for result in results:
        logging.info(f"Conformer {result['conformer']}:")
        logging.info(f"  Nodes along pathway: {result['n_nodes']}")
        logging.info(f"  Starting energy: {result['start_energy']:.6f} Ha")
        logging.info(f"  Barrier height: {result['barrier_height']:.2f} kcal/mol")
        logging.info(f"  TS at node: {result['ts_node']}")
        logging.info(f"  Final energy: {result['final_energy']:.6f} Ha")
        logging.info(f"  Reaction energy: {result['reaction_energy']:.2f} kcal/mol")
        logging.info('')
    
    # Find lowest barrier
    barriers = [r['barrier_height'] for r in results]
    lowest_idx = np.argmin(barriers)
    
    logging.info('='*70)
    logging.info(f'Lowest barrier pathway: Conformer {results[lowest_idx]["conformer"]}')
    logging.info(f'Barrier: {results[lowest_idx]["barrier_height"]:.2f} kcal/mol')
    logging.info('='*70)
    
    # Print energy profiles
    logging.info('\nEnergy Profiles (kcal/mol):')
    logging.info('-'*70)
    logging.info(f'{"Node":>6}  {"Conf 0":>10}  {"Conf 1":>10}  {"Conf 2":>10}')
    logging.info('-'*70)
    
    max_nodes = max(r['n_nodes'] for r in results)
    for i in range(max_nodes):
        row = f'{i:>6}'
        for result in results:
            if i < result['n_nodes']:
                row += f'  {result["profile"][i]:>10.2f}'
            else:
                row += f'  {"---":>10}'
        logging.info(row)
    
    logging.info('='*70)
