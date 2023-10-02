# %%
%reload_ext autoreload
%autoreload 2
from dataclasses import dataclass
import os
import re
from IPython.display import display
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolToPDBBlock
from rdkit.Chem.rdMolAlign import AlignMol
from rdkit.Chem import rdMolTransforms

import param
import hvplot.pandas # noqa
import holoviews as hv
from holoviews import opts
from holoviews.streams import Selection1D
import panel as pn
from panel_chemistry.pane import \
    NGLViewer  # panel_chemistry needs to be imported before you run pn.extension()
from panel_chemistry.pane.ngl_viewer import EXTENSIONS
pn.extension('bokeh')
pn.extension(comms='vscode')
pn.extension('tabulator')
pn.extension("ngl_viewer", sizing_mode="stretch_width")
from pathlib import Path
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock
import openbabel as ob
from openbabel import pybel as pb

from pygsm.utilities.units import KCAL_MOL_PER_AU
from conformational_sampling.analyze import ts_node, truncate_string_at_bond_formation
from conformational_sampling.utils import free_energy_diff

    
# setup_mol()
@dataclass
class System:
    reductive_elim_torsion: tuple
    pro_dis_torsion: tuple

systems = {
    'ligand_l1': System(reductive_elim_torsion=(56, 55, 79, 78), pro_dis_torsion=(21, 11, 55, 65)),
    'ligand_l8': System(reductive_elim_torsion=(74, 73, 97, 96), pro_dis_torsion=(47, 9, 73, 83)),
    'ligand_achiral': System(reductive_elim_torsion=(36, 35, 59, 58), pro_dis_torsion=(21, 11, 35, 45)),
}
ligand_name = 'ligand_achiral'
system = systems[ligand_name]

@dataclass
class Conformer:
    string_path: Path
    
    def __post_init__(self):
        string_nodes = list(pb.readfile('xyz', str(self.string_path)))
        self.string_nodes = [MolFromMolBlock(node.write('mol'), removeHs=False)
                             for node in string_nodes] # convert to rdkit
        
        # previous way of getting reactant energy
        # raw_dft_energy_path = self.string_path.parent / 'scratch' / '000' / 'E_0.txt'
        # raw_dft_energy_au = float(raw_dft_energy_path.read_text().split()[2])
        # self.string_energies = [(raw_dft_energy_au + float(MolToMolBlock(node).split()[0])) * KCAL_MOL_PER_AU
        #                         for node in self.string_nodes]
        
        # getting reactant energy
        raw_dft_energy_path = self.string_path.parent / 'de_output.txt'
        with open(raw_dft_energy_path) as file:
            while line := file.readline():
                if line.startswith(' Energy of the end points are'):
                    # extracting energy from output file formatting
                    raw_dft_energy_kcal_mol = float(line.split()[-2][:-1])
        self.string_energies = [raw_dft_energy_kcal_mol + float(MolToMolBlock(node).split()[0]) * KCAL_MOL_PER_AU
                                for node in self.string_nodes]
        
        self.truncated_string = truncate_string_at_bond_formation(
            self.string_nodes,
            *system.reductive_elim_torsion[1:3] # middle atoms of torsion
        )
        if not self.truncated_string:
            return
        max_diff, self.ts_energy, self.activation_energy = ts_node(self.string_energies[:len(self.truncated_string)])
        self.ts_node_num = self.string_energies.index(self.ts_energy)
        self.ts_rdkit_mol = self.string_nodes[self.ts_node_num]
        self.pdt_rdkit_mol = self.string_nodes[-1]
        
        # compute properties of the transition state
        self.forming_bond_torsion = rdMolTransforms.GetDihedralDeg(
            self.ts_rdkit_mol.GetConformer(),
            *system.reductive_elim_torsion
        )
        self.pro_dis_torsion = rdMolTransforms.GetDihedralDeg(
            self.ts_rdkit_mol.GetConformer(),
            *system.pro_dis_torsion
        )
        #compute properties of the product 
        self.formed_bond_torsion = rdMolTransforms.GetDihedralDeg(
            self.pdt_rdkit_mol.GetConformer(),
            *system.reductive_elim_torsion
        )

        self.pro_dis = 'proximal' if -90 <= self.pro_dis_torsion <= 90 else 'distal'
        # ts is exo if the torsion of the bond being formed is positive and the ts is proximal
        # if distal, the relationship is reversed
        self.endo_exo = 'exo' if (self.forming_bond_torsion >= 0) ^ (self.pro_dis == 'distal') else 'endo'
        self.syn_anti = 'syn' if -90 <= self.forming_bond_torsion <= 90 else 'anti'
        self.pdt_stereo = 'R' if self.formed_bond_torsion <= 0 else 'S'


class ConformationalSamplingDashboard(param.Parameterized):

    refresh = param.Action(lambda x: x.param.trigger('refresh'), label='Refresh')

    def __init__(self):
        super().__init__()
        self.setup_mols()
        self.dataframe()
        self.stream = Selection1D()
        self.stream_string = Selection1D()
        self.attribute_error = 0
    
    @param.depends('refresh', watch=True)
    def setup_mols(self):
        # extract the conformers for a molecule from an xyz file
        
        # self.mol_path = Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l8_degsm')
        self.mol_path = Path('/export/zimmerman/soumikd/py-conformational-sampling/example_l1_symm_xtb')
        xtb_string_paths = tuple(self.mol_path.glob('scratch/pystring_*/opt_converged_001.xyz'))
        dft_string_paths = tuple(self.mol_path.glob('scratch/pystring_*/opt_converged_002.xyz'))
        xtb_mol_confs = self.get_mol_confs(xtb_string_paths)
        dft_mol_confs = self.get_mol_confs(dft_string_paths)
        
        self.mols = { # molecule name -> conformer index -> Conformer
            f'{ligand_name}_xtb': xtb_mol_confs,
            f'{ligand_name}_dft': dft_mol_confs
        }
        
        # try out getting distances between xtb and dft conformers
        conformer_rows = []
        for conf_idx, dft_conf in dft_mol_confs.items():
            xtb_conf = xtb_mol_confs[conf_idx]
            id_list = [[i,i] for i in range(dft_conf.ts_rdkit_mol.GetNumAtoms())]
            rms = AlignMol(
                dft_conf.ts_rdkit_mol,
                xtb_conf.ts_rdkit_mol,
                atomMap=id_list
            )
            print(rms)
            
            conformer_rows.append({
                'conf_idx': conf_idx,
                # 'dft activation energy (kcal/mol)': dft_conf.activation_energy,
                # 'xtb activation energy (kcal/mol)': xtb_conf.activation_energy,
                'absolute_reactant_xtb_energy (kcal/mol)': xtb_conf.string_energies[0],
                'absolute_ts_xtb_energy (kcal/mol)': xtb_conf.ts_energy,
                'absolute_reactant_dft_energy (kcal/mol)': dft_conf.string_energies[0],
                'absolute_ts_dft_energy (kcal/mol)': dft_conf.ts_energy,
                'rmsd': rms,
            })

        df = pd.DataFrame(conformer_rows)
        df['relative_ts_xtb_energy (kcal/mol)'] = (df['absolute_ts_xtb_energy (kcal/mol)']
                                                    - df['absolute_reactant_xtb_energy (kcal/mol)'].min())
        df['relative_ts_dft_energy (kcal/mol)'] = (df['absolute_ts_dft_energy (kcal/mol)']
                                                    - df['absolute_reactant_dft_energy (kcal/mol)'].min())
        df['relative_ts_energy_diff (kcal/mol)'] = (df['relative_ts_dft_energy (kcal/mol)']
                                                       - df['relative_ts_xtb_energy (kcal/mol)'])
      
        plot = df.hvplot.scatter(y='rmsd', x='relative_ts_energy_diff (kcal/mol)', hover_cols='all')
        # plot.opts(
        #     opts.Scatter(tools=['tap', 'hover'], active_tools=['wheel_zoom'],
        #                 # width=600, height=600,
        #                 marker='triangle', size=10, fontsize={'labels': 14}),
        # )
        self.rmsd_plot = plot

        print(df)
           

    
    def get_mol_confs(self, string_paths):
        return {
            # get the conformer index for this string
            int(re.search(r"pystring_(\d+)", str(string_path)).groups()[0]):
            mol_conf
            for string_path in string_paths
            # remove any reactant structures that have already optimized to the product
            if (mol_conf := Conformer(string_path)).truncated_string # ignore conf if reactant optimized to product
        }
        
        
    @param.depends('setup_mols', watch=True)
    def dataframe(self):
        conformer_rows = []
        for mol_name, mol_confs in self.mols.items():
            for conf_idx, conformer in mol_confs.items():
                conformer_rows.append({
                    'mol_name': mol_name,
                    'conf_idx': conf_idx,
                    'activation energy (kcal/mol)': conformer.activation_energy,
                    'absolute_reactant_energy (kcal/mol)': conformer.string_energies[0],
                    'absolute_ts_energy (kcal/mol)': conformer.ts_energy,
                    'forming_bond_torsion (deg)': conformer.forming_bond_torsion,
                    'formed_bond_torsion (deg)': conformer.formed_bond_torsion,
                    'pro_dis_torsion': conformer.pro_dis_torsion,
                    'pro_dis': conformer.pro_dis,
                    'exo_endo': conformer.endo_exo,
                    'syn_anti': conformer.syn_anti,
                    'pdt_stereo': conformer.pdt_stereo,
                })
        self.df = pd.DataFrame(conformer_rows)
        self.df['relative_ts_energy (kcal/mol)'] = (self.df['absolute_ts_energy (kcal/mol)']
                                                    - self.df['absolute_reactant_energy (kcal/mol)'].min())
        
        # self.df.sort_values(axis=0, by='relative_ts_energy (kcal/mol)')['conf_idx'].to_csv(Path.home() / 'confs_df.csv')
        # self.to_csv(Path.home() / 'df.csv')
        return pn.widgets.Tabulator(self.df)
    
    
    @param.depends('stream.index', watch=True)
    def current_conformer(self):
        index = self.stream.index
        if not index:
            index = [0]
            return None
        index = index[0]
        mol_name = self.df.iloc[index]['mol_name']
        conf_index = int(self.df.iloc[index]['conf_idx'])
        return self.mols[mol_name][conf_index]
    

    @param.depends('stream_string.index', watch=True)
    def current_string_mol(self) -> Mol:
        index = self.stream_string.index
        if not index:
            index = [0]
            return None
        index = index[0]
        conformer = self.current_conformer()
        return conformer.string_nodes[index]
    
    
    @param.depends('current_conformer', watch=True)
    def conf_dataframe(self):
        conformer = self.current_conformer()
        if not conformer:
            return None
        self.conf_df = pd.DataFrame(
            {'node_num': i, 'energy (kcal/mol)': energy}
            for i, energy in enumerate(conformer.string_energies)
        )
        self.conf_df['relative_energy (kcal/mol)'] = (self.conf_df['energy (kcal/mol)']
                                                      - conformer.string_energies[0])
    
    @param.depends('dataframe')
    def scatter_plot(self):
        df = self.df
        plot = df.hvplot.box(by='mol_name', y='relative_ts_energy (kcal/mol)', c='cyan', title='Conformer Energies', height=500, width=400, legend=False) 
        plot *= df.hvplot.scatter(y='relative_ts_energy (kcal/mol)', x='mol_name', c='pdt_stereo', ylim=(0,100), hover_cols='all').opts(jitter=0.5)
        plot += self.rmsd_plot
        plot.opts(
            opts.Scatter(tools=['tap', 'hover'], active_tools=['wheel_zoom'],
                        # width=600, height=600,
                        marker='triangle', size=10, fontsize={'labels': 14}),
        )
        self.stream.update(index=[])
        self.stream.source = plot
        return plot
    
    @param.depends('conf_dataframe', watch=True)
    def string_plot(self):
        try:
            self.current_conformer()
            self.conf_dataframe()
            plot = self.conf_df.hvplot.scatter(
                y='relative_energy (kcal/mol)',
                x='node_num',
                c='blue'
            ).opts(
                opts.Scatter(tools=['tap', 'hover'], active_tools=['wheel_zoom'],
                            marker='circle', size=10, fontsize={'labels': 14}),                
            )
            self.stream_string.update(index=[])
            self.stream_string.source = plot
            return plot
        except AttributeError as error:
            self.attribute_error += 1
            return f'Attribute Error {self.attribute_error}:\n{error}'


    @param.depends('current_conformer', 'string_plot', 'current_string_mol', watch=True)
    def display_string_mol(self):
        mol = self.current_string_mol()
        if not mol:
            return None
        
        pdb_block = MolToPDBBlock(mol)
        viewer = NGLViewer(object=pdb_block, extension='pdb', background="#F7F7F7", min_height=500, sizing_mode="stretch_both")
        return viewer

    
    param.depends('display_mol', 'dataframe', 'scatter_plot', 'debug_data', 'string_plot', 'conf_dataframe', 'display_string_mol', 'stream_string.index')
    def app(self):
        return pn.Column(
            self.param.refresh,
            pn.Row(self.scatter_plot, self.display_mol),
            pn.Row(self.string_plot, self.display_string_mol),
            self.debug_data,
            self.dataframe
        )
    
    
    @param.depends('current_conformer', 'scatter_plot', 'stream_string.index', watch=True)
    def display_mol(self):
        conformer = self.current_conformer()
        if not conformer:
            return None
        pdb_block = MolToPDBBlock(conformer.ts_rdkit_mol)
        viewer = NGLViewer(object=pdb_block, extension='pdb', background="#F7F7F7", min_height=500, sizing_mode="stretch_both")
        return viewer


    @param.depends('refresh', 'stream.index', watch=True)
    def debug_data(self):
        # index = self.stream.index
        # return index
        return (f'{self.stream.index = }\n{repr(dashboard) = }\n'
                + f'{self.free_energy_R_minus_S() = }')


    def free_energy_R_minus_S(self):
        group_by = self.df.groupby('pdt_stereo')['relative_ts_energy (kcal/mol)'].apply(list)
        return free_energy_diff(group_by['S'], group_by['R'], temperature=358.15)
    

dashboard = ConformationalSamplingDashboard()
try: # reboot server if already running in interactive mode
    bokeh_server.stop()
except (NameError, AssertionError):
    pass
bokeh_server = dashboard.app().show()

# dashboard.app()

# %%
