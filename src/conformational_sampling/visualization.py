# %%
# %reload_ext autoreload
# %autoreload 2
import os
import re
import pickle
# from IPython.display import display
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.rdmolfiles import MolToPDBBlock, MolToXYZBlock
from rdkit.Chem.rdMolAlign import AlignMol

print('Made it halfway through imports!')

import param
import hvplot.pandas # noqa
import holoviews as hv
from holoviews import opts, Slope
from holoviews.streams import Selection1D
import panel as pn
from panel_chemistry.pane import \
    NGLViewer  # panel_chemistry needs to be imported before you run pn.extension()
from panel_chemistry.pane.ngl_viewer import EXTENSIONS
pn.extension('bokeh')
hv.extension('bokeh')
pn.extension(comms='vscode')
pn.extension('tabulator')
pn.extension("ngl_viewer", sizing_mode="stretch_width")
from pathlib import Path
import openbabel as ob

from conformational_sampling.analyze import Conformer, systems, exclude_confs
from conformational_sampling.utils import free_energy_diff

print('Finished imports!')

READ_FROM_PICKLE = True
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
        
        pickle_path = Path.home() / 'mols.pkl'
        if pickle_path.exists() and READ_FROM_PICKLE:
            with pickle_path.open('rb') as file:
                self.mols = pickle.load(file)
        else:
            self.mols = {}
            for ligand_name, system in systems.items():
                # string_paths = tuple(system.mol_path.glob('scratch/pystring_*/opt_converged_001.xyz'))
                string_paths = tuple(
                    system.mol_path.glob("scratch/pystring_*/xtb_singlepoints_001.xyz")
                )
                self.mols[ligand_name] = self.get_mol_confs(system, string_paths)
            with pickle_path.open("wb") as file:
                pickle.dump(self.mols, file)
        
    
    def get_mol_confs(self, system, string_paths):
        return {
            # get the conformer index for this string
            int(re.search(r"pystring_(\d+)", str(string_path)).groups()[0]):
            mol_conf
            for string_path in string_paths
            # remove any reactant structures that have already optimized to the product
            if ((mol_conf := Conformer(system, string_path, singlepoints=True)).truncated_string
                and mol_conf.activation_energy <= 200)
        }
        
        
    @param.depends('setup_mols', watch=True)
    def dataframe(self):
        conformer_rows = []
        for mol_name, mol_confs in self.mols.items():
            for conf_idx, conformer in mol_confs.items():
                # remove conformers manually identified as problematic
                if mol_name in exclude_confs and conf_idx in exclude_confs[mol_name]:
                    continue
                        
                if ((atrop_torsion := conformer.system.atrop_torsion)
                    and rdMolTransforms.GetDihedralDeg(conformer.string_nodes[0].GetConformer(), *atrop_torsion) < 0):
                        continue

                # wrap forming_bond_torsion differently based on stereochemistry
                forming_bond_torsion = conformer.forming_bond_torsion
                if conformer.forming_bond_torsion < -90 and conformer.pdt_stereo == 'S':
                    forming_bond_torsion += 360
                elif conformer.forming_bond_torsion > 90 and conformer.pdt_stereo == 'R':
                    forming_bond_torsion -= 360
                    
                conformer_rows.append({
                    'mol_name': mol_name,
                    'conf_idx': conf_idx,
                    'activation energy (kcal/mol)': conformer.activation_energy,
                    'absolute_reactant_energy (kcal/mol)': conformer.min_reac_energy,
                    'absolute_ts_energy (kcal/mol)': conformer.ts_energy,
                    'forming_bond_torsion (deg)': forming_bond_torsion,
                    'formed_bond_torsion (deg)': conformer.formed_bond_torsion,
                    'pro_dis_torsion': conformer.pro_dis_torsion,
                    'improper_torsion': conformer.improper_torsion,
                    'improper_torsion_ts': conformer.improper_torsion_ts,
                    'ligands_angle': conformer.ligands_angle,
                    'out_of_plane_angle': conformer.out_of_plane_angle,
                    'pro_dis': conformer.pro_dis,
                    'exo_endo': conformer.endo_exo,
                    'syn_anti': conformer.syn_anti,
                    'Product stereochemistry': conformer.pdt_stereo,
                    'tau_4_prime': conformer.tau_4_prime,
                    'tau_4_prime_ts':conformer.tau_4_prime_ts,
                })
        self.df = pd.DataFrame(conformer_rows)
        self.df['relative_ts_energy (kcal/mol)'] = (
            self.df['absolute_ts_energy (kcal/mol)']
            - self.df.groupby('mol_name')['absolute_reactant_energy (kcal/mol)'].transform('min')
        )
        
        # self.df.sort_values(axis=0, by='relative_ts_energy (kcal/mol)')['conf_idx'].to_csv(Path.home() / 'confs_df.csv')
        self.df.to_csv(Path.home() / 'df.csv')
        # return pn.widgets.Tabulator(self.df)
        
        plot = self.df.hvplot(
            by=['Product stereochemistry'],
            fontscale=1.2,
            groupby=['mol_name'],
            kind='scatter',
            x='forming_bond_torsion (deg)',
            xlabel='Reductive elimination torsion at TS (deg)',
            y=['relative_ts_energy (kcal/mol)'],
            ylabel='Relative TS energy (kcal/mol)',
            title='',
            width=600,
            height=500,
            xlim=(-220,220),
            ylim=(0,40),
        ).opts(legend_position='top_left')
        # plot *= hv.VLine(-180).opts(color='black', line_dash='dashed', line_width=3)
        # plot *= hv.VLine(180).opts(color='black', line_dash='dashed', line_width=3)
        # plot.opts(opts.Scatter())
        # hvplot.save(plot, Path.home() / 'Lilly' / 'py-conformational-sampling' / 'test_plot.svg', fmt='svg')
        
        # testing out generating descriptive statistics
        pd.options.display.width = 200
        # self.df.groupby('mol_name').describe().to_csv(Path.home() / 'Lilly' / 'py-conformational-sampling' / 'summary_statistics.csv')
        sum_stats_file = Path.home() / 'Lilly' / 'py-conformational-sampling' / 'summary_statistics.csv'
        group_by = self.df.groupby(['mol_name', 'Product stereochemistry'])
        sum_stats = group_by[[
            'activation energy (kcal/mol)',
            'relative_ts_energy (kcal/mol)',
            'forming_bond_torsion (deg)',
        ]].describe(percentiles=[0.5]).T.to_csv(sum_stats_file, float_format=lambda x: f' {x:.6f}')
        # sum_stats_file.write_text(sum_stats)

        return self.df.hvplot.explorer(
            by=['Product stereochemistry'],
            groupby=['mol_name'],
            kind='scatter',
            x='forming_bond_torsion (deg)',
            y=['relative_ts_energy (kcal/mol)'],
            width=700,
            ylim=(0,40),
        )
    
    
    @param.depends('stream.index', watch=True)
    def current_conformers(self):
        index = self.stream.index
        if not index:
            index = [0]
            return None
        # index = index[0]
        mols = []
        for idx in index:
            mol_name = self.df.iloc[idx]['mol_name']
            conf_index = int(self.df.iloc[idx]['conf_idx'])
            mols.append(self.mols[mol_name][conf_index])
        return mols
    

    @param.depends('stream_string.index', watch=True)
    def current_string_mol(self) -> Mol:
        index = self.stream_string.index
        if not index:
            index = [0]
            return None
        index = index[0]
        conformer = self.current_conformers()[0]
        return conformer.string_nodes[index]
    
    
    @param.depends('current_conformers', watch=True)
    def conf_dataframe(self):
        conformers = self.current_conformers()
        if not conformers:
            return None
        conformer = conformers[0] # use the first conformer selected
        self.conf_df = pd.DataFrame(
            {'node_num': i, 'energy (kcal/mol)': energy}
            for i, energy in enumerate(conformer.string_energies)
        )
        self.conf_df['relative_energy (kcal/mol)'] = (self.conf_df['energy (kcal/mol)']
                                                      - conformer.string_energies[0])
    
    @param.depends('dataframe')
    def scatter_plot(self):
        df = self.df
        plot = df.hvplot.scatter(
            y='relative_ts_energy (kcal/mol)', x='forming_bond_torsion (deg)',
            c='mol_name', #groupby='mol_name',
            ylim=(0,100), hover_cols='all'
        )
        # plot += self.rmsd_plot
        plot.opts(
            opts.Scatter(tools=['tap', 'hover'], active_tools=['wheel_zoom'],
                        # width=600, height=600,
                        marker='triangle', size=10, fontsize={'labels': 14}),
        )
        self.stream.update(index=[])
        self.stream.source = plot
        return plot
    
        df = self.df
        plot = df.hvplot.box(
            by='mol_name', y='relative_ts_energy (kcal/mol)', c='cyan',
            title='Conformer Energies', height=500, width=400, legend=False
        ) 
        plot *= df.hvplot.scatter(
            y='relative_ts_energy (kcal/mol)', x='mol_name', c='Product stereochemistry',
            ylim=(0,100), hover_cols='all'
        ).opts(jitter=0.5)
        # plot += self.rmsd_plot
        plot.opts(
            opts.Scatter(tools=['tap', 'hover'], active_tools=['wheel_zoom'],
                        # width=600, height=600,
                        marker='triangle', size=10, fontsize={'labels': 14}),
        )
        self.stream.update(index=[])
        self.stream.source = plot
        return plot
    
    @param.depends('conf_dataframe', 'stream.index', watch=True)
    def string_plot(self):
        try:
            self.current_conformers()
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


    @param.depends('current_conformers', 'string_plot', 'current_string_mol', watch=True)
    def display_string_mol(self):
        mol = self.current_string_mol()
        if not mol:
            return None
        
        pdb_block = MolToPDBBlock(mol)
        viewer = NGLViewer(object=pdb_block, extension='pdb', background="#F7F7F7", min_height=400, sizing_mode="stretch_both")
        print(MolToXYZBlock(mol))
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
    
    
    @param.depends('current_conformers', 'scatter_plot', 'stream_string.index', watch=True)
    def display_mol(self):
        conformers = self.current_conformers()
        if not conformers:
            return None
        viewers = []
        for conformer in conformers:
            pdb_block = MolToPDBBlock(conformer.ts_rdkit_mol)
            viewer = NGLViewer(
                object=pdb_block,
                extension='pdb',
                background="#F7F7F7",
                min_height=400,
                sizing_mode="stretch_both"
            )
            viewers.append(viewer)
        return pn.Row(*viewers)


    @param.depends('refresh', 'stream.index', watch=True)
    def debug_data(self):
        # index = self.stream.index
        # return index
        return (f'{self.stream.index = }\n{repr(dashboard) = }\n'
                + f'{self.free_energy_R_minus_S() = }')


    def free_energy_R_minus_S(self):
        group_by = self.df.groupby(['mol_name', 'Product stereochemistry'])['relative_ts_energy (kcal/mol)'].apply(list)
        free_energy_diffs = {}
        for mol_name in self.df['mol_name'].unique():
            free_energy_diffs[mol_name] = free_energy_diff(
                group_by[(mol_name, 'S')],
                group_by[(mol_name, 'R')],
                temperature=358.15
            )
        return free_energy_diffs
    

dashboard = ConformationalSamplingDashboard()
dashboard.app().servable()  # for running dashboard with panel serve

# below code block for running panel in interactive vscode window
# try: # reboot server if already running in interactive mode
#     bokeh_server.stop()
# except (NameError, AssertionError):
#     pass
# bokeh_server = dashboard.app().show()

# %%
# test_df = pd.read_csv(Path.home() / 'df.csv')
# groupby = test_df.value_counts(['mol_name', 'exo_endo', 'syn_anti', 'Product stereochemistry']).reset_index(name='count')
# display(groupby)

# pn.panel(hvplot.explorer(groupby)).show()

# pass
# %%
# 'fontsize': {'labels': 10}, 'title': '', 'group_label': 'Prod'
