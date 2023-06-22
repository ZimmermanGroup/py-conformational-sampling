# %%
%reload_ext autoreload
%autoreload 2
from dataclasses import dataclass
import re
from IPython.display import display
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolToPDBBlock

import stk
# from conformer_rl.analysis.lignin_contacts import CONF_ID, FUNC_GROUP_ID_1, setup_mol
# from conformer_rl.analysis.lignin_pericyclic import \
#     LigninPericyclicCalculator, LigninPericyclicFunctionalGroupFactory, \
#     LigninMaccollCalculator, LigninMaccollFunctionalGroupFactory, init_stk_from_rdkit

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
pn.extension("ngl_viewer", sizing_mode="stretch_width")
from pathlib import Path
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock
import openbabel as ob
from openbabel import pybel as pb
import nglview

from pygsm.utilities.units import KCAL_MOL_PER_AU
from conformational_sampling.analyze import ts_node

    
# setup_mol()
@dataclass
class Conformer:
    string_path: Path
    
    def __post_init__(self):
        string_nodes = list(pb.readfile('xyz', str(self.string_path)))
        self.string_nodes = [MolFromMolBlock(node.write('mol'), removeHs=False)
                             for node in string_nodes] # convert to rdkit
        self.string_energies = [float(MolToMolBlock(node).split()[0]) * KCAL_MOL_PER_AU
                                for node in self.string_nodes]
        max_diff, self.ts_energy, activation_barrier = ts_node(self.string_energies)
        self.ts_node_num = self.string_energies.index(self.ts_energy)
        self.ts_rdkit_mol = self.string_nodes[self.ts_node_num]

class ConformationalSamplingDashboard(param.Parameterized):

    refresh = param.Action(lambda x: x.param.trigger('refresh'), label='Refresh')

    def __init__(self):
        super().__init__()
        self.setup_mols()
        self.dataframe()
        self.stream = Selection1D()
    
    @param.depends('refresh', watch=True)
    def setup_mols(self):
        # extract the conformers for a molecule from an xyz file
        
        mol_path = Path('/export/zimmerman/soumikd/py-conformational-sampling/')
        string_paths = tuple(mol_path.glob('scratch/pystring_*/opt_converged_000.xyz'))
        self.mol_confs = {
            # get the conformer index for this string
            int(re.search(r"pystring_(\d+)", str(string_path)).groups()[0]):
            Conformer(string_path)
            for string_path in string_paths
        }
        
        self.mols = {'ligand_l8': self.mol_confs} # molecule name -> conformer index -> Conformer
        
    @param.depends('setup_mols', watch=True)
    def dataframe(self):
        conformer_rows = []
        for mol_name, mol_confs in self.mols.items():
            for conf_idx, conformer in mol_confs.items():
                conformer_rows.append({
                    'mol_name': mol_name,
                    'conf_idx': conf_idx,
                    'ts energy (kcal/mol)': conformer.ts_energy,
                })
        self.df = pd.DataFrame(conformer_rows)
    
    @param.depends('dataframe')
    def scatter_plot(self):
        df = self.df
        plot = df.hvplot.box(by='mol_name', y='ts energy (kcal/mol)', c='orange', title='Conformer Energies', height=800, width=400, legend=False) 
        plot *= df.hvplot.scatter(y='ts energy (kcal/mol)', x='mol_name', c='blue').opts(jitter=0.5)
        plot.opts(
            opts.Scatter(tools=['tap', 'hover'], active_tools=['wheel_zoom'],
                        # width=600, height=600,
                        marker='triangle', size=10, fontsize={'labels': 14}),
        )
        self.stream.update(index=[])
        self.stream.source = plot
        return plot
    
    
    param.depends('display_mol', 'dataframe', 'scatter_plot', 'index_conf')
    def app(self):
        return pn.Row(pn.Column(self.param.refresh, self.scatter_plot, self.index_conf), self.display_mol)
    
    
    @param.depends('stream.index', 'scatter_plot', watch=True)
    def display_mol(self):
        index = self.stream.index
        if not index: # abort if nothing is selected
            index = [0]
            # return None
        index = index[0]
        mol_name = self.df.iloc[index]['mol_name']
        conf_index = int(self.df.iloc[index]['conf_idx'])
        pdb_block = MolToPDBBlock(self.mols[mol_name][conf_index].ts_rdkit_mol)
        viewer = NGLViewer(object=pdb_block, extension='pdb', background="#F7F7F7", min_height=800, sizing_mode="stretch_both")
        return viewer


    @param.depends('stream.index', watch=True)
    def index_conf(self):
        # index = self.stream.index
        # return index
        return f'{self.stream.index = }\n{repr(dashboard) = }'
    

dashboard = ConformationalSamplingDashboard()
try: # reboot server if already running in interactive mode
    bokeh_server.stop()
except (NameError, AssertionError):
    pass
bokeh_server = dashboard.app().show(port=45350)
# dashboard.app()
