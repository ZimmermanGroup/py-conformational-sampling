# %%
%reload_ext autoreload
%autoreload 2
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

def setup_mols():
    # extract the conformers for a molecule from an xyz file
    def setup_mol(path):
        mol_confs = list(pb.readfile('xyz', str(path)))
        mol_confs = [MolFromMolBlock(conf.write('mol'), removeHs=False) for conf in mol_confs]
        return mol_confs
    
    # paths = tuple(Path('/export/zimmerman/joshkamm/Lilly/ConformationalSampling-1/examples/').glob('*/all_openbabel5.xyz'))
    paths = tuple(Path('/export/zimmerman/joshkamm/Lilly/py-conformational-sampling/examples/').glob('*/conformers_3_xtb.xyz'))
    mols = {path.parts[-2] : setup_mol(path) for path in paths}
    return mols
    
# setup_mol()

class ConformationalSamplingDashboard(param.Parameterized):
    # mechanism = param.Selector(['Pericylic', 'Maccoll'])
    # mechanism = param.ObjectSelector(default=LigninPericyclicCalculator,
    #                                  objects=[LigninPericyclicCalculator, LigninMaccollCalculator])
    # func_group_factories = {LigninPericyclicCalculator: LigninPericyclicFunctionalGroupFactory,
    #                         LigninMaccollCalculator: LigninMaccollFunctionalGroupFactory}

    def __init__(self):
        super().__init__()
        self.mols = setup_mols()
        self.df = self.dataframe()
        self.stream = Selection1D()
    
    # @param.depends('mechanism', watch=True)
    def dataframe(self):
        def mol_dataframe(mol_confs):
            energies = [MolToMolBlock(conf).split()[0] for conf in mol_confs]
            energies = pd.Series(energies, name='energies (kcal/mol)')
            energies = pd.to_numeric(energies, errors='coerce')
            energies -= energies.min()
            energies = energies.where(energies <= 100) * 627.5 #hartree -> kcal/mol
            # indices = pd.Series(range(len(energies), name='conf_index'))
            return energies
        df = pd.concat({name : mol_dataframe(mol_confs) for (name, mol_confs) in self.mols.items()},
                       names=['mol_name', 'idx'])
        return df.reset_index()
    #     distances = self.mechanism().calculate_distances(self.mol)
    #     self.df = distances.to_dataframe().reset_index().astype({FUNC_GROUP_ID_1: 'str'})
    
    @param.depends('dataframe')
    def scatter_plot(self):
        df = self.dataframe()
        plot = df.hvplot.box(by='mol_name', y='energies (kcal/mol)', c='orange', title='Conformer Energies', height=400, width=400, legend=False) 
        plot *= df.hvplot.scatter(y='energies (kcal/mol)', x='mol_name', c='blue').opts(jitter=0.5)
        plot.opts(
            opts.Scatter(tools=['tap', 'hover'], active_tools=['wheel_zoom'],
                        # width=600, height=600,
                        marker='triangle', size=10, fontsize={'labels': 14}),
        )
        self.stream.update(index=[])
        self.stream.source = plot
        return plot
    #     if self.mechanism == LigninMaccollCalculator:
    #         points = self.df.hvplot.scatter(x='Lignin Maccoll mechanism distances', y='Energies', c='func_group_id_1')
    #     elif self.mechanism == LigninPericyclicCalculator:
    #         points = self.df.hvplot.scatter(x='Lignin retro-ene mechanism distances',
    #                                         y='Inhibition distance differences')
    #     return points
    
    # param.depends('display_mol', 'scatter_plot', 'index_conf', 'disp_mechanism', 'mechanism')
    param.depends('display_mol', 'dataframe', 'scatter_plot', 'index_conf')
    def app(self):
        # return pn.Row(pn.Column(self.param.mechanism, self.scatter_plot, self.index_conf), self.display_mol)
        return pn.Row(pn.Column(self.index_conf, self.dataframe), self.scatter_plot, self.display_mol)
    
    # @param.depends('mechanism', watch=True)
    # def setup_stk_mol(self):
    #     return init_stk_from_rdkit(
    #         self.mol,
    #         functional_groups=(self.func_group_factories[self.mechanism](),),
    #     )
        
    @param.depends('stream.index', 'scatter_plot', watch=True)
    def display_mol(self):
        index = self.stream.index
        if not index: # abort if nothing is selected
            return None
        index = index[0]
        mol_name = self.df.iloc[index]['mol_name']
        conf_index = int(self.df.iloc[index]['idx'])
        pdb_block = MolToPDBBlock(self.mols[mol_name][conf_index])
        # pdb_block = MolToPDBBlock(self.highlighted_mol(self.df.iloc[index][FUNC_GROUP_ID_1]), confId=conf_id)
        viewer = NGLViewer(object=pdb_block, extension='pdb', background="#F7F7F7", min_height=800, sizing_mode="stretch_both")
        return viewer

    @param.depends('stream.index', watch=True)
    def index_conf(self):
        # index = self.stream.index
        # return index
        return f'{self.stream.index = }\n{repr(dashboard) = }'
    
    # @param.depends('mechanism', watch=True)
    # def disp_mechanism(self):
    #     return f'{self.param.mechanism = }'

    # def highlighted_mol(self, func_group_id):
    #     mol = setup_mol()
    #     for i, f_group in enumerate(self.setup_stk_mol().get_functional_groups()):
    #         if self.mechanism == LigninMaccollCalculator:
    #             atom_1 = mol.GetAtomWithIdx(f_group.H.get_id())
    #             atom_2 = mol.GetAtomWithIdx(f_group.O.get_id())
    #             if i == int(func_group_id):
    #                 atom_1.SetAtomicNum(9)
    #                 atom_2.SetAtomicNum(15)
    #             else: # gently highlight other functional groups in the same conformer
    #                 atom_1.SetAtomicNum(2)
    #         elif self.mechanism == LigninPericyclicCalculator:
    #             atom_1 = mol.GetAtomWithIdx(f_group.H_phenyl.get_id())
    #             atom_2 = mol.GetAtomWithIdx(f_group.H_alkyl.get_id())
    #             atom_1.SetAtomicNum(10)
    #             atom_2.SetAtomicNum(10)
    #     return mol

dashboard = ConformationalSamplingDashboard()
try:
    bokeh_server.stop()
except (NameError, AssertionError):
    pass
bokeh_server = dashboard.app().show(port=45350)
# dashboard.app()
