import logging
from dataclasses import dataclass, field

import stk

from conformational_sampling.config import Config
from conformational_sampling.main import (
    ConformerEnsembleOptimizer,
    bind_ligands,
    gen_confs_openbabel,
)


@dataclass
class CatalyticReactionComplex:
    metal: stk.BuildingBlock
    ancillary_ligand: stk.BuildingBlock
    reactive_ligand_1: stk.BuildingBlock
    reactive_ligand_2: stk.BuildingBlock
    config: Config = field(default_factory=Config)

    def __post_init__(self) -> None:
        """Create a complex with one ancillary ligand and two reactive ligands"""
        self.complex = bind_ligands(
            self.metal,
            self.ancillary_ligand,
            self.reactive_ligand_1,
            self.reactive_ligand_2,
        )

    def gen_conformers(self):
        reactive_ligand_1_conformers = gen_confs_openbabel(
            self.reactive_ligand_1, self.config
        )
        reactive_ligand_2_conformers = gen_confs_openbabel(
            self.reactive_ligand_2, self.config
        )
        ancillary_ligand_conformers = gen_confs_openbabel(
            self.ancillary_ligand, self.config
        )
        self.unoptimized_conformers = [
            bind_ligands(self.metal, ancillary, ligand_1, ligand_2)
            for ancillary in ancillary_ligand_conformers
            for ligand_1 in reactive_ligand_1_conformers
            for ligand_2 in reactive_ligand_2_conformers
        ]
        self.optimized_conformers = ConformerEnsembleOptimizer(
            self.unoptimized_conformers, self.config
        ).optimize()
        logging.debug('Finished generating CatalyticReactionComplex conformers')
        return self.optimized_conformers

    def _breaking_bonds(self):
        breaking_bonds = []

        # Identify metal / reactive ligand bonds
        for bond_info in self.complex.get_bond_infos():
            if bond_info.get_building_block() is None:
                bond = bond_info.get_bond()
                atom_ids = {bond.get_atom1().get_id(), bond.get_atom2().get_id()}
                atom_infos = set(self.complex.get_atom_infos(atom_ids))
                building_blocks = {
                    atom_info.get_building_block() for atom_info in atom_infos
                }
                if self.metal in building_blocks:
                    building_blocks.remove(self.metal)
                    other_building_block = building_blocks.pop()
                    if other_building_block in {
                        self.reactive_ligand_1,
                        self.reactive_ligand_2,
                    }:
                        # This bond is broken in reductive elimination
                        breaking_bonds.append(atom_ids)

        # There should be one bond from the metal to each reactive ligand
        assert len(breaking_bonds) == 2
        return breaking_bonds

    def _forming_bond(self):
        breaking_bonds = self._breaking_bonds()
        return breaking_bonds[0] ^ breaking_bonds[1]

    def gen_reductive_elim_drive_coords(self):
        """Generate reductive elimination driving coordinates for this complex"""
        breaking_bonds = self._breaking_bonds()

        # Form a bond between the distinct atom from each metal-ligand bond
        forming_bond = self._forming_bond()

        driving_coordinates = [
            ('BREAK', *breaking_bond) for breaking_bond in breaking_bonds
        ]
        driving_coordinates.append(('ADD', *forming_bond))

        # Convert from zero-indexed to one-indexed for pyGSM
        return [(type, i + 1, j + 1) for (type, i, j) in driving_coordinates]

    def forming_bond_torsion(self):
        forming_bond = list(self._forming_bond())

        def connected_atom_within_building_block(atom_id):
            for bond_info in self.complex.get_bond_infos():
                if bond_info.get_building_block() is not None:
                    bond = bond_info.get_bond()
                    bond_atom_ids = {
                        bond.get_atom1().get_id(),
                        bond.get_atom2().get_id(),
                    }
                    if atom_id in bond_atom_ids:
                        bond_atom_ids.remove(atom_id)
                        return bond_atom_ids.pop()

        return (
            connected_atom_within_building_block(forming_bond[0]),
            forming_bond[0],
            forming_bond[1],
            connected_atom_within_building_block(forming_bond[1]),
        )

        # for torsion in (
        #     stko.TorsionCalculator().get_results(self.complex).get_torsions()
        # ):
        #     torsion_bond = set(list(torsion.get_atom_ids())[1:3])
        #     if torsion_bond == forming_bond:
        #         return torsion
        assert False  # didn't find forming bond torsion
