from dataclasses import dataclass

import stk

from conformational_sampling.main import bind_ligands


@dataclass
class CatalyticReactionComplex:
    metal: stk.Molecule
    ancillary_ligand: stk.Molecule
    reactive_ligand_1: stk.Molecule
    reactive_ligand_2: stk.Molecule

    def __post_init__(self) -> None:
        "create a complex with ane ancillary ligand and two reactive ligands, like bind_ligands"
        self.complex = bind_ligands(
            self.metal,
            self.ancillary_ligand,
            self.reactive_ligand_1,
            self.reactive_ligand_2,
        )

    def gen_reductive_elim_drive_coords(self):
        "generate reductive elimination driving coordinates for this complex"
        atom_infos = list(self.complex.get_atom_infos())
        driving_coordinates = []

        def is_metal_reactive_ligand_bond(
            atom_info_1: stk.AtomInfo, atom_info_2: stk.AtomInfo
        ) -> bool:
            building_block_1 = atom_info_1.get_building_block()
            building_block_2 = atom_info_2.get_building_block()
            return building_block_1 is self.metal and building_block_2 in (
                self.reactive_ligand_1,
                self.reactive_ligand_2,
            )

        # determine driving coordinates to break
        for bond_info in self.complex.get_bond_infos():
            if bond_info.get_building_block() is None:
                bond = bond_info.get_bond()
                atom_id_1, atom_id_2 = (
                    bond.get_atom1().get_id(),
                    bond.get_atom2().get_id(),
                )
                atom_info_1 = next(self.complex.get_atom_infos(atom_id_1))
                atom_info_2 = next(self.complex.get_atom_infos(atom_id_2))
                if is_metal_reactive_ligand_bond(
                    atom_info_1, atom_info_2
                ) or is_metal_reactive_ligand_bond(atom_info_2, atom_info_1):
                    # this bond is broken in reductive elimination
                    driving_coordinates += ('BREAK', atom_id_1, atom_id_2)

        # CHECK THAT THERE ARE TWO BREAKING DRIVING COORDS

        # go from zero-indexed to one-indexed
        return [(type, i + 1, j + 1) for (type, i, j) in driving_coordinates]