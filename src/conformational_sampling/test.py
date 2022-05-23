# %%

import openbabel as ob
from openbabel import pybel as pb
test = list(pb.readfile('xyz', '/home/paulzim/zstruct2/test/formwater_mopac/stringfile.xyz11204'))
for mol in pb.readfile('xyz', '/home/paulzim/zstruct2/test/formwater_mopac/stringfile.xyz11204'):
    print(len(mol.atoms))
# mol = ob.OBMol()
# print(mol.NumAtoms())