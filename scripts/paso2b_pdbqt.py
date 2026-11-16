from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

mol = Chem.SDMolSupplier("FA_dock_input.sdf", removeHs=False)[0]
setups = MoleculePreparation().prepare(mol)
pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
print("ok:", ok, err)
open("FA_dock_input.pdbqt", "w").write(pdbqt)
print("torsiones:", pdbqt.count("BRANCH") // 2)