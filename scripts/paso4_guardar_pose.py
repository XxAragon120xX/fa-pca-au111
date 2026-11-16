from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate

pm = PDBQTMolecule.from_file("redock_caja18.pdbqt", skip_typing=True)
mol = RDKitMolCreate.from_pdbqt_mol(pm)[0]
Chem.MolToMolFile(mol, "FA_pose1_caja18.sdf", confId=0)
print("Guardada la pose 1:", mol.GetNumAtoms(), "atomos")
