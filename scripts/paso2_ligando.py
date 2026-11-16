from rdkit import Chem
from rdkit.Chem import AllChem

# FA neutro (plantilla para los órdenes de enlace del cristal)
SMI_NEUTRO = "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CCC(=O)O)C(=O)O)nc2c(=O)[nH]1"
# FA a pH 7.4: ambos carboxilatos desprotonados (carga -2)
SMI_PH74   = "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CCC(=O)[O-])C(=O)[O-])nc2c(=O)[nH]1"

# --- Referencia cristalográfica (átomos pesados, con órdenes de enlace) ---
pdb  = Chem.MolFromPDBFile("ligando_cristal.pdb", removeHs=True, proximityBonding=True)
tmpl = Chem.MolFromSmiles(SMI_NEUTRO)
print("Átomos pesados cristal / plantilla:", pdb.GetNumHeavyAtoms(), tmpl.GetNumHeavyAtoms())

ref = AllChem.AssignBondOrdersFromTemplate(tmpl, pdb)
for m in ref.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[OD1]")):
    a = ref.GetAtomWithIdx(m[2])
    a.SetFormalCharge(-1); a.SetNumExplicitHs(0); a.SetNoImplicit(True)
Chem.SanitizeMol(ref)
Chem.AssignStereochemistryFrom3D(ref)
print("Carga del ligando cristal:", Chem.GetFormalCharge(ref))
print("Estereocentros (esperado S):", Chem.FindMolChiralCenters(ref, useLegacyImplementation=False))
Chem.MolToMolFile(ref, "FA_cristal_ref.sdf")

# --- Entrada para docking: conformero generado desde cero (no la pose del cristal) ---
mol = Chem.AddHs(Chem.MolFromSmiles(SMI_PH74))
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
mol.SetProp("_Name", "FA_pH74")
Chem.MolToMolFile(mol, "FA_dock_input.sdf")
print("Listo: FA_cristal_ref.sdf y FA_dock_input.sdf")