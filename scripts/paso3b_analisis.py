import numpy as np
from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate
from Bio.PDB import PDBParser, NeighborSearch
import re

def esqueleto(m):
    m = Chem.RemoveHs(m, sanitize=False)
    rw = Chem.RWMol(m)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0); a.SetIsAromatic(False); a.SetNumExplicitHs(0); a.SetNoImplicit(True)
    return rw.GetMol(), m

ref = Chem.SDMolSupplier("FA_cristal_ref.sdf", removeHs=False)[0]
refS, refH = esqueleto(ref)
refxyz = refH.GetConformer().GetPositions()

core = Chem.MolFromSmarts("[#7]~[#6]1~[#7]~[#6]2~[#7]~[#6]~[#6]~[#7]~[#6]2~[#6](~[#8])~[#7]~1")
core_ref = sorted(refS.GetSubstructMatch(core))
print("Átomos de la pterina en el ligando cristal:", len(core_ref))

pm = PDBQTMolecule.from_file("redock_out.pdbqt", skip_typing=True)
probe = RDKitMolCreate.from_pdbqt_mol(pm)[0]
probeS, probeH = esqueleto(probe)
matches = refS.GetSubstructMatches(probeS, uniquify=False, useChirality=False, maxMatches=100000)
scores = [float(x) for x in re.findall(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", open("redock_out.pdbqt").read())]

prot = PDBParser(QUIET=True).get_structure("p", "proteina.pdb")
ns = NeighborSearch(list(prot.get_atoms()))
def contactos(xyz, d=4.0):
    s = set()
    for p in xyz:
        for a in ns.search(p, d):
            r = a.get_parent(); s.add((r.get_resname(), r.id[1]))
    return s

cref = contactos(refxyz[core_ref])
print("Residuos a <4 A de la pterina en el cristal:")
print("  ", ", ".join(f"{n}{i}" for n, i in sorted(cref, key=lambda x: x[1])))

print("\nPose  Score  RMSD_total  RMSD_pterina  Contactos_recuperados")
for i in range(probeH.GetNumConformers()):
    pxyz = probeH.GetConformer(i).GetPositions()
    mejor = None
    for m in matches:
        tot = np.sqrt(np.mean(np.sum((refxyz[list(m)] - pxyz) ** 2, axis=1)))
        sel = [k for k, x in enumerate(m) if x in core_ref]
        pt = np.sqrt(np.mean(np.sum((refxyz[[m[k] for k in sel]] - pxyz[sel]) ** 2, axis=1)))
        if mejor is None or pt < mejor[1]:
            mejor = (tot, pt, sel)
    cp = contactos(pxyz[mejor[2]])
    frac = len(cp & cref) / len(cref)
    print(f"{i+1:4d}  {scores[i]:6.2f}  {mejor[0]:9.2f}  {mejor[1]:11.2f}  {frac:12.0%}")
