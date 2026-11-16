import subprocess, re
import numpy as np
from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate

CENTRO = (6.43, 21.76, 2.22)
CAJA = 25

cmd = [r".\vina_1.2.7_win.exe",
       "--receptor", "proteina.pdbqt", "--ligand", "FA_dock_input.pdbqt",
       "--center_x", str(CENTRO[0]), "--center_y", str(CENTRO[1]), "--center_z", str(CENTRO[2]),
       "--size_x", str(CAJA), "--size_y", str(CAJA), "--size_z", str(CAJA),
       "--exhaustiveness", "32", "--num_modes", "20", "--energy_range", "4",
       "--seed", "42", "--out", "redock_out.pdbqt"]
print("Corriendo Vina...")
subprocess.run(cmd, check=True)

scores = [float(x) for x in re.findall(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", open("redock_out.pdbqt").read())]

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

pm = PDBQTMolecule.from_file("redock_out.pdbqt", skip_typing=True)
probe = RDKitMolCreate.from_pdbqt_mol(pm)[0]
probeS, probeH = esqueleto(probe)
matches = refS.GetSubstructMatches(probeS, uniquify=False, useChirality=False, maxMatches=100000)
print("Átomos pesados ref/docked:", refH.GetNumAtoms(), probeH.GetNumAtoms(), "| mapeos simétricos:", len(matches))

res = []
for i in range(probeH.GetNumConformers()):
    pxyz = probeH.GetConformer(i).GetPositions()
    best = min(np.sqrt(np.mean(np.sum((refxyz[list(m)] - pxyz) ** 2, axis=1))) for m in matches)
    res.append(best)

print("\nPose  Score(kcal/mol)  RMSD vs cristal (A)")
for i, r in enumerate(res):
    print(f"{i+1:4d}  {scores[i] if i < len(scores) else float('nan'):10.2f}  {r:12.2f}")
print("\nPose 1:", round(res[0], 2), "A | mejor RMSD:", round(min(res), 2), "A | poses < 2 A:", sum(r < 2 for r in res))
