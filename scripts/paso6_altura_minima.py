import numpy as np
from rdkit import Chem
from Bio.PDB import PDBParser

CONTACTO = 3.0   # A: distancia minima aceptable proteina - capa superior de Au
ALTURAS_DFT = {"V-Z-Pt": 3.07, "V-X-A": 4.21, "F-A": 6.20, "V-X-B": 6.40}  # A (paso5)

pose = Chem.SDMolSupplier("FA_pose1_caja18.sdf", removeHs=False)[0]
core = Chem.MolFromSmarts("[#7]~[#6]1~[#7]~[#6]2~[#7]~[#6]~[#6]~[#7]~[#6]2~[#6](~[#8])~[#7]~1")
idx = list(pose.GetSubstructMatch(core))
if len(idx) != 12:
    raise SystemExit("No encontre los 12 atomos de la pterina en la pose.")
X = pose.GetConformer().GetPositions()[idx]
c = X.mean(axis=0)                                   # centro de la pterina
nrm = np.linalg.svd(X - c)[2][2]                     # normal al plano de la pterina

atoms = list(PDBParser(QUIET=True).get_structure("p", "proteina.pdb").get_atoms())
R = np.array([a.coord for a in atoms], dtype=float)
print("Atomos de proteina: %d | contacto minimo con el metal: %.1f A" % (len(R), CONTACTO))

def altura(n):
    # n: normal de la superficie apuntando del metal hacia la proteina.
    # altura del centro de la pterina sobre la capa superior de Au cuando
    # el punto mas bajo de la proteina queda a CONTACTO de esa capa.
    return n @ c - (R @ n).min() + CONTACTO

# barrido de orientaciones (esfera de Fibonacci)
N = 60000
k = np.arange(N) + 0.5
phi = np.arccos(1 - 2 * k / N)
th = np.pi * (1 + 5 ** 0.5) * k
dirs = np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)
h = np.empty(N)
for s in range(0, N, 4000):
    D = dirs[s:s + 4000]
    h[s:s + 4000] = D @ c - (D @ R.T).min(axis=1) + CONTACTO

i = int(np.argmin(h))
ang = np.degrees(np.arccos(abs(dirs[i] @ nrm)))
print("\n--- Altura minima de la pterina sobre el Au para que FRalfa (rigida) no choque ---")
print("En CUALQUIER orientacion: %.1f A  (el plano de la pterina forma %.0f grados con la normal de la superficie)" % (h[i], ang))
print("Con la pterina paralela a la superficie (como en V-Z-Pt): %.1f A (cara +) | %.1f A (cara -)" % (altura(nrm), altura(-nrm)))
print("Alturas reales en tus geometrias DFT:", ", ".join("%s %.2f A" % kv for kv in ALTURAS_DFT.items()))
for nombre, hh in ALTURAS_DFT.items():
    print("  Fraccion de orientaciones posibles con altura <= %.2f A (%s): %.2f %%" % (hh, nombre, 100.0 * (h <= hh).mean()))
print("Altura minima en el 5%% mejor de las orientaciones: %.1f A" % np.percentile(h, 5))
