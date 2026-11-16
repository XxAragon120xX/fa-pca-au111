import sys, re, os
import numpy as np
from scipy.spatial import cKDTree
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from Bio.PDB import PDBParser, PDBIO

QE_IN = sys.argv[1] if len(sys.argv) > 1 else "FA_Au111_flat_A_down.in"
POSE = "FA_pose1_caja18.sdf"
PROT = "proteina.pdb"
STEM = os.path.splitext(os.path.basename(QE_IN))[0]

# ---------- 1) leer el input de QE ----------
lines = open(QE_IN, encoding="utf-8", errors="ignore").read().splitlines()
cell = None
pos = []
for i, ln in enumerate(lines):
    s = ln.strip().upper()
    if s.startswith("CELL_PARAMETERS"):
        cell = np.array([[float(x) for x in lines[i + k].split()[:3]] for k in (1, 2, 3)])
    if s.startswith("ATOMIC_POSITIONS"):
        if "angstrom" not in ln.lower():
            print("AVISO: ATOMIC_POSITIONS no esta en angstrom:", ln)
        j = i + 1
        while j < len(lines) and re.match(r"^\s*[A-Za-z]{1,2}\s+-?\d", lines[j]):
            t = lines[j].split()
            pos.append((t[0], float(t[1]), float(t[2]), float(t[3])))
            j += 1
if cell is None or not pos:
    raise SystemExit("No pude leer CELL_PARAMETERS / ATOMIC_POSITIONS de " + QE_IN)

el = [p[0] for p in pos]
xyz = np.array([[p[1], p[2], p[3]] for p in pos])
is_au = np.array([e == "Au" for e in el])
au = xyz[is_au]
fa = xyz[~is_au]
fa_el = [e for e, a in zip(el, is_au) if not a]
print("Au: %d | atomos de FA: %d" % (len(au), len(fa)))
ztop = au[:, 2].max()
top = au[au[:, 2] > ztop - 1.0]
print("Capa superior de Au: z medio = %.2f A (%d atomos)" % (top[:, 2].mean(), len(top)))

# ---------- 2) FA del DFT y de la pose: correspondencia atomo a atomo ----------
def esqueleto(m):
    m.UpdatePropertyCache(strict=False)
    m = Chem.RemoveHs(m, sanitize=False)
    rw = Chem.RWMol(m)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0)
        a.SetIsAromatic(False)
        a.SetNumExplicitHs(0)
        a.SetNoImplicit(True)
    return rw.GetMol(), m

blk = "%d\nFA\n" % len(fa) + "\n".join(
    "%s %.6f %.6f %.6f" % (e, p[0], p[1], p[2]) for e, p in zip(fa_el, fa))
qe = Chem.MolFromXYZBlock(blk)
rdDetermineBonds.DetermineConnectivity(qe)
pose = Chem.SDMolSupplier(POSE, removeHs=False)[0]

poseS, poseH = esqueleto(pose)
qeS, qeH = esqueleto(qe)
print("Atomos pesados pose/DFT: %d / %d | enlaces: %d / %d" % (
    poseH.GetNumAtoms(), qeH.GetNumAtoms(), poseS.GetNumBonds(), qeS.GetNumBonds()))

matches = qeS.GetSubstructMatches(poseS, uniquify=False, useChirality=False, maxMatches=100000)
if not matches:
    raise SystemExit("No hay correspondencia entre la pose y el FA del DFT (revisa los enlaces).")
m = matches[0]

core = Chem.MolFromSmarts("[#7]~[#6]1~[#7]~[#6]2~[#7]~[#6]~[#6]~[#7]~[#6]2~[#6](~[#8])~[#7]~1")
core_pose = list(poseS.GetSubstructMatch(core))
if len(core_pose) != 12:
    raise SystemExit("No encontre los 12 atomos de la pterina en la pose.")

P = poseH.GetConformer().GetPositions()[core_pose]
Q = qeH.GetConformer().GetPositions()[[m[k] for k in core_pose]]

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Qc - R @ Pc

R, t = kabsch(P, Q)
rmsd = np.sqrt(np.mean(np.sum((P @ R.T + t - Q) ** 2, axis=1)))
print("RMSD del alineamiento de la pterina (docking vs DFT): %.2f A" % rmsd)
print("Altura media de la pterina sobre la capa superior de Au: %.2f A" % (Q[:, 2].mean() - top[:, 2].mean()))

# ---------- 3) mover la proteina y contar choques con el slab ----------
prot = PDBParser(QUIET=True).get_structure("p", PROT)
atoms = list(prot.get_atoms())
pxyz = np.array([a.coord for a in atoms], dtype=float) @ R.T + t
for a, c in zip(atoms, pxyz):
    a.set_coord(c.astype("float32"))

aus = np.vstack([au + i * cell[0] + j * cell[1] for i in range(-3, 4) for j in range(-3, 4)])
d, idx = cKDTree(aus).query(pxyz)
res = [(a.get_parent().get_resname(), a.get_parent().id[1]) for a in atoms]

def lista(mask):
    s = sorted({res[i] for i in np.where(mask)[0]}, key=lambda x: x[1])
    return ", ".join("%s%d" % (n, k) for n, k in s) if s else "ninguno"

iz = int(np.argmin(d))
print("\n--- Proteina vs slab de Au (%d atomos de proteina) ---" % len(atoms))
print("Distancia minima proteina-Au: %.2f A (%s%d, atomo %s)" % (
    d[iz], res[iz][0], res[iz][1], atoms[iz].get_id()))
print("Atomos con d(Au) < 2.5 A (choque duro): %d" % int((d < 2.5).sum()))
print("Atomos con d(Au) < 3.5 A (contacto):    %d" % int((d < 3.5).sum()))
print("Atomos por debajo de la capa superior de Au (dentro del metal): %d" % int((pxyz[:, 2] < ztop).sum()))
print("Altura del punto mas bajo de la proteina sobre la capa superior: %.2f A" % (pxyz[:, 2].min() - ztop))
print("Residuos en choque duro (<2.5 A):", lista(d < 2.5))
print("Residuos a menos de 3.5 A:", lista(d < 3.5))

fah = qeH.GetConformer().GetPositions()
dfa, _ = cKDTree(pxyz).query(fah)
print("\n--- FA adsorbido vs proteina (solo la pterina esta alineada) ---")
print("Distancia minima FA-proteina: %.2f A | atomos del FA a menos de 2.0 A de la proteina: %d de %d" % (
    dfa.min(), int((dfa < 2.0).sum()), len(dfa)))

# ---------- 4) archivos para visualizar ----------
io = PDBIO()
io.set_structure(prot)
io.save(STEM + "_proteina.pdb")

def pdb_line(n, name, resn, resi, x, y, z, elem):
    return "HETATM%5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s" % (
        n % 100000, name, resn, resi % 10000, x, y, z, elem)

out = []
n = 0
for i in (-1, 0, 1):
    for j in (-1, 0, 1):
        for p in au + i * cell[0] + j * cell[1]:
            n += 1
            out.append(pdb_line(n, "AU", "AUS", 1, p[0], p[1], p[2], "Au"))
for k, (e, p) in enumerate(zip(fa_el, fa)):
    n += 1
    out.append(pdb_line(n, ("%s%d" % (e, k + 1))[:4], "FA", 2, p[0], p[1], p[2], e))
open(STEM + "_Au_FA.pdb", "w").write("\n".join(out) + "\nEND\n")
print("\nListo: %s_proteina.pdb y %s_Au_FA.pdb (abrelos juntos en PyMOL o ChimeraX)" % (STEM, STEM))
print("\nRESUMEN | %s | altura pterina %.2f A | atomos en metal %d de %d (%.0f%%) | choque duro %d | d_min %.2f A" % (
    STEM, Q[:, 2].mean() - top[:, 2].mean(), int((pxyz[:, 2] < ztop).sum()), len(atoms),
    100.0 * (pxyz[:, 2] < ztop).sum() / len(atoms), int((d < 2.5).sum()), d.min()))
