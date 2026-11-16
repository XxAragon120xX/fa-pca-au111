import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser

ANCLA = 2.2        # A: altura supuesta del O del carboxilato anclado sobre la capa superior de Au (supuesto)
N_CONF = 500       # conformeros para estimar la extension maxima
CORE = Chem.MolFromSmarts("[#7]~[#6]1~[#7]~[#6]2~[#7]~[#6]~[#6]~[#7]~[#6]2~[#6](~[#8])~[#7]~1")
SMI_ANION = "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CCC(=O)[O-])C(=O)[O-])nc2c(=O)[nH]1"

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

def partes(mol):
    S, H = esqueleto(mol)
    core = list(S.GetSubstructMatch(CORE))
    gamma, alfa = [], []
    for a in S.GetAtoms():
        nb = list(a.GetNeighbors())
        if a.GetSymbol() == "C" and len(nb) == 3:
            oxis = [x.GetIdx() for x in nb if x.GetSymbol() == "O"]
            cs = [x for x in nb if x.GetSymbol() == "C"]
            if len(oxis) == 2 and len(cs) == 1:
                es_alfa = any(y.GetSymbol() == "N" for y in cs[0].GetNeighbors())
                (alfa if es_alfa else gamma).extend(oxis)
    if len(core) != 12 or len(gamma) != 2 or len(alfa) != 2:
        raise SystemExit("No pude identificar la pterina o los carboxilatos (core %d, gamma %d, alfa %d)" % (len(core), len(gamma), len(alfa)))
    return core, gamma, alfa, H

def dist(H, i, core, oxis):
    X = H.GetConformer(i).GetPositions()
    c = X[core].mean(axis=0)
    return [float(np.linalg.norm(X[o] - c)) for o in oxis]

# ---- (a) pose docked (la misma que usamos en el analisis de oclusion) ----
pose = Chem.SDMolSupplier("FA_pose1_caja18.sdf", removeHs=False)[0]
core, gamma, alfa, H = partes(pose)
dg = dist(H, 0, core, gamma)
da = dist(H, 0, core, alfa)
print("--- Pose docked (extendida en el bolsillo) ---")
print("Centro de la pterina -> O del carboxilato gamma (A): %.1f y %.1f" % tuple(dg))
print("Centro de la pterina -> O del carboxilato alfa  (A): %.1f y %.1f" % tuple(da))

# ---- (b) todos los conformeros posibles: extension maxima ----
m = Chem.AddHs(Chem.MolFromSmiles(SMI_ANION))
cids = AllChem.EmbedMultipleConfs(m, numConfs=N_CONF, randomSeed=42)
core2, gamma2, alfa2, H2 = partes(m)
dmax = np.array([max(dist(H2, i, core2, gamma2)) for i in range(H2.GetNumConformers())])
print("\n--- %d conformeros generados (sin proteina) ---" % len(dmax))
print("Distancia centro de la pterina -> O gamma: mediana %.1f | percentil 95 %.1f | maxima %.1f A" % (
    np.median(dmax), np.percentile(dmax, 95), dmax.max()))

# ---- (c) altura minima exigida por FRalfa (misma logica que paso6) ----
atoms = list(PDBParser(QUIET=True).get_structure("p", "proteina.pdb").get_atoms())
R = np.array([a.coord for a in atoms], dtype=float)
c = H.GetConformer(0).GetPositions()[core].mean(axis=0)
N = 60000
k = np.arange(N) + 0.5
phi = np.arccos(1 - 2 * k / N)
th = np.pi * (1 + 5 ** 0.5) * k
dirs = np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)
base = np.empty(N)
for s in range(0, N, 4000):
    D = dirs[s:s + 4000]
    base[s:s + 4000] = D @ c - (D @ R.T).min(axis=1)

print("\n--- Altura minima de la pterina que exige FRalfa segun el contacto minimo con el metal ---")
for contacto in (2.5, 3.0, 3.5):
    print("  contacto %.1f A -> %.1f A" % (contacto, base.min() + contacto))

# ---- (d) comparacion: FA erguido y anclado por el carboxilato gamma ----
hmin = base.min() + 3.0
print("\n--- FA erguido, anclado por el carboxilato gamma (O a %.1f A sobre la capa de Au) ---" % ANCLA)
for nombre, d in (("pose docked", max(dg)), ("conformero mas extendido", dmax.max()), ("percentil 95 de conformeros", np.percentile(dmax, 95))):
    h = d + ANCLA
    print("  %-28s altura maxima de la pterina = %.1f A  (exigida: %.1f A) -> %s" % (
        nombre, h, hmin, "condicion necesaria CUMPLIDA" if h >= hmin else "NO alcanza"))
print("\nNota: es una cota superior (cadena totalmente erguida); que se cumpla es necesario, no suficiente.")
