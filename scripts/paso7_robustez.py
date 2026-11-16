import os, re, subprocess, csv
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import PDBQTMolecule, RDKitMolCreate
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch

# ---------------- parametros (ajusta si hace falta) ----------------
PDB_FILE = "pdb4lrh.ent"
LIG = "FOL"
CHAINS = list(os.environ.get("CHAINS", "ABCDEFGH"))
SEEDS_A = [int(x) for x in os.environ.get("SEEDS", "1,2,3,4,5").split(",")]   # semillas extra, solo cadena A
CAJA = 18
EXH = int(os.environ.get("EXH", "32"))
LIGANDO_PDBQT = os.path.abspath("FA_dock_input.pdbqt")
CARPETA = "robustez"
SMI_NEUTRO = "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CCC(=O)O)C(=O)O)nc2c(=O)[nH]1"
CORE = Chem.MolFromSmarts("[#7]~[#6]1~[#7]~[#6]2~[#7]~[#6]~[#6]~[#7]~[#6]2~[#6](~[#8])~[#7]~1")

VINA = next((p for p in ("vina_1.2.7_win.exe", "vina_linux") if os.path.exists(p)), None)
if VINA is None:
    raise SystemExit("No encuentro vina_1.2.7_win.exe en esta carpeta.")
VINA = os.path.abspath(VINA)
os.makedirs(CARPETA, exist_ok=True)
P = lambda nombre: os.path.join(CARPETA, nombre)

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

def contactos(ns, xyz, d=4.0):
    s = set()
    for p in xyz:
        for a in ns.search(p, d):
            r = a.get_parent()
            s.add((r.get_resname(), r.id[1]))
    return s

def evaluar(out_pdbqt, refS, refxyz, core_ref, ns):
    txt = open(out_pdbqt).read()
    scores = [float(x) for x in re.findall(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", txt)]
    pm = PDBQTMolecule.from_file(out_pdbqt, skip_typing=True)
    probeS, probeH = esqueleto(RDKitMolCreate.from_pdbqt_mol(pm)[0])
    matches = refS.GetSubstructMatches(probeS, uniquify=False, useChirality=False, maxMatches=100000)
    if not matches:
        raise RuntimeError("sin correspondencia entre la pose y el ligando de referencia")
    cref = contactos(ns, refxyz[core_ref])
    filas = []
    for i in range(probeH.GetNumConformers()):
        pxyz = probeH.GetConformer(i).GetPositions()
        mejor = None
        for m in matches:
            tot = np.sqrt(np.mean(np.sum((refxyz[list(m)] - pxyz) ** 2, axis=1)))
            sel = [k for k, x in enumerate(m) if x in core_ref]
            pt = np.sqrt(np.mean(np.sum((refxyz[[m[k] for k in sel]] - pxyz[sel]) ** 2, axis=1)))
            if mejor is None or pt < mejor[1]:
                mejor = (tot, pt, sel)
        frac = len(contactos(ns, pxyz[mejor[2]]) & cref) / len(cref)
        filas.append((scores[i], mejor[0], mejor[1], frac))
    return filas

struct = PDBParser(QUIET=True).get_structure("x", PDB_FILE)
io = PDBIO()
io.set_structure(struct)
resumen = []

for ch in CHAINS:
    class Prot(Select):
        def accept_chain(self, c): return c.id == ch
        def accept_residue(self, r): return r.id[0] == " "
    class Lig(Select):
        def accept_chain(self, c): return c.id == ch
        def accept_residue(self, r): return r.resname == LIG
    io.save(P("prot_%s.pdb" % ch), Prot())
    io.save(P("lig_%s.pdb" % ch), Lig())

    pdbmol = Chem.MolFromPDBFile(P("lig_%s.pdb" % ch), removeHs=True, proximityBonding=True)
    if pdbmol is None:
        print("Cadena %s: sin ligando %s, la salto" % (ch, LIG)); continue
    try:
        ref = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(SMI_NEUTRO), pdbmol)
    except Exception as e:
        print("Cadena %s: no pude asignar enlaces al ligando (%s), la salto" % (ch, e)); continue
    refS, refH = esqueleto(ref)
    refxyz = refH.GetConformer().GetPositions()
    core_ref = list(refS.GetSubstructMatch(CORE))
    centro = refxyz.mean(axis=0)

    subprocess.run(["obabel", P("prot_%s.pdb" % ch), "-O", P("prot_%s.pdbqt" % ch), "-xr", "-p", "7.4"],
                   check=True, capture_output=True)
    ns = NeighborSearch(list(PDBParser(QUIET=True).get_structure("p", P("prot_%s.pdb" % ch)).get_atoms()))

    semillas = [42] + (SEEDS_A if ch == "A" else [])
    for sd in semillas:
        out = P("out_%s_s%d.pdbqt" % (ch, sd))
        cmd = [VINA, "--receptor", P("prot_%s.pdbqt" % ch), "--ligand", LIGANDO_PDBQT,
               "--center_x", "%.3f" % centro[0], "--center_y", "%.3f" % centro[1], "--center_z", "%.3f" % centro[2],
               "--size_x", str(CAJA), "--size_y", str(CAJA), "--size_z", str(CAJA),
               "--exhaustiveness", str(EXH), "--num_modes", "20", "--energy_range", "4",
               "--seed", str(sd), "--out", out]
        if not os.path.exists(out):          # si ya existe (corrida interrumpida), se reutiliza
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                print("Vina fallo en cadena %s semilla %d:\n%s" % (ch, sd, r.stderr[-500:])); continue
        f = evaluar(out, refS, refxyz, core_ref, ns)
        nat = [i + 1 for i, x in enumerate(f) if x[2] < 1.0 and x[3] >= 0.99]
        resumen.append((ch, sd, len(f), f[0][0], f[0][1], f[0][2], f[0][3], nat[0] if nat else 0, len(nat)))
        print("cadena %s semilla %2d | pose1: %.2f kcal/mol, RMSD total %.2f, pterina %.2f A, contactos %3.0f%% | primera nativa: rango %s | nativas: %d de %d" % (
            ch, sd, f[0][0], f[0][1], f[0][2], 100 * f[0][3], nat[0] if nat else "-", len(nat), len(f)))

with open("robustez.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["cadena", "semilla", "n_poses", "score_pose1", "rmsd_total_pose1", "rmsd_pterina_pose1",
                "contactos_pose1", "rango_primera_nativa", "n_nativas"])
    w.writerows(resumen)

n = len(resumen)
if n:
    print("\n=== RESUMEN (%d corridas) ===" % n)
    print("Pose 1 con RMSD total < 2 A:            %d de %d" % (sum(x[4] < 2.0 for x in resumen), n))
    print("Pose 1 con pterina < 1 A del cristal:   %d de %d" % (sum(x[5] < 1.0 for x in resumen), n))
    print("Corridas con alguna pose nativa en las 20: %d de %d" % (sum(x[7] > 0 for x in resumen), n))
    rangos = [x[7] for x in resumen if x[7] > 0]
    if rangos:
        print("Rango de la primera pose nativa: mediana %d, peor %d" % (np.median(rangos), max(rangos)))
    print("Guardado: robustez.csv")
