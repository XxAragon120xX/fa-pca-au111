#!/usr/bin/env python3
"""
analyze_pca.py — Análisis electrónico de adsorbato/Au(111) a partir de las salidas
de post_electronic.sh (Quantum ESPRESSO). Reproduce la Sección 3.2 usada para FA:

  1. DOS total de A (complejo), B (slab) y C (molécula)           -> fig_dos.png
  2. PDOS por fragmentos (pterina, COOH, capa Au superior),
     picos HOMO/LUMO, estados en el gap                          -> fig_pdos.png
  3. Bader: Δq slab, por capa, por átomo                         -> bader_*.csv
  4. Δρ(z): integral, acumulación/depleción, carga desplazada    -> fig_drho_z.png
  5. Función trabajo: ajuste lineal del vacío (sin corrección dipolar) -> fig_potential.png

Uso (desde ~/Au/SCF_final):
    python analyze_pca.py PCA_*            # analiza cada carpeta + summary_PCA.csv
    python analyze_pca.py PCA_Au111_flat_A_down

Las salidas por sistema quedan en <carpeta>/analysis/.
Requiere numpy y matplotlib (scipy opcional).
"""
import sys, os, re, glob, json, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RY = 13.605693122994          # eV
BOHR = 0.529177210903         # Å
ZSYM = {1: "H", 6: "C", 7: "N", 8: "O", 79: "Au"}
ZVAL_DEFAULT = {"Au": 11.0, "C": 4.0, "N": 5.0, "O": 6.0, "H": 1.0}

# ----------------------------- parámetros -----------------------------------
GAP_WIN = (0.4, 1.0)      # ventana de estados en el gap respecto a E_F(A), eV
VAC_MARGIN_A = 3.0        # Å excluidos a cada lado del vacío para el ajuste
LAYER_TOL = 0.5           # Å, tolerancia para agrupar capas de Au
PEAK_REL = 0.05           # prominencia mínima de picos (fracción del máximo)
E_PLOT = (-7.0, 4.0)      # rango de energía de las figuras (E − E_F, eV)


def log(*a):
    print(*a, flush=True)


# ============================ lectores ======================================
def read_cube_header(path):
    """Átomos (Z, xyz en bohr), origen, vectores de voxel y forma de la malla."""
    with open(path) as f:
        f.readline(); f.readline()
        t = f.readline().split()
        nat = abs(int(t[0])); origin = np.array(t[1:4], float)
        n, vox = [], []
        for _ in range(3):
            t = f.readline().split()
            n.append(int(t[0])); vox.append([float(x) for x in t[1:4]])
        Z, xyz = [], []
        for _ in range(nat):
            t = f.readline().split()
            Z.append(int(float(t[0]))); xyz.append([float(x) for x in t[2:5]])
    vox = np.array(vox)
    cell = vox * np.array(n)[:, None]
    return dict(Z=np.array(Z), xyz=np.array(xyz), origin=origin, n=np.array(n),
                vox=vox, cell=cell, sym=[ZSYM.get(z, str(z)) for z in Z])


def read_acf(path):
    pop, tail = [], {}
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) >= 5 and t[0].isdigit():
                pop.append(float(t[4]))
            elif ":" in line:
                k, v = line.split(":", 1)
                try:
                    tail[k.strip()] = float(v.split()[0])
                except (ValueError, IndexError):
                    pass
    return np.array(pop), tail


def read_pwout(path):
    d = dict(efermi=None, zval={}, nelec=None, etot=None)
    if not os.path.isfile(path):
        return d
    txt = open(path, errors="ignore").read()
    m = re.findall(r"the Fermi energy is\s+([-\d.]+)\s+ev", txt)
    if m:
        d["efermi"] = float(m[-1])
    m = re.findall(r"number of electrons\s*=\s*([\d.]+)", txt)
    if m:
        d["nelec"] = float(m[0])
    m = re.findall(r"!\s+total energy\s+=\s+([-\d.]+)\s+Ry", txt)
    if m:
        d["etot"] = float(m[-1])
    blk = re.search(r"atomic species\s+valence\s+mass.*?\n(.*?)\n\s*\n", txt, re.S)
    if blk:
        for line in blk.group(1).splitlines():
            t = line.split()
            if len(t) >= 2:
                try:
                    d["zval"][re.sub(r"[\d_-].*$", "", t[0])] = float(t[1])
                except ValueError:
                    pass
    return d


def read_dos(path):
    ef = None
    with open(path) as f:
        head = f.readline()
    m = re.search(r"EFermi\s*=\s*([-\d.]+)", head)
    if m:
        ef = float(m.group(1))
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 1], data[:, 2], ef


def read_pdos_atoms(prefix_glob):
    """{índice_atómico(1-based): (E, ldos)} sumando todos los wfc de cada átomo."""
    out = {}
    for fn in glob.glob(prefix_glob + ".pdos_atm#*"):
        m = re.search(r"pdos_atm#(\d+)\(", fn)
        if not m:
            continue
        i = int(m.group(1))
        dat = np.loadtxt(fn, comments="#")
        E, ld = dat[:, 0], dat[:, 1]
        if i in out:
            out[i] = (E, out[i][1] + ld)
        else:
            out[i] = (E, ld)
    return out


def read_spilling(path):
    if not os.path.isfile(path):
        return None
    m = re.findall(r"Spilling Parameter:\s*([-\d.Ee+]+)", open(path, errors="ignore").read())
    return float(m[-1]) if m else None


def read_avg(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1], d[:, 2]          # z (bohr), planar, macroscópico


# ============================ geometría =====================================
def map_atoms(src, ref, cell, tol=0.05):
    """Para cada átomo de src, índice del átomo de ref en la misma posición (bohr, PBC)."""
    inv = np.linalg.inv(cell.T)
    idx = []
    for r in src:
        d = ref - r
        f = d @ inv.T
        f -= np.round(f)
        dist = np.linalg.norm(f @ cell, axis=1)
        j = int(np.argmin(dist))
        idx.append(j if dist[j] < tol else -1)
    return np.array(idx)


def fragments(geo):
    """Índices (0-based) de Au, molécula, COOH, pterina, capa superior/inferior."""
    sym = np.array(geo["sym"]); xyz = geo["xyz"]
    au = np.where(sym == "Au")[0]; mol = np.where(sym != "Au")[0]
    cooh = []
    if len(mol):
        P = xyz[mol]; S = sym[mol]
        D = np.linalg.norm(P[:, None] - P[None], axis=2) * BOHR
        for a in range(len(mol)):
            if S[a] != "C":
                continue
            Os = [b for b in range(len(mol)) if S[b] == "O" and D[a, b] < 1.45]
            if len(Os) == 2:
                grp = [a] + Os
                for o in Os:
                    grp += [h for h in range(len(mol)) if S[h] == "H" and D[o, h] < 1.15]
                cooh = sorted(mol[grp])
                break
    pter = sorted(set(mol) - set(cooh))
    top, bot, layers = [], [], []
    if len(au):
        z = xyz[au, 2] * BOHR
        order = np.argsort(z)
        cur = [order[0]]
        for k in order[1:]:
            if z[k] - z[cur[-1]] < LAYER_TOL:
                cur.append(k)
            else:
                layers.append(au[cur]); cur = [k]
        layers.append(au[cur])
        zmol = xyz[mol, 2].mean() * BOHR if len(mol) else z.max() + 1
        if zmol > z.mean():
            top, bot = layers[-1], layers[0]
        else:
            top, bot = layers[0], layers[-1]
            layers = layers[::-1]          # siempre: capas[0]=inferior ... capas[-1]=superior
    return dict(au=au, mol=mol, cooh=np.array(cooh, int), pter=np.array(pter, int),
                top=np.array(top, int), bot=np.array(bot, int), layers=layers)


def vacuum_fit(z, v, atoms_z, L, margin_bohr):
    """Ajuste lineal de V(z) en el hueco de vacío más grande. Devuelve dict con bordes."""
    za = np.sort(np.mod(atoms_z, L))
    gaps = np.diff(np.concatenate([za, [za[0] + L]]))
    k = int(np.argmax(gaps))
    z_lo, z_hi = za[k], za[k] + gaps[k]          # z_hi puede ser > L (envuelto)
    a, b = z_lo + margin_bohr, z_hi - margin_bohr
    if b - a < 2.0:
        a, b = z_lo + 0.3 * gaps[k], z_hi - 0.3 * gaps[k]
    zu = np.concatenate([z, z + L])
    vu = np.concatenate([v, v])
    sel = (zu >= a) & (zu <= b)
    p = np.polyfit(zu[sel], vu[sel], 1)
    return dict(z_lo=z_lo, z_hi=z_hi, a=a, b=b, slope=p[0], icpt=p[1],
                f=lambda zz: np.polyval(p, zz), sel_z=zu[sel], sel_v=vu[sel],
                center=np.polyval(p, 0.5 * (z_lo + z_hi)))


def find_peaks_simple(E, y, rel=PEAK_REL):
    try:
        from scipy.signal import find_peaks
        idx, _ = find_peaks(y, prominence=rel * y.max())
        return idx
    except ImportError:
        thr = rel * y.max()
        return np.array([i for i in range(1, len(y) - 1)
                         if y[i] >= y[i - 1] and y[i] > y[i + 1] and y[i] > thr])


def interp_cross(x, y, target):
    """Primer x donde y (monótona creciente) cruza target."""
    i = np.searchsorted(y, target)
    if i <= 0 or i >= len(y):
        return None
    return x[i - 1] + (target - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1])


# ============================ análisis ======================================
def analyze(folder):
    log(f"\n================ {folder} ================")
    os.makedirs(os.path.join(folder, "analysis"), exist_ok=True)
    P = lambda *p: os.path.join(folder, *p)
    OUT = lambda *p: os.path.join(folder, "analysis", *p)
    R = dict(system=os.path.basename(os.path.normpath(folder)))
    warn = []

    pw = {s: read_pwout(P(f"{s}.out")) for s in "ABC"}
    EF = {s: pw[s]["efermi"] for s in "ABC"}
    for s in "ABC":
        R[f"EF_{s}_eV"] = EF[s]
        R[f"Etot_{s}_Ry"] = pw[s]["etot"]
    if all(pw[s]["etot"] is not None for s in "ABC"):
        R["Eint_frozen_eV"] = (pw["A"]["etot"] - pw["B"]["etot"] - pw["C"]["etot"]) * RY
    zval = dict(ZVAL_DEFAULT); zval.update(pw["A"]["zval"])

    # ---------- geometría desde los cubes ----------
    geo = {s: read_cube_header(P(f"{s}_rho.cube")) if os.path.isfile(P(f"{s}_rho.cube")) else None
           for s in "ABC"}
    if geo["A"] is None:
        log("  ERROR: falta A_rho.cube; no se puede continuar con este sistema")
        return R
    gA = geo["A"]; cell = gA["cell"]; L = cell[2, 2]
    fr = fragments(gA)
    log(f"  átomos: {len(gA['Z'])}  Au={len(fr['au'])}  molécula={len(fr['mol'])}  "
        f"COOH={[int(i) + 1 for i in fr['cooh']]}  capas Au={len(fr['layers'])}")
    if len(fr["cooh"]) != 4:
        warn.append(f"COOH detectado con {len(fr['cooh'])} átomos (esperado 4)")
    mapB = map_atoms(gA["xyz"][fr["au"]], geo["B"]["xyz"], cell) if geo["B"] else None
    mapC = map_atoms(gA["xyz"][fr["mol"]], geo["C"]["xyz"], cell) if geo["C"] else None
    if mapB is not None and (mapB < 0).any():
        warn.append(f"{(mapB < 0).sum()} átomos de Au de A no coinciden con B (¿geometría distinta?)")
    if mapC is not None and (mapC < 0).any():
        warn.append(f"{(mapC < 0).sum()} átomos de la molécula de A no coinciden con C")
    z_top = gA["xyz"][fr["top"], 2].mean()
    R["n_Au_layers"] = len(fr["layers"])
    zmol = gA["xyz"][fr["mol"], 2]
    R["h_min_mol_above_top_A"] = float((zmol.min() - z_top) * BOHR)

    # ---------- spilling ----------
    for s in "ABC":
        R[f"spilling_{s}"] = read_spilling(P(f"projwfc_{s}.out"))

    # ---------- función trabajo / vacío ----------
    vac = {}
    for s in "ABC":
        f = P(f"{s}_vavg.dat")
        if not (os.path.isfile(f) and geo[s]):
            warn.append(f"falta {s}_vavg.dat o {s}_rho.cube: sin vacío para {s}")
            continue
        z, vp, vm = read_avg(f)
        vp_eV = vp * RY
        vf = vacuum_fit(z, vp_eV, geo[s]["xyz"][:, 2], geo[s]["cell"][2, 2], VAC_MARGIN_A / BOHR)
        vac[s] = dict(z=z, v=vp_eV, vm=vm * RY, fit=vf)
        R[f"vac_slope_{s}_eV_per_A"] = vf["slope"] / BOHR
        R[f"Vvac_center_{s}_eV"] = vf["center"]

    def side_values(s, g, f_ids):
        """Vacío del lado cubierto (molécula) y del lado limpio para A/B."""
        vf = vac[s]["fit"]
        Ls = g["cell"][2, 2]
        za = np.mod(g["xyz"][:, 2], Ls)
        # el átomo justo debajo del vacío es z_lo; el de arriba (envuelto) es z_hi
        i_lo = int(np.argmin(np.abs(za - np.mod(vf["z_lo"], Ls))))
        i_hi = int(np.argmin(np.abs(za - np.mod(vf["z_hi"], Ls))))
        return dict(lo=dict(atom=i_lo, edge=vf["f"](vf["a"]), plane=vf["f"](vf["z_lo"])),
                    hi=dict(atom=i_hi, edge=vf["f"](vf["b"]), plane=vf["f"](vf["z_hi"])))

    if "A" in vac and EF["A"] is not None:
        sv = side_values("A", gA, fr)
        molset = set(fr["mol"].tolist())
        cov = "lo" if sv["lo"]["atom"] in molset else ("hi" if sv["hi"]["atom"] in molset else None)
        if cov is None:
            warn.append("la molécula no está en el borde del vacío; lado cubierto ambiguo")
            cov = "lo" if zmol.mean() > gA["xyz"][fr["au"], 2].mean() else "hi"
        cln = "hi" if cov == "lo" else "lo"
        R["WF_A_covered_edge_eV"] = sv[cov]["edge"] - EF["A"]
        R["WF_A_covered_extrap_eV"] = sv[cov]["plane"] - EF["A"]
        R["WF_A_clean_edge_eV"] = sv[cln]["edge"] - EF["A"]
        R["WF_A_clean_extrap_eV"] = sv[cln]["plane"] - EF["A"]
        R["Evac_ref_A_eV"] = sv[cln]["plane"]         # vacío de la cara limpia, extrapolado (como FA)
    if "B" in vac and EF["B"] is not None and geo["B"]:
        sv = side_values("B", geo["B"], None)
        vals = [sv["lo"]["edge"] - EF["B"], sv["hi"]["edge"] - EF["B"]]
        R["WF_B_side1_eV"], R["WF_B_side2_eV"] = vals
    if "C" in vac:
        R["Evac_ref_C_eV"] = vac["C"]["fit"]["center"]

    # ---------- DOS ----------
    dos = {}
    for s in "ABC":
        if os.path.isfile(P(f"{s}.dos")):
            E, d, idos, ef = read_dos(P(f"{s}.dos"))
            dos[s] = (E, d, idos, ef if ef is not None else EF[s])
    if dos:
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.6), sharey=False)
        names = {"A": "PCA/Au(111)", "B": "Au(111) (congelado)", "C": "PCA aislada (congelada)"}
        for ax, s in zip(axs, "BCA"):
            if s in dos:
                E, d, _, ef = dos[s]
                ax.plot(E - ef, d, lw=0.8, color="k")
                ax.axvline(0, ls="--", lw=0.8, color="gray")
                ax.set_xlim(*E_PLOT)
            ax.set_title(names[s]); ax.set_xlabel(r"$E-E_F$ (eV)")
        axs[0].set_ylabel("DOS (estados/eV)")
        fig.tight_layout(); fig.savefig(OUT("fig_dos.png"), dpi=300); plt.close(fig)

    # ---------- PDOS por fragmentos ----------
    pd = {s: read_pdos_atoms(P(f"pdos_{s}", s)) for s in "ABC"}
    frag_pdos = {}
    if pd["A"]:
        E = next(iter(pd["A"].values()))[0]
        def fsum(dct, ids):
            return sum((dct[i + 1][1] for i in ids if (i + 1) in dct), np.zeros_like(E))
        frag_pdos["A"] = dict(E=E, pter=fsum(pd["A"], fr["pter"]), cooh=fsum(pd["A"], fr["cooh"]),
                              top=fsum(pd["A"], fr["top"]))
        frag_pdos["A"]["mol"] = frag_pdos["A"]["pter"] + frag_pdos["A"]["cooh"]
    if pd["C"] and mapC is not None:
        EC = next(iter(pd["C"].values()))[0]
        cp = [int(mapC[list(fr["mol"]).index(i)]) for i in fr["pter"]]
        cc = [int(mapC[list(fr["mol"]).index(i)]) for i in fr["cooh"]]
        def fsumC(ids):
            return sum((pd["C"][i + 1][1] for i in ids if i >= 0 and (i + 1) in pd["C"]), np.zeros_like(EC))
        frag_pdos["C"] = dict(E=EC, pter=fsumC(cp), cooh=fsumC(cc))
        frag_pdos["C"]["mol"] = frag_pdos["C"]["pter"] + frag_pdos["C"]["cooh"]

    def homo_lumo_peaks(E, y, ef, window=4.0):
        sel = (E > ef - window) & (E < ef + window)
        Es, ys = E[sel], y[sel]
        if len(ys) < 5 or ys.max() <= 0:
            return None, None
        pk = find_peaks_simple(Es, ys)
        occ = [Es[i] for i in pk if Es[i] < ef]
        unocc = [Es[i] for i in pk if Es[i] > ef]
        return (max(occ) if occ else None), (min(unocc) if unocc else None)

    if "A" in frag_pdos and EF["A"] is not None:
        h, l = homo_lumo_peaks(frag_pdos["A"]["E"], frag_pdos["A"]["mol"], EF["A"])
        R["HOMO_A_minus_EF_eV"] = None if h is None else h - EF["A"]
        R["LUMO_A_minus_EF_eV"] = None if l is None else l - EF["A"]
        if h is not None and l is not None:
            R["gap_A_eV"] = l - h
        if "Evac_ref_A_eV" in R:
            R["HOMO_A_vac_eV"] = None if h is None else h - R["Evac_ref_A_eV"]
            R["LUMO_A_vac_eV"] = None if l is None else l - R["Evac_ref_A_eV"]
            R["EF_A_vac_eV"] = EF["A"] - R["Evac_ref_A_eV"]
    if "C" in frag_pdos and EF["C"] is not None:
        # HOMO/LUMO del aislado por conteo de electrones en la DOS integrada (más robusto)
        hC = lC = None
        if "C" in dos and pw["C"]["nelec"]:
            E, d, idos, _ = dos["C"]
            N = pw["C"]["nelec"]
            hC, lC = interp_cross(E, idos, N - 1.0), interp_cross(E, idos, N + 1.0)
        if hC is None or lC is None:
            hC, lC = homo_lumo_peaks(frag_pdos["C"]["E"], frag_pdos["C"]["mol"], EF["C"])
        R["HOMO_C_minus_EF_eV"] = None if hC is None else hC - EF["C"]
        R["LUMO_C_minus_EF_eV"] = None if lC is None else lC - EF["C"]
        if hC is not None and lC is not None:
            R["gap_C_eV"] = lC - hC
        if "Evac_ref_C_eV" in R:
            R["HOMO_C_vac_eV"] = None if hC is None else hC - R["Evac_ref_C_eV"]
            R["LUMO_C_vac_eV"] = None if lC is None else lC - R["Evac_ref_C_eV"]

    # estados en el gap (ventana sobre E_F(A)); para C la misma ventana alineada al vacío
    if "A" in frag_pdos and EF["A"] is not None:
        E = frag_pdos["A"]["E"]
        w = (E >= EF["A"] + GAP_WIN[0]) & (E <= EF["A"] + GAP_WIN[1])
        dE = np.mean(np.diff(E))
        ip = frag_pdos["A"]["pter"][w].sum() * dE
        ic = frag_pdos["A"]["cooh"][w].sum() * dE
        R["gapstates_pter_A"] = ip
        R["gapstates_cooh_A"] = ic
        npt, nco = max(len(fr["pter"]), 1), max(len(fr["cooh"]), 1)
        R["gapstates_ratio_per_atom_pter_over_cooh"] = (ip / npt) / (ic / nco) if ic > 0 else None
        if "C" in frag_pdos and "Evac_ref_A_eV" in R and "Evac_ref_C_eV" in R:
            EC = frag_pdos["C"]["E"]
            lo = EF["A"] + GAP_WIN[0] - R["Evac_ref_A_eV"] + R["Evac_ref_C_eV"]
            hi = EF["A"] + GAP_WIN[1] - R["Evac_ref_A_eV"] + R["Evac_ref_C_eV"]
            wc = (EC >= lo) & (EC <= hi)
            R["gapstates_pter_C_same_vac_window"] = frag_pdos["C"]["pter"][wc].sum() * np.mean(np.diff(EC))

    if frag_pdos:
        fig, axs = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
        if "A" in frag_pdos:
            fa = frag_pdos["A"]; x = fa["E"] - EF["A"]
            axs[0].plot(x, fa["pter"], label="pterina", lw=1)
            axs[0].plot(x, fa["cooh"], label="COOH", lw=1)
            axs[0].plot(x, fa["top"] / max(len(fr["top"]), 1) * 5, label="Au capa sup. (×5/átomo)",
                        lw=0.8, color="goldenrod", alpha=0.7)
            axs[0].axvspan(*GAP_WIN, color="gray", alpha=0.15)
            axs[0].set_title("PCA/Au(111) — PDOS"); axs[0].legend(fontsize=8)
        if "C" in frag_pdos:
            fc = frag_pdos["C"]; x = fc["E"] - EF["C"]
            axs[1].plot(x, fc["pter"], label="pterina", lw=1)
            axs[1].plot(x, fc["cooh"], label="COOH", lw=1)
            axs[1].set_title("PCA aislada — PDOS"); axs[1].legend(fontsize=8)
        for ax in axs:
            ax.axvline(0, ls="--", lw=0.8, color="gray"); ax.set_ylabel("PDOS (estados/eV)")
        axs[1].set_xlim(*E_PLOT); axs[1].set_xlabel(r"$E-E_F$ (eV)")
        fig.tight_layout(); fig.savefig(OUT("fig_pdos.png"), dpi=300); plt.close(fig)

    # ---------- Bader ----------
    acf = {s: read_acf(P(f"bader_{s}", "ACF.dat")) if os.path.isfile(P(f"bader_{s}", "ACF.dat")) else None
           for s in "ABC"}
    if acf["A"] is not None:
        popA, tA = acf["A"]
        nom = sum(zval.get(s, 0) for s in gA["sym"])
        R["bader_Ne_A"] = float(popA.sum()); R["bader_Ne_err_A"] = float(popA.sum() - nom)
        rows = []
        qA = np.array([zval.get(s, 0) for s in gA["sym"]]) - popA
        slabA = popA[fr["au"]].sum()
        R["bader_slab_pop_A"] = float(slabA)
        if acf["B"] is not None and mapB is not None and (mapB >= 0).all():
            popB, _ = acf["B"]
            R["bader_Ne_err_B"] = float(popB.sum() - 11.0 * len(popB))
            R["bader_slab_pop_B"] = float(popB.sum())
            dpop_au = popA[fr["au"]] - popB[mapB]        # >0: el Au gana electrones
            R["dq_mol_to_Au_e"] = float(dpop_au.sum())
            pos = {a: k for k, a in enumerate(fr["au"])}
            R["dq_top_layer_e"] = float(sum(dpop_au[pos[a]] for a in fr["top"]))
            R["dq_bottom_layer_e"] = float(sum(dpop_au[pos[a]] for a in fr["bot"]))
            topd = np.array([dpop_au[pos[a]] for a in fr["top"]])
            R["dq_top_atom_max_gain_e"] = float(topd.max()); R["dq_top_atom_max_loss_e"] = float(topd.min())
            for li, lay in enumerate(fr["layers"]):
                R[f"dq_layer{li + 1}_e"] = float(sum(dpop_au[pos[a]] for a in lay))
            for k, a in enumerate(fr["au"]):
                rows.append(dict(atom=a + 1, el="Au", z_A=gA["xyz"][a, 2] * BOHR,
                                 pop_A=popA[a], pop_ref=popB[mapB[k]], dpop=dpop_au[k],
                                 q_A=qA[a]))
        if acf["C"] is not None and mapC is not None and (mapC >= 0).all():
            popC, _ = acf["C"]
            symC = geo["C"]["sym"]
            nomC = sum(zval.get(s, 0) for s in symC)
            R["bader_Ne_err_C"] = float(popC.sum() - nomC)
            dpop_m = popA[fr["mol"]] - popC[mapC]
            R["dq_mol_total_e"] = float(-dpop_m.sum())   # >0: la molécula pierde electrones
            R["dq_mol_atom_maxabs_e"] = float(np.abs(dpop_m).max())
            R["dq_pter_e"] = float(-sum(dpop_m[list(fr["mol"]).index(i)] for i in fr["pter"]))
            R["dq_cooh_e"] = float(-sum(dpop_m[list(fr["mol"]).index(i)] for i in fr["cooh"]))
            for k, a in enumerate(fr["mol"]):
                dmin = np.min(np.linalg.norm(
                    (gA["xyz"][fr["top"]] - gA["xyz"][a]), axis=1)) * BOHR
                rows.append(dict(atom=a + 1, el=gA["sym"][a], z_A=gA["xyz"][a, 2] * BOHR,
                                 pop_A=popA[a], pop_ref=popC[mapC[k]], dpop=dpop_m[k],
                                 q_A=qA[a], d_to_top_Au_A=dmin,
                                 frag="COOH" if a in fr["cooh"] else "pterina"))
        if rows:
            keys = ["atom", "el", "frag", "z_A", "d_to_top_Au_A", "pop_A", "pop_ref", "dpop", "q_A"]
            with open(OUT("bader_atoms.csv"), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
                w.writeheader()
                for r in rows:
                    w.writerow({k: (f"{v:.5f}" if isinstance(v, float) else v) for k, v in r.items()})

    # ---------- Δρ(z) ----------
    if os.path.isfile(P("drho_zavg.dat")):
        z, pl, ma = read_avg(P("drho_zavg.dat"))
        area = np.linalg.norm(np.cross(cell[0], cell[1]))          # bohr²
        dz = np.mean(np.diff(z))
        R["drho_integral_e"] = float(pl.sum() * dz * area)
        lam = pl * area                                            # e/bohr
        zt = z_top
        zA = (z - zt) * BOHR                                       # Å respecto a la capa superior
        # carga por encima de cada plano (desde la mitad del vacío hacia abajo)
        vf = vacuum_fit(z, pl, gA["xyz"][:, 2], L, VAC_MARGIN_A / BOHR)
        zc = np.mod(0.5 * (vf["z_lo"] + vf["z_hi"]), L)
        side = 1 if zmol.mean() > zt else -1
        if side > 0:
            # Q(z) = ∫_z^{zc} λ dz' (envolviendo si zc < z_top)
            zz = np.where(z < zt - 2.0 / BOHR, z + L, z)
            order = np.argsort(zz); zs, ls_ = zz[order], lam[order]
            zcu = zc if zc > zt else zc + L
            m = zs <= zcu
            Q = np.cumsum(ls_[m][::-1])[::-1] * dz
            zQ = (zs[m] - zt) * BOHR
        else:
            zz = np.where(z > zt + 2.0 / BOHR, z - L, z)
            order = np.argsort(zz); zs, ls_ = zz[order], lam[order]
            zcu = zc if zc < zt else zc - L
            m = zs >= zcu
            Q = np.cumsum(ls_[m]) * dz
            zQ = (zt - zs[m]) * BOHR
        win = (zQ > -0.5) & (zQ < (zmol.max() - zt) * side * BOHR + 3.0)
        if win.any():
            R["drho_displaced_charge_e"] = float(np.max(np.abs(Q[win])))
        # alturas de acumulación/depleción máximas en la interfaz (entre la capa superior y la molécula)
        hmax = abs(zmol.min() - zt) * BOHR if side > 0 else abs(zt - zmol.max()) * BOHR
        hz = (z - zt) * BOHR * side
        iw = (hz > 0.3) & (hz < hmax + 1.0)
        if iw.any():
            ii = np.where(iw)[0]
            R["drho_accum_height_A"] = float(hz[ii[np.argmax(pl[ii])]])
            R["drho_deplet_height_A"] = float(hz[ii[np.argmin(pl[ii])]])
        fig, ax = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
        ax[0].plot(zA * side, lam, lw=1); ax[0].axhline(0, color="gray", lw=0.6)
        ax[0].set_ylabel(r"$\Delta\lambda(z)$ (e/bohr)")
        ax[1].plot(zQ, Q, lw=1, color="C3"); ax[1].axhline(0, color="gray", lw=0.6)
        ax[1].set_ylabel("Q por encima de z (e)")
        ax[1].set_xlabel("z respecto a la capa de Au superior (Å)")
        for a in ax:
            a.axvline(0, ls="--", color="goldenrod", lw=0.8)
            a.set_xlim(-12, 15)
        fig.tight_layout(); fig.savefig(OUT("fig_drho_z.png"), dpi=300); plt.close(fig)
        np.savetxt(OUT("drho_z.dat"), np.c_[zA * side, lam, pl],
                   header="z_rel_top_Au(A)  lambda(e/bohr)  planar_avg(e/bohr^3)")

    # ---------- figura de potencial ----------
    if vac:
        fig, ax = plt.subplots(len(vac), 1, figsize=(7, 2.6 * len(vac)), squeeze=False)
        for a, s in zip(ax[:, 0], vac):
            v = vac[s]; vf = v["fit"]
            ef = EF[s] if EF[s] is not None else 0
            a.plot(v["z"] * BOHR, v["v"] - ef, lw=0.8, label="planar")
            a.plot(v["z"] * BOHR, v["vm"] - ef, lw=0.8, label="macroscópico")
            zz = np.linspace(vf["z_lo"], vf["z_hi"], 50)
            Ls = geo[s]["cell"][2, 2]
            a.plot(np.mod(zz, Ls) * BOHR, vf["f"](zz) - ef, "k.", ms=2, label="ajuste vacío")
            a.set_ylabel(r"$V-E_F$ (eV)"); a.set_title(s); a.legend(fontsize=7)
        ax[-1, 0].set_xlabel("z (Å)")
        fig.tight_layout(); fig.savefig(OUT("fig_potential.png"), dpi=300); plt.close(fig)

    R["warnings"] = "; ".join(warn)
    with open(OUT("summary.json"), "w") as f:
        json.dump(R, f, indent=2, default=lambda o: float(o))
    for k, v in R.items():
        if isinstance(v, float):
            log(f"  {k:40s} {v: .4f}")
        elif v not in (None, ""):
            log(f"  {k:40s} {v}")
    for w_ in warn:
        log(f"  AVISO: {w_}")
    return R


def main():
    folders = [f for f in sys.argv[1:] if os.path.isdir(f)] or ["."]
    allR = []
    for f in folders:
        try:
            allR.append(analyze(f))
        except Exception as e:
            log(f"  ERROR en {f}: {type(e).__name__}: {e}")
            allR.append(dict(system=f, warnings=f"ERROR {e}"))
    if len(folders) > 1:
        keys = []
        for r in allR:
            keys += [k for k in r if k not in keys]
        with open("summary_PCA.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            for r in allR:
                w.writerow({k: (f"{v:.5f}" if isinstance(v, float) else v) for k, v in r.items()})
        log("\nResumen de todas las configuraciones -> summary_PCA.csv")


if __name__ == "__main__":
    main()
