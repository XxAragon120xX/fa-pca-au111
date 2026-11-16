"""Análisis de una corrida de FA o PCA sobre Au(111) en agua.
Uso (entorno mda), desde ~/MD/Au_MD:   python scripts/analisis_au.py fa3_sol1 [fa3_sol2 ...]
Lee <corrida>/npt.gro (posiciones fijas del Au y orden de átomos), <corrida>/md.xtc (solo el ligando)
y prep/<LIG>/<LIG>.json (fragmentos). Las alturas se miden respecto a la cara del slab más cercana
(capa superior o inferior: con agua a ambos lados, el ligando puede adsorberse en cualquiera).
Salida en <corrida>/an_au/: resumen.txt y serie.csv"""
import os, sys, json, warnings
import numpy as np
import MDAnalysis as mda
warnings.filterwarnings('ignore')

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H_CONTACT = 3.5     # Å: átomo pesado "en contacto" con el oro
H_RECOG = 15.0      # Å: altura mínima del pterín compatible con el reconocimiento por FRa (Sección 3.4.3)

def gro_atoms(fn):
    G = open(fn).read().split('\n'); n = int(G[1])
    A = [(l[5:10].strip(), l[10:15].strip(), float(l[20:28]) * 10, float(l[28:36]) * 10, float(l[36:44]) * 10) for l in G[2:2 + n]]
    box = np.array([float(v) for v in G[2 + n].split()[:3]]) * 10
    return A, box

def analizar(run):
    D = os.path.join(B, run); out = os.path.join(D, 'an_au'); os.makedirs(out, exist_ok=True)
    A, box = gro_atoms(os.path.join(D, 'npt.gro' if os.path.exists(os.path.join(D, 'npt.gro')) else 'sys.gro'))
    zau = np.array([a[4] for a in A if a[0] == 'AU'])
    zs = np.sort(zau); cuts = np.where(np.diff(zs) > 1.0)[0]
    capas = np.split(zs, cuts + 1)                                   # capas de Au separadas por > 1 Å
    layers = capas; ztop, zbot = capas[-1].mean(), capas[0].mean()   # plano medio de la capa superior / inferior
    lig = [a for a in A if a[0] not in ('AU', 'WAT', 'Na+', 'Cl-')]
    resn = lig[0][0]; info = json.load(open(os.path.join(B, 'prep', resn, f'{resn}.json')))
    # topología del ligando para leer el xtc (que solo contiene el ligando)
    pdb = os.path.join(out, 'lig.pdb')
    with open(pdb, 'w') as f:
        for i, (rn, an, x, y, z) in enumerate(lig, 1):
            f.write(f'HETATM{i:5d} {an:<4} {rn:<3} A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n')
        f.write('END\n')
    u = mda.Universe(pdb, os.path.join(D, 'md.xtc'))
    heavy = np.array(info['heavy']); pt = np.array(info['pterin'])
    frag = {'pterina': pt}
    if info['pABA_ring']: frag['anillo_pABA'] = np.array(info['pABA_ring'])
    if info['gamma_COO']: frag['gamma_COO'] = np.array(info['gamma_COO'])
    if info['alpha_COO']: frag['alpha_COO'] = np.array(info['alpha_COO'])
    if info['ring_COO']: frag['COO'] = np.array(info['ring_COO'])
    T, Hfr, Hmin, ANG, FACE = [], {k: [] for k in frag}, [], [], []
    Cfr = {k: [] for k in frag}
    for ts in u.trajectory:
        X = u.atoms.positions.copy(); L = ts.dimensions[:3] if ts.dimensions is not None else box
        X = X[0] + (X - X[0]) - L * np.round((X - X[0]) / L)          # molécula entera
        cz = X[heavy, 2].mean()
        dtop = (cz - ztop) % L[2]; dbot = (zbot - cz) % L[2]
        if dtop <= dbot:
            h = ((X[:, 2] - ztop + L[2] / 2) % L[2]) - L[2] / 2; FACE.append('sup')
        else:
            h = ((zbot - X[:, 2] + L[2] / 2) % L[2]) - L[2] / 2; FACE.append('inf')
        T.append(ts.time / 1000); Hmin.append(h[heavy].min())
        for k, idx in frag.items():
            Hfr[k].append(h[idx].mean()); Cfr[k].append(h[idx].min() < H_CONTACT)
        P = X[pt] - X[pt].mean(0); nrm = np.linalg.svd(P)[2][2]
        ANG.append(np.degrees(np.arccos(abs(nrm[2]))))                 # 0° = pterina plana sobre el oro, 90° = de canto
    T = np.array(T); Hmin = np.array(Hmin); ANG = np.array(ANG)
    cols = [T, Hmin] + [np.array(Hfr[k]) for k in frag] + [ANG]
    np.savetxt(os.path.join(out, 'serie.csv'), np.column_stack(cols), delimiter=',', fmt='%.3f',
               header='t_ns,h_min_ligando,' + ','.join(f'h_{k}' for k in frag) + ',angulo_pterina', comments='')
    ads = Hmin < 4.0
    t_ads = T[np.argmax(ads)] if ads.any() else None
    last = T >= T[-1] - 100 if T[-1] > 150 else T >= T[-1] / 2
    hp = np.array(Hfr['pterina'])
    L = [f'== {run} ({resn}, {len(T)} frames, {T[-1]:.1f} ns; capas de Au: {len(layers)}; cara usada: '
         f'{100*np.mean(np.array(FACE)=="sup"):.0f} % superior)',
         f'primera adsorción (algún átomo pesado < 4 Å): ' + (f'{t_ads:.2f} ns' if t_ads is not None else 'nunca'),
         f'fracción adsorbida (todo el tiempo / último tramo): {100*ads.mean():.1f} % / {100*ads[last].mean():.1f} %',
         f'altura del centroide de la pterina (último tramo, t >= {T[last][0]:.0f} ns): {hp[last].mean():.2f} ± {hp[last].std():.2f} Å'
         f'   (mín {hp[last].min():.1f}, máx {hp[last].max():.1f})',
         f'fracción con la pterina a >= {H_RECOG:.0f} Å (último tramo): {100*np.mean(hp[last] >= H_RECOG):.1f} %',
         f'ángulo pterina–superficie (último tramo): {ANG[last].mean():.0f} ± {ANG[last].std():.0f}°  (0° = plana sobre el oro)',
         '', f'fragmento        altura media (Å)   en contacto con Au (< {H_CONTACT} Å), último tramo']
    for k in frag:
        a = np.array(Hfr[k]); c = np.array(Cfr[k])
        L.append(f'  {k:14s}  {a[last].mean():6.2f} ± {a[last].std():4.2f}      {100*c[last].mean():5.1f} %')
    txt = '\n'.join(L); open(os.path.join(out, 'resumen.txt'), 'w').write(txt + '\n'); print(txt + '\n')

for r in sys.argv[1:]:
    analizar(r)
