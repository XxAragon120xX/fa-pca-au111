"""Análisis de una corrida FRa–PCA. Se ejecuta dentro de <corrida>/an con ref_nowat.pdb y md_fit.xtc.
Uso: python analisis_pca.py <corrida>"""
import os, sys, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA
from MDAnalysis.lib.distances import distance_array
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from refs_pca import LAB, PCA_PT, PCA_COO, PCA_ELEM, CRY, NAT, TOP, O4_SITE, rmsd, crystal_fit
warnings.filterwarnings('ignore')

run = sys.argv[1]
T0 = 5000.0   # ps descartados
u = mda.Universe('ref_nowat.pdb', 'md_fit.xtc')
cry = mda.Universe(os.path.expanduser('~/MD/FRa_FA/input/rec_A.pdb'))
cmap = {a.resid: a.position for a in cry.select_atoms('protein and name CA')}

lig = u.select_atoms('resname PCA')
heavy = lig.select_atoms('not name H*')
elem = ''.join(a.name[0] for a in lig[:15])
assert len(lig) == 19 and len(heavy) == 15 and elem == PCA_ELEM, f'PCA inesperado: {len(lig)} átomos, {elem}'
pter, coo = lig[PCA_PT], lig[PCA_COO]
prot = u.select_atoms('protein and not name H*')
ca, keep, ref, rc = crystal_fit(u, cmap)

print(f'== {run}   ({len(u.trajectory)} frames, {u.trajectory[-1].time/1000:.1f} ns)')
print('átomos del pterín (topología -> cristal):', ' '.join(f'{a.name}={l}' for a, l in zip(pter, LAB)))
print('carboxilato:', ' '.join(a.name for a in coo))

# backbone RMSD respecto al frame 0
R = rms.RMSD(u, select='backbone').run()
t_all = R.results.rmsd[:, 1]
bb = R.results.rmsd[:, 2]

rows, fitq, Xmean, occ_p, occ_c = [], [], [], {}, {}
nfr = 0
resids = prot.resids
for ts in u.trajectory:
    mob = ca.positions[keep]; mc = mob.mean(0)
    Rm, q = align.rotation_matrix(mob - mc, ref - rc)
    X = (lig.positions - mc) @ Rm.T + rc
    Xp = X[PCA_PT]
    dmin = distance_array(heavy.positions, prot.positions, box=ts.dimensions).min()
    rows.append([ts.time, rmsd(Xp, CRY), rmsd(Xp, NAT), rmsd(Xp, TOP),
                 np.linalg.norm(Xp.mean(0) - CRY.mean(0)),
                 np.linalg.norm(X[PCA_PT[LAB.index('O4')]] - O4_SITE),
                 min(np.linalg.norm(X[i] - O4_SITE) for i in PCA_COO[1:]), dmin])
    fitq.append(q)
    if ts.time >= T0:
        nfr += 1; Xmean.append(Xp)
        for grp, occ in ((pter, occ_p), (coo, occ_c)):
            d = distance_array(grp.positions, prot.positions, box=ts.dimensions)
            for r in set(resids[(d < 4.0).any(0)]):
                occ[r] = occ.get(r, 0) + 1

A = np.array(rows)
np.savetxt(f'serie_{run}.csv', np.column_stack([A, bb]), delimiter=',', fmt='%.3f',
           header='t_ps,rmsd_cristal,rmsd_nat,rmsd_top,desplaz_centroide,dO4_O4sitio,dCOO_O4sitio,dmin_prot,rmsd_backbone', comments='')
np.save(f'pterin_media_{run}.npy', np.mean(Xmean, 0))
m = A[:, 0] >= T0
ms = lambda v: f'{v[m].mean():.2f} ± {v[m].std():.2f}'

print(f'ajuste Cα al cristal: {np.mean(fitq):.2f} Å (debe ser ~1.3–1.8)')
print(f'\nt >= 5 ns ({m.sum()} frames)')
print(f'  RMSD backbone            {ms(bb)} Å')
print(f'  pterín vs cristal        {ms(A[:,1])} Å')
print(f'  pterín vs pose nativa    {ms(A[:,2])} Å')
print(f'  pterín vs pose top       {ms(A[:,3])} Å')
print(f'  desplazamiento centroide {ms(A[:,4])} Å  (máx {A[m,4].max():.1f})')
print(f'  distancia mínima PCA–proteína: media {A[m,7].mean():.2f}, máx {A[m,7].max():.2f} Å')
fr = lambda c: f'{100*c[m].mean():5.1f} %'
print(f'  frames tipo nativa (RMSD nat < 2 Å):    {fr(A[:,2] < 2)}')
print(f'  frames tipo top    (RMSD top < 2 Å):    {fr(A[:,3] < 2)}')
print(f'  sitio del O4 ocupado por O4 de PCA (<1.5 Å):        {fr(A[:,5] < 1.5)}')
print(f'  sitio del O4 ocupado por el carboxilato (<1.5 Å):   {fr(A[:,6] < 1.5)}')
print(f'  PCA fuera del bolsillo (centroide > 8 Å):           {fr(A[:,4] > 8)}')
print(f'  PCA sin contacto con la proteína (dmin > 4 Å):      {fr(A[:,7] > 4)}')

print('\nserie (cada 5 ns): t  RMSDcris  RMSDnat  RMSDtop  despl  dO4  dCOO  dmin')
for row in A:
    if abs(row[0] % 5000) < 1e-3 or row[0] in (500.0, 1000.0):
        print(f'  {row[0]/1000:5.1f} ns  {row[1]:5.2f}  {row[2]:5.2f}  {row[3]:5.2f}  {row[4]:5.2f}  {row[5]:5.2f}  {row[6]:5.2f}  {row[7]:5.2f}')

for nombre, occ in (('pterín', occ_p), ('carboxilato', occ_c)):
    print(f'\ncontactos {nombre} (<4.0 Å, ≥10 %, numeración PDB):')
    for r in sorted(occ):
        f = 100 * occ[r] / nfr
        if f >= 10:
            rn = u.select_atoms(f'protein and resid {r}')[0].resname
            print(f'  {rn}{r+7:<4} {f:5.1f} %')

# puentes de H proteína <-> PCA (t >= 5 ns)
lab = {a.ix: l for a, l in zip(pter, LAB)}
lab.update({coo[1].ix: 'OC1', coo[2].ix: 'OC2'})
h = HBA(u, between=['protein', 'resname PCA'],
        donors_sel='(protein or resname PCA) and (name N* O*)',
        hydrogens_sel='(protein or resname PCA) and name H*',
        acceptors_sel='(protein or resname PCA) and (name N* O*)',
        d_h_cutoff=1.2, d_a_cutoff=3.5, d_h_a_angle_cutoff=150)
start = int(np.argmax(A[:, 0] >= T0))
h.run(start=start)
hb = h.results.hbonds
def name(ix):
    a = u.atoms[ix]
    return f'PCA:{lab.get(ix, a.name)}' if a.resname == 'PCA' else f'{a.resname}{a.resid+7}:{a.name}'
cnt = {}
for fr_, d, _, acc, *_ in hb:
    k = (name(int(d)), name(int(acc)))
    cnt.setdefault(k, set()).add(int(fr_))
print(f'\npuentes de H (d < 3.5 Å, ángulo > 150°, t >= 5 ns), ≥ 5 %:')
for (dn, ac), s in sorted(cnt.items(), key=lambda kv: -len(kv[1])):
    f = 100 * len(s) / m.sum()
    if f >= 5:
        print(f'  {dn:>18} -> {ac:<18} {f:5.1f} %')
print(f'\nlisto: {run}')
