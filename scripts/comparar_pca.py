"""Compara la posición media del pterín de cada corrida de PCA con el cristal, las poses del docking
y la posición media del pterín de FA en sus seis corridas (todas en el marco del cristal).
Uso (desde ~/MD/FRa_PCA, entorno mda): python scripts/comparar_pca.py"""
import os, sys, glob, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from refs_pca import FOA_PT, CRY, NAT, TOP, rmsd, crystal_fit
warnings.filterwarnings('ignore')

B = os.path.expanduser('~/MD/FRa_PCA'); FA = os.path.expanduser('~/MD/FRa_FA')
cache = f'{B}/fa_pterin_medias.npz'
if os.path.exists(cache):
    fa = dict(np.load(cache))
else:
    cry = mda.Universe(f'{FA}/input/rec_A.pdb')
    cmap = {a.resid: a.position for a in cry.select_atoms('protein and name CA')}
    fa = {}
    for r in ['rep1', 'rep2', 'cry1', 'hie1', 'hie2', 'hie3']:
        u = mda.Universe(f'{FA}/{r}/an/ref_nowat.pdb', f'{FA}/{r}/an/md_fit.xtc')
        lig = u.select_atoms('resname FOA')
        pt = sum((lig.select_atoms(f'name {n}') for n in FOA_PT[1:]), lig.select_atoms(f'name {FOA_PT[0]}'))
        assert len(pt) == 12
        ca, keep, ref, rc = crystal_fit(u, cmap); acc = []
        for ts in u.trajectory[250:]:
            mob = ca.positions[keep]; mc = mob.mean(0)
            R, _ = align.rotation_matrix(mob - mc, ref - rc)
            acc.append((pt.positions - mc) @ R.T + rc)
        fa[r] = np.mean(acc, 0); print(f'FA {r}: listo', file=sys.stderr)
    np.savez(cache, **fa)
fa['FA_HIE_media'] = np.mean([fa[k] for k in ('hie1', 'hie2', 'hie3')], 0)

refs = {'cristal': CRY, 'pose_nat': NAT, 'pose_top': TOP, 'FA_HIE_media': fa['FA_HIE_media'],
        'FA_rep1': fa['rep1'], 'FA_cry1': fa['cry1']}
runs = sorted(glob.glob(f'{B}/pca_*/an/pterin_media_*.npy'))
print('RMSD (Å) entre la posición media del pterín de cada corrida de PCA (t >= 5 ns) y las referencias')
print(f'{"":10s}' + ''.join(f'{k:>14s}' for k in refs))
print(f'{"FA_HIE":10s}' + ''.join(f'{rmsd(fa["FA_HIE_media"], v):14.2f}' for v in refs.values()))
P = {}
for f in runs:
    n = os.path.basename(f)[len('pterin_media_'):-4]; P[n] = np.load(f)
    print(f'{n:10s}' + ''.join(f'{rmsd(P[n], v):14.2f}' for v in refs.values()))
if len(P) > 1:
    print('\nentre corridas de PCA:')
    ks = list(P)
    for i, a in enumerate(ks):
        for b in ks[i + 1:]:
            print(f'  {a} vs {b}: {rmsd(P[a], P[b]):.2f}')
