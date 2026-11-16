"""Tras gmx solvate: elimina las aguas dentro o pegadas al slab de Au y renombra SOL -> WAT (O, H1, H2).
Uso: python clean_water.py solvatado.gro salida.gro [margen_A=2.6]"""
import sys
gin, gout = sys.argv[1], sys.argv[2]
marg = float(sys.argv[3]) / 10 if len(sys.argv) > 3 else 0.26
G = open(gin).read().split('\n')
n = int(G[1]); atoms = G[2:2 + n]; box = G[2 + n]
zau = [float(l[36:44]) for l in atoms if l[5:10].strip() == 'AU']
lo, hi = min(zau) - marg, max(zau) + marg
keep, i, nw, drop = [], 0, 0, 0
while i < n:
    l = atoms[i]
    if l[5:10].strip() == 'SOL':
        w = atoms[i:i + 3]
        if any(lo <= float(a[36:44]) <= hi for a in w):
            drop += 1
        else:
            for a, nm in zip(w, ('O', 'H1', 'H2')):
                keep.append((a[:5], 'WAT', nm, a[20:44])); nw += 1 if nm == 'O' else 0
        i += 3
    else:
        keep.append((l[:5], l[5:10].strip(), l[10:15].strip(), l[20:44])); i += 1
# renumerar residuos de agua y átomos
out, res, last = [], 0, None
with open(gout, 'w') as f:
    f.write(G[0] + '\n' + f'{len(keep)}\n')
    rcount = 0
    for k, (rid, rn, an, xyz) in enumerate(keep, 1):
        if rn == 'WAT':
            if an == 'O':
                rcount += 1
            rnum = rcount
        else:
            rnum = int(rid)
        f.write(f'{rnum % 100000:5d}{rn:<5}{an:>5}{k % 100000:5d}{xyz}\n')
    f.write(box + '\n')
print(f'{gout}: {nw} aguas conservadas, {drop} eliminadas cerca del Au')
