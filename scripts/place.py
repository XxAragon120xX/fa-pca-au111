"""Coloca el ligando sobre la cara superior del slab.
Uso: python place.py slab.gro lig.mol2 lig.json MODO SEMILLA salida.gro
MODO: sol  -> orientación aleatoria, átomo más bajo a 12 Å de la capa superior de Au
      flat -> plano del pterín paralelo a la superficie, a 3.4 Å; el resto de la molécula hacia arriba
      tail -> FA de pie: γ-carboxilato abajo (O más bajo a 3.0 Å), pterín arriba
      head -> de pie sobre el grupo amino del pterín (amino abajo a 3.0 Å)"""
import sys, json
import numpy as np

slab, mol2, js, mode, seed, out = sys.argv[1:7]
seed = int(seed)
info = json.load(open(js))
rng = np.random.default_rng(seed)

G = open(slab).read().split('\n')
n_au = int(G[1]); au = G[2:2 + n_au]; box = np.array([float(v) for v in G[2 + n_au].split()])
zau = np.array([float(l[36:44]) for l in au]); ztop = zau.max() * 10

L = open(mol2).read().split('\n')
i0 = L.index('@<TRIPOS>ATOM'); i1 = next(i for i in range(i0 + 1, len(L)) if L[i].startswith('@<TRIPOS>'))
A = [l.split() for l in L[i0 + 1:i1]]
names = [a[1] for a in A]; resn = A[0][7][:5]
X = np.array([[float(a[2]), float(a[3]), float(a[4])] for a in A])
heavy = info['heavy']; pt = info['pterin']

def rot_to(u, v):
    """matriz que lleva el vector u a v"""
    u = u / np.linalg.norm(u); v = v / np.linalg.norm(v)
    c = np.dot(u, v); w = np.cross(u, v); s = np.linalg.norm(w)
    if s < 1e-8:
        return np.eye(3) if c > 0 else np.diag([1, -1, -1])
    K = np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]]) / s
    return np.eye(3) + s * K + (1 - c) * K @ K

def rand_rot():
    q = rng.normal(size=4); q /= np.linalg.norm(q); a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])

X = X - X[heavy].mean(0)
if mode == 'sol':
    X = X @ rand_rot().T
    X[:, 2] += ztop + 12.0 - X[heavy, 2].min()
elif mode == 'flat':
    P = X[pt] - X[pt].mean(0)
    nrm = np.linalg.svd(P)[2][2]
    X = X @ rot_to(nrm, np.array([0, 0, 1.0])).T
    rest = [i for i in heavy if i not in pt]
    if rest and X[rest, 2].mean() < X[pt, 2].mean():
        X = X @ np.diag([1, -1, -1]).T          # el resto de la molécula hacia arriba
    X = X @ rot_to(np.array([1, 0, 0]), np.array([np.cos(rng.uniform(0, 2*np.pi)), np.sin(rng.uniform(0, 2*np.pi)), 0])).T
    X[:, 2] += ztop + 3.4 - X[pt, 2].mean()
    low = X[heavy, 2].min()
    if low < ztop + 2.8:
        X[:, 2] += ztop + 2.8 - low
elif mode in ('tail', 'head'):
    if mode == 'tail':
        anchor = info['gamma_COO'] or info['ring_COO']
        assert anchor, 'el ligando no tiene carboxilato para anclar'
        base = X[anchor].mean(0); low_idx = anchor[1:]
        v = X[pt].mean(0) - base                     # pterín lo más lejos posible de la superficie
    else:
        base = X[info['amino_N']]; low_idx = [info['amino_N']]
        v = X[heavy].mean(0) - base                  # el resto de la molécula hacia arriba
    X = X @ rot_to(v, np.array([0, 0, 1.0])).T
    X[:, 2] += ztop + 3.0 - X[low_idx, 2].min()
    low = X[heavy, 2].min()
    if low < ztop + 2.8:
        X[:, 2] += ztop + 2.8 - low
else:
    raise SystemExit(f'modo desconocido: {mode}')

# centrar en x, y (con un desplazamiento aleatorio pequeño)
X[:, :2] += box[:2] * 10 / 2 - X[heavy, :2].mean(0) + rng.uniform(-5, 5, 2)
span = X.max(0) - X.min(0)
assert span[0] < box[0] * 10 - 15 and span[1] < box[1] * 10 - 15, 'el ligando no cabe lateralmente'
assert X[:, 2].max() < box[2] * 10 - 15, 'el ligando queda demasiado cerca del borde superior'

with open(out, 'w') as f:
    f.write(f'Au(111) + {resn} ({mode}, semilla {seed})\n{n_au + len(X)}\n')
    f.write('\n'.join(au) + '\n')
    for k, (nm, (x, y, z)) in enumerate(zip(names, X / 10), n_au + 1):
        f.write(f'{(n_au + 1) % 100000:5d}{resn:<5}{nm:>5}{k % 100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n')
    f.write(f'{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n')
pc = X[pt].mean(0)[2] - ztop
print(f'{out}: {mode}, altura del centroide del pterín {pc:.1f} Å, átomo pesado más bajo {X[heavy,2].min()-ztop:.1f} Å')
