"""Escribe un mol2 por pose: cargas y tipos de PCA_q.mol2 (antechamber), coordenadas del SDF del docking.
Uso: python mol2_pose.py PCA_q.mol2 pose.sdf salida.mol2
Ambos archivos deben tener los mismos átomos en el mismo orden (antechamber conserva el orden del SDF)."""
import sys

q, sdf, out = sys.argv[1:4]
L = open(q).read().split('\n')
i0 = L.index('@<TRIPOS>ATOM')
i1 = next(i for i in range(i0 + 1, len(L)) if L[i].startswith('@<TRIPOS>'))
atoms = L[i0 + 1:i1]

s = open(sdf).read().split('\n')
n = int(s[3][:3])
xyz = [(float(l[0:10]), float(l[10:20]), float(l[20:30]), l[31:34].strip()) for l in s[4:4 + n]]
assert n == len(atoms), f'átomos: SDF {n} vs mol2 {len(atoms)}'

new = []
for k, (a, (x, y, z, el)) in enumerate(zip(atoms, xyz)):
    c = a.split()
    name, typ = c[1], c[5]
    # comprobar elemento: primera letra del nombre (antechamber) vs elemento del SDF
    assert name[0].upper() == el[0].upper(), f'átomo {k+1}: {name} vs {el}'
    new.append(f'{int(c[0]):>7} {name:<8}{x:>10.4f}{y:>10.4f}{z:>10.4f} {typ:<6}{int(c[6]):>4} {c[7]:<8}{float(c[8]):>10.6f}')

open(out, 'w').write('\n'.join(L[:i0 + 1] + new + L[i1:]))
qt = sum(float(a.split()[8]) for a in atoms)
print(f'{out}: {n} átomos, carga total {qt:+.3f}')
