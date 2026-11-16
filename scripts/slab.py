"""Slab de Au(111) en celda ortorrómbica (apilamiento ABC), a = 4.078 Å.
Uso: python slab.py NX NY NCAPAS BOXZ_nm salida.gro
Celda: x = NX·d, y = NY·√3·d, con d = a/√2. Capa inferior en z = 0.5 nm."""
import sys
import numpy as np

nx, ny, nl = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
boxz, out = float(sys.argv[4]), sys.argv[5]
a = 4.078
d = a / np.sqrt(2)                 # distancia Au–Au en el plano
dz = a / np.sqrt(3)                # distancia entre capas (111)
Lx, Ly = nx * d, ny * np.sqrt(3) * d
X = []
for k in range(nl):
    sx, sy = (k % 3) * d / 2, (k % 3) * d * np.sqrt(3) / 6
    for j in range(2 * ny):
        for i in range(nx):
            X.append([(i * d + (j % 2) * d / 2 + sx) % Lx, (j * d * np.sqrt(3) / 2 + sy) % Ly, 5.0 + k * dz])
X = np.array(X) / 10.0
with open(out, 'w') as f:
    f.write(f'Au(111) {nx}x{2*ny} rows x {nl} layers, a=4.078 A\n{len(X)}\n')
    for n, (x, y, z) in enumerate(X, 1):
        f.write(f'{n % 100000:5d}{"AU":<5}{"AU":>5}{n % 100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n')
    f.write(f'{Lx/10:10.5f}{Ly/10:10.5f}{boxz:10.5f}\n')
print(f'{out}: {len(X)} Au, caja {Lx/10:.4f} x {Ly/10:.4f} x {boxz:.2f} nm, capa superior z = {X[:,2].max():.4f} nm')
