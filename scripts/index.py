"""Escribe index.ndx con los grupos System, AU, LIG, nonAU, WAT e Iones a partir de un .gro.
Uso: python index.py sys.gro LIG index.ndx"""
import sys
gro, lig, out = sys.argv[1:4]
G = open(gro).read().split('\n'); n = int(G[1])
groups = {'System': [], 'AU': [], 'LIG': [], 'nonAU': [], 'WAT': [], 'Iones': []}
for k, l in enumerate(G[2:2 + n], 1):
    rn = l[5:10].strip()
    groups['System'].append(k)
    if rn == 'AU':
        groups['AU'].append(k)
    else:
        groups['nonAU'].append(k)
        if rn == lig: groups['LIG'].append(k)
        elif rn == 'WAT': groups['WAT'].append(k)
        else: groups['Iones'].append(k)
with open(out, 'w') as f:
    for g, ids in groups.items():
        if not ids: continue
        f.write(f'[ {g} ]\n')
        for i in range(0, len(ids), 15):
            f.write(' '.join(f'{x:6d}' for x in ids[i:i + 15]) + '\n')
print(f'{out}: ' + ', '.join(f'{g} {len(v)}' for g, v in groups.items()))
