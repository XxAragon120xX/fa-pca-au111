# DFT inputs

One directory per adsorption configuration. Each holds the three single-point
calculations that enter the frozen-fragment adsorption energy:

```
ΔE_ads = E(A) − E(B) − E(C)
```

| File | Calculation |
| --- | --- |
| `A.in` | ligand adsorbed on the slab (the complex) |
| `B.in` | clean slab, frozen in the adsorbed geometry |
| `C.in` | isolated ligand, frozen in the adsorbed geometry, in the same cell |

The four FA directories of the first series use the older names of the same three
calculations: `*_SCF.in` for the complex, `Slab.in` for the slab and `Mol.in` for the
ligand. The `FAfix_*` directories are the re-run of the FA systems with the corrected
slab described in the README.

## Pseudopotentials

The inputs expect a `pseudo/` directory next to them. These are the five files used,
all from PSlibrary 1.0.0, as distributed with Quantum ESPRESSO:

```
Au.pbe-n-rrkjus_psl.1.0.0.UPF
C.pbe-n-rrkjus_psl.1.0.0.UPF
H.pbe-rrkjus_psl.1.0.0.UPF
N.pbe-n-rrkjus_psl.1.0.0.UPF
O.pbe-n-rrkjus_psl.1.0.0.UPF
```

They are not redistributed here. Download them from
<https://pseudopotentials.quantum-espresso.org/legacy_tables> or from the `pslibrary`
distribution, place them in `DFT/pseudo/`, and the inputs run unchanged.

## Settings common to every calculation

| | |
| --- | --- |
| Functional | PBE |
| Dispersion | D3 with Becke–Johnson damping (`vdw_corr='grimme-d3'`, `dftd3_version=4`) |
| Cutoffs | 60 Ry (wavefunction) / 600 Ry (charge density) |
| k-mesh | 4 × 4 × 1 |
| Cell | a = b = 14.4186 Å (5 × 5 surface cell), c = 47.0636 Å |
| Slab | four Au(111) layers, 7.4 Å thick; the two bottom layers fixed |
| Vacuum | ≈ 40 Å along the surface normal |
| Smearing | Marzari–Vanderbilt, degauss = 0.01 Ry |

The dispersion term reported in Table S3 is the `DFT-D3 Dispersion` line of each
output, combined in the same way as the total energies.
