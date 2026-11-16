#!/bin/bash
# Uso: ./analizar_pca.sh pca_nat1 [tiempo_final_ps]
RUN=$1; END=${2:-100000}
B=~/MD/FRa_PCA
module load EESSI/2023.06 GROMACS/2024.4-foss-2023b; unset PYTHONPATH
source ~/miniconda3/etc/profile.d/conda.sh
cd $B/$RUN && mkdir -p an && cd an
TPR=../md.tpr; XTC=../md.xtc
[ -f $XTC ] || { echo "no existe $RUN/md.xtc"; exit 1; }
echo q | gmx make_ndx -f $TPR -o index0.ndx > make_ndx0.log 2>&1
NP=$(grep '^\[' index0.ndx | grep -n -x '\[ Protein \]' | head -n1 | cut -d: -f1)
NL=$(grep '^\[' index0.ndx | grep -n -x '\[ PCA \]' | head -n1 | cut -d: -f1)
[ -n "$NP" ] && [ -n "$NL" ] || { echo "no encontré los grupos Protein/PCA en index0.ndx"; exit 1; }
NP=$((NP-1)); NL=$((NL-1))
printf "%s | %s\nq\n" $NP $NL | gmx make_ndx -f $TPR -n index0.ndx -o index.ndx > make_ndx.log 2>&1
echo -e "Protein\nProtein_PCA"  | gmx trjconv -s $TPR -f $XTC -n index.ndx -o md_c.xtc -pbc mol -center -ur compact -e $END > trj1.log 2>&1
echo -e "Backbone\nProtein_PCA" | gmx trjconv -s $TPR -f md_c.xtc -n index.ndx -o md_fit.xtc -fit rot+trans > trj2.log 2>&1
echo "Protein_PCA" | gmx trjconv -s $TPR -f md_fit.xtc -n index.ndx -o ref_nowat.pdb -dump 0 > trj3.log 2>&1
conda activate mda
python $B/scripts/analisis_pca.py $RUN > resumen.txt 2>&1
tail -n 1 resumen.txt
