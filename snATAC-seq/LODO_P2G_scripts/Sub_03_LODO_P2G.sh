cd /public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/LODO_P2G/LODO_P2G_scripts
qsub -l nodes=1:ppn=8,vmem=190gb -m ae -o 03_LODO_P2G.o -e 03_LODO_P2G.e -N 03_LODO_P2G 03_LODO_P2G.sh