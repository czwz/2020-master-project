module purge
module load gcc/7.4.0
module load mvapich2/2.3.1
module load gromacs/2019.2-mpi

echo "0" | gmx_mpi trjconv -s step4.1_equilibration.tpr -f step4.1_equilibration.xtc -o step4.1_equilibration-noPBC.xtc -pbc mol -ur compact
echo "0" | gmx_mpi trjconv -s step5_1.tpr -f step5_1.xtc -o step5_1-noPBC.xtc -pbc mol -ur compact
echo "0" | gmx_mpi trjconv -s step5_1.tpr -f step5_1.part0002.xtc -o step5_1.part0002-noPBC.xtc -pbc mol -ur compact

rsync -t -v -v -P *PBC.xtc slwang@lbmpc4:/data/slwang/MD_input_structures/3E2H/3E2H_enrich_best/charmm-gui-8374552026/gromacs
