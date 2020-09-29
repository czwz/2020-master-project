#LOAD REUIQRED MODULE
module purge
module load gcc mvapich2 gromacs/2019.2-mpi

#SET UP INPUT AND OUTPUT NAME; ONLY SUPPORT PDB
#THE INPUT FILE SHOULD BE ALREADY DEPROTONATED
FILENAME=$1
INPUT=$FILENAME    # NAME-OF-INPUT.PDB
OUTPUT=$FILENAME   # NAME-OF-OUTPUT.PDB

PDB_DIR='pdb'	         # NAME-OF-PDB-FILE-DIRECTORY
OUTPUT_DIR=$FILENAME     # NAME-OF-OUTPUT-DIR
OUTLOG="${FILENAME}/log" # NAME-OF-OUTPUT-LOG

if [ -d $OUTPUT_DIR ]; then
    rm -rf $OUTPUT_DIR
    mkdir ${OUTPUT_DIR}
else
    mkdir ${OUTPUT_DIR}
fi

##GET THE FORCE FIELD FROM INTERNET
##UNCOMMENT IT IF THERE IS NO FORCE-FIELD DIRECTORY CALLED "charmm36-mar2019"
#wget http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-mar2019.ff.tgz
#mv download.php?filename=CHARMM_ff_params_files%2Fcharmm36-mar2019.ff.tgz charmm36-mar2019.ff.tgz
#tar -xvf charmm36-mar2019.ff.tgz
#rm charmm36-mar2019.ff.tgz
cp -r charmm36-mar2019.ff -t $OUTPUT_DIR

cat > ions.mdp << EOF
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform
; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.2       ; Short-range electrostatic cut-off
rvdw            = 1.2       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
EOF

#START TO PROCESS FILES
#SOLVATE SYSTEM
echo "yes\nyes" | gmx_mpi pdb2gmx -f ${PDB_DIR}/${INPUT}.pdb -o ${OUTPUT}.gro -ff charmm36-mar2019 -water tip3p -ss yes -ignh yes >> $OUTLOG |& >/dev/null
gmx_mpi editconf -f ${OUTPUT}.gro -o ${OUTPUT}.gro -c -d 1.2 -bt cubic >> $OUTLOG |& >/dev/null #-d: sets the cubic edge to be (protein+12) Ang -c: center
gmx_mpi solvate -cp ${OUTPUT}.gro -cs spc216.gro -o ${OUTPUT}.gro -p topol.top >> $OUTLOG &>/dev/null #-cp/cs: solute/solvent -p:input topology
gmx_mpi grompp -f ions.mdp -c ${OUTPUT}.gro -p topol.top -o ions.tpr >> $OUTLOG |& >/dev/null #-c: structure file
echo "13" |  gmx_mpi genion -s ions.tpr -o ${OUTPUT}.gro -p topol.top -pname K -nname CL -neutral >> $OUTLOG |& >/dev/null
gmx_mpi editconf -f ${OUTPUT}.gro -o step3_charmm2gmx.pdb >> $OUTLOG |& >/dev/null
echo "q" |   gmx_mpi make_ndx -f ${OUTPUT}.gro -o index.ndx >> $OUTLOG |& >/dev/null
rm \#* mdout.mdp


#CREATE INPUT FILE FOR MINIMIZATION, EQUILIBRATION AND PRODUCTION
cat > step4.0_minimization.mdp << EOF
define                  = -DREST_ON -DSTEP4_0 ;No idea about this option
integrator              = steep
emtol                   = 1000.0
nsteps                  = 5000
nstlist                 = 10
cutoff-scheme           = Verlet
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = pme
rcoulomb                = 1.2
;
constraints             = h-bonds
constraint_algorithm    = LINCS
EOF

cat > step4.1_equilibration.mdp << EOF
define                  = -DREST_ON -DSTEP4_1
integrator              = md
dt                      = 0.001
nsteps                  = 125000
nstxtcout               = 5000
nstvout                 = 5000
nstfout                 = 5000
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = pme
rcoulomb                = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
;
tcoupl                  = Nose-Hoover
tc_grps                 = Protein Water_and_ions
tau_t                   = 1.0     1.0
ref_t                   = 303.15    303.15
;
constraints             = h-bonds
constraint_algorithm    = LINCS
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = Protein Water_and_ions
;
gen-vel                 = yes
gen-temp                = 303.15
gen-seed                = -1
;
refcoord_scaling        = com
EOF

cat > step5_production.mdp << EOF
integrator              = md
dt                      = 0.002
nsteps                  = 5000000000
; Output control
nstxout                 = 0         ; suppress bulky .trr file by specifying
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 100000      ; save energies every 10.0 ps
nstlog                  = 100000      ; update log file every 10.0 ps
nstxout-compressed      = 100000      ; save compressed coordinates every 200 ps
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
coulombtype             = pme
rcoulomb                = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 0.8
rvdw                    = 1.2
;
tcoupl                  = Nose-Hoover
tc_grps                 = Protein Water_and_ions
tau_t                   = 1.0     1.0
ref_t                   = 303.15 303.15
;
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 5.0
compressibility         = 4.5e-5
ref_p                   = 1.0
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = Protein Water_and_ions
;
refcoord_scaling        = com
EOF

mv step* topol.top posre.itp ions* index.ndx *.gro -t $OUTPUT_DIR
