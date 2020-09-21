#LOAD REUIQRED MODULE
module purge
module load gcc mvapich2 gromacs/2019.2-mpi

#SET UP INPUT AND OUTPUT NAME; ONLY SUPPORT PDB FILE
#THE INPUT FILE SHOULD BE ALREADY DEPROTONATED
INPUT='dH_3E2H_enrich_best' #NO NEED TO ADD .PDB
OUTPUT='3E2H_enrich_best'   #NO NEED TO ADD .PDB

#GET THE FORCE FIELD FROM INTERNET
wget http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-mar2019.ff.tgz
mv download.php?filename=CHARMM_ff_params_files%2Fcharmm36-mar2019.ff.tgz charmm36-mar2019.ff.tgz
tar -xvf charmm36-mar2019.ff.tgz
rm charmm36-mar2019.ff.tgz

#CREATE INPUT FILE FOR PREPARATION
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
echo "yes" | gmx_mpi pdb2gmx -f ${INPUT}.pdb -o ${OUTPUT}.gro -ff charmm36-mar2019 -water tip3p -ss yes
gmx_mpi editconf -f ${OUTPUT}.gro -o ${OUTPUT}.gro -c -d 0.5 -bt cubic
gmx_mpi solvate -cp ${OUTPUT}.gro -cs spc216.gro -o ${OUTPUT}.gro -p topol.top
gmx_mpi grompp -f ions.mdp -c ${OUTPUT}.gro -p topol.top -o ions.tpr
echo "13" | gmx_mpi genion -s ions.tpr -o ${OUTPUT}.gro -p topol.top -pname K -nname CL -neutral
gmx_mpi editconf -f ${OUTPUT}.gro -o step3_charmm2gmx.pdb
echo "q" | gmx_mpi make_ndx -f ${OUTPUT}.gro -o index.ndx
rm \#* mdout.mdp posre.itp ions.*


#CREATE INPUT FILE FOR MINIMIZATION, EQUILIBRATION AND PRODUCTION
cat > step4.0_minimization.mdp << EOF
define                  = -DREST_ON -DSTEP4_0
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
