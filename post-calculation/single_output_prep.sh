module load daint-gpu
module load GROMACS

echo "0" | gmx trjconv -s step4.1_equilibration.tpr -f step4.1_equilibration.xtc -o step4.1_equilibration-noPBC.xtc -pbc mol -ur compact >> vmd.log 2>&1
echo "0" | gmx trjconv -s step5_1.tpr -f step5_1.xtc -o step5_1-noPBC.xtc -pbc mol -ur compact >> vmd.log 2>&1
echo "0" | gmx trjconv -s step5_1.tpr -f step5_1.part0002.xtc -o step5_1.part0002-noPBC.xtc -pbc mol -ur compact >> vmd.log 2>&1
echo "0" | gmx trjconv -s step5_1.tpr -f step5_1.part0003.xtc -o step5_1.part0003-noPBC.xtc -pbc mol -ur compact >> vmd.log 2>&1
echo "0" | gmx trjconv -s step5_1.tpr -f step5_1.part0004.xtc -o step5_1.part0004-noPBC.xtc -pbc mol -ur compact >> vmd.log 2>&1
echo "0" | gmx trjconv -s step5_1.tpr -f step5_1.part0005.xtc -o step5_1.part0005-noPBC.xtc -pbc mol -ur compact >> vmd.log 2>&1
echo "0" | gmx trjconv -s step5_1.tpr -f step5_1.part0006.xtc -o step5_1.part0006-noPBC.xtc -pbc mol -ur compact >> vmd.log 2>&1
echo "0" | gmx trjconv -s step5_1.tpr -f step5_1.part0007.xtc -o step5_1.part0007-noPBC.xtc -pbc mol -ur compact >> vmd.log 2>&1
echo "0" | gmx trjconv -s step5_1.tpr -f step5_1.part0008.xtc -o step5_1.part0008-noPBC.xtc -pbc mol -ur compact >> vmd.log 2>&1

NAME=$1
VMD_EXE='/users/swang/program/vmd-1.9.4a43/bin/vmd'

if [[ $NAME == *"4E1H"* ]]; then
    epitope='11 to 16'
elif [[ $NAME == *"3E2H"* ]]; then
    epitope='12 to 17'
elif [[ $NAME == *"4E2H"* ]]; then
    epitope='45 to 50'
elif [[ $NAME == *"4H"* ]]; then
    epitope='11 to 39 57 to 71'
fi

cat > vmd.tcl << EOF
################
#load trj files#
################

mol new step3_charmm2gmx.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile step4.1_equilibration-noPBC.xtc type xtc first 0 last -1 step 20 filebonds 1 autobonds 1 waitfor all
mol addfile step5_1-noPBC.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile step5_1.part0002-noPBC.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile step5_1.part0003-noPBC.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile step5_1.part0004-noPBC.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile step5_1.part0005-noPBC.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile step5_1.part0006-noPBC.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile step5_1.part0007-noPBC.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile step5_1.part0008-noPBC.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

#################################
#clip till max duration and save#
#################################

animate delete beg 20000 end 99999 skip -1
save_state view_state.vmd

############################
#alignment CA w.r.t CA(t=0)#
############################

set ref [atomselect top "protein and name CA" frame 0]
set sel [atomselect top "protein and name CA"]
set all [atomselect top "all"]
set n [molinfo top get numframes]

for {set i 0} {\$i < \$n} {incr i} {
    \$sel frame \$i   
    \$all frame \$i
    \$all move [measure fit \$sel \$ref]
    }

#################################################################
#prmsd: calculate & output overall CA rmsd aligned w.r.t CA(t=0)#
#################################################################

set outfile [open "$NAME.pprmsd" w]
set ref [atomselect top "protein and name CA" frame 0]
set sel [atomselect top "protein and name CA"]
set n [molinfo top get numframes]
for {set i 0} {\$i < \$n} {incr i} {
  \$sel frame \$i
  set rmsd [measure rmsd \$sel \$ref]
  puts \$outfile "\$i \$rmsd"
  }
close \$outfile

##############################################################################
#prmsd-side: calculate & output overall side chain rmsd aligned w.r.t CA(t=0)#
##############################################################################

set outfile [open "$NAME.sprmsd" w]
set ref [atomselect top "protein and noh and not name CA" frame 0]
set sel [atomselect top "protein and noh and not name CA"]
set n [molinfo top get numframes]
for {set i 0} {\$i < \$n} {incr i} {
  \$sel frame \$i
  set rmsd [measure rmsd \$sel \$ref]
  puts \$outfile "\$i \$rmsd"
  }
close \$outfile

##################################################################
#ermsd: calculate & output epitope CA rmsd aligned w.r.t  CA(t=0)#
##################################################################

set outfile [open "$NAME.eprmsd" w]
set ref [atomselect top "protein and name CA and (resid $epitope)" frame 0]
set sel [atomselect top "protein and name CA and (resid $epitope)"]
for {set i 0} {\$i < \$n} {incr i} {
  \$sel frame \$i
  set rmsd [measure rmsd \$sel \$ref]
  puts \$outfile "\$i \$rmsd"
  }
close \$outfile

########################################################
#rsmf: calculate & output CA rmsf aligned w.r.t CA(t=0)#
########################################################

set outfile [open "$NAME.rmsf" w]
set sel [atomselect top "protein and name CA"]
set mol [\$sel molindex]
for {set i 0} {\$i < [\$sel num]} {incr i} {
    set rmsf [measure rmsf \$sel]
    puts \$outfile "[expr {\$i+1}] [lindex \$rmsf \$i]"
    }
close \$outfile

##########################################
#rgyr: calculate & output rgyr of protein#
##########################################

set mol [molinfo top]
set outfile [open "$NAME.rgyr" w]
set sel [atomselect top "protein"]
set frames [molinfo \$mol get numframes]

for {set i 0} {\$i < \$frames} {incr i} {
    \$sel frame \$i
    \$sel update
    set a [measure rgyr \$sel]
    puts \$outfile "\$i \$a"
    }
close \$outfile

#######################################################################
#pnc: calculate & output overall native contact fraction w/ cutoff 4.5#
#######################################################################

source VMDextensions.tcl

set mol [molinfo top]
set outfile [open "$NAME.pnc" w]
set frames [molinfo \$mol get numframes]
set ca [atomselect top "protein and noh"]

\$ca frame 0
\$ca update
set ref [ prepareNativeContacts 4.5 \$ca];
set nnc [ llength \$ref]

for {set i 0} {\$i < \$frames} {incr i 10} {
    \$ca frame \$i
    \$ca update
    set nc [measureNativeContacts \$ref 4.5 \$ca]
    set qnc [ expr 1.0 * \$nc / \$nnc ]
    puts \$outfile  "\$i \$qnc"
}
close \$outfile

##########################################################################
#Recall the visualization state wo/ alignment and realign wrt. epitope CA#
##########################################################################

mol delete top
play view_state.vmd

####################################
#alignment CA w.r.t epitope CA(t=0)#
####################################

set ref [atomselect top "protein and name CA and (resid $epitope)" frame 0]
set sel [atomselect top "protein and name CA and (resid $epitope)"]
set all [atomselect top "all"]
set n [molinfo top get numframes]

for {set i 1} {\$i < \$n} {incr i} {
    \$sel frame \$i
    \$all frame \$i
    \$all move [measure fit \$sel \$ref]
    }

#########################################################################
#ermsd: calculate & output epitope CA rmsd aligned w.r.t epitope CA(t=0)#
#########################################################################

set outfile [open "$NAME.eermsd" w]
set ref [atomselect top "protein and name CA and (resid $epitope)" frame 0]
set sel [atomselect top "protein and name CA and (resid $epitope)"]
for {set i 0} {\$i < \$n} {incr i} {
  \$sel frame \$i
  set rmsd [measure rmsd \$sel \$ref]
  puts \$outfile "\$i \$rmsd"
  }
close \$outfile

######################################################################################
#ermsd-side: calculate & output epitope side chain rmsd aligned w.r.t epitope CA(t=0)#
######################################################################################

set outfile [open "$NAME.sermsd" w]
set ref [atomselect top "protein and noh and not name CA and (resid $epitope)" frame 0]
set sel [atomselect top "protein and noh and not name CA and (resid $epitope)"]
for {set i 0} {\$i < \$n} {incr i} {
  \$sel frame \$i
  set rmsd [measure rmsd \$sel \$ref]
  puts \$outfile "\$i \$rmsd"
  }
close \$outfile

#######################################################################
#enc: calculate & output epitope native contact fraction w/ cutoff 4.5#
#######################################################################

source VMDextensions.tcl

set mol [molinfo top]
set outfile [open "$NAME.enc" w]
set frames [molinfo \$mol get numframes]
set ca [atomselect top "protein and (resid $epitope) and noh"]

\$ca frame 0
\$ca update
set ref [ prepareNativeContacts 4.5 \$ca];
set nnc [ llength \$ref]

for {set i 0} {\$i < \$frames} {incr i 10} {
    \$ca frame \$i
    \$ca update
    set nc [measureNativeContacts \$ref 4.5 \$ca]
    set qnc [ expr 1.0 * \$nc / \$nnc ]
    puts \$outfile  "\$i \$qnc"
}
close \$outfile

exit
EOF

eval "$VMD_EXE -dispdev text -eofexit -e vmd.tcl >> vmd.log 2>&1";
