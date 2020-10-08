module load daint-gpu
module load GROMACS

echo "0" | gmx trjconv -s step4.1_equilibration.tpr -f step4.1_equilibration.xtc -o step4.1_equilibration-noPBC.xtc -pbc mol -ur compact
echo "0" | gmx trjconv -s step5_1.tpr -f step5_1.xtc -o step5_1-noPBC.xtc -pbc mol -ur compact
echo "0" | gmx trjconv -s step5_1.tpr -f step5_1.part0002.xtc -o step5_1.part0002-noPBC.xtc -pbc mol -ur compact
echo "0" | gmx trjconv -s step5_1.tpr -f step5_1.part0003.xtc -o step5_1.part0003-noPBC.xtc -pbc mol -ur compact

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

#################################
#clip till max duration and save#
#################################

animate delete beg 5000 end 10000 skip -1
save_state view_state.vmd

############################
#alignment CA w.r.t CA(t=0)#
############################

set ref [atomselect top "protein and name CA" frame 0]
set sel [atomselect top "protein and name CA"]
set all [atomselect top "all"]
set n [molinfo top get numframes]

for {set i 1} {\$i < \$n} {incr i} {
    \$sel frame \$i   
    \$all frame \$i
    \$all move [measure fit \$sel \$ref]
    }

##########################################################
#calculate & output protein CA rmsd aligned w.r.t CA(t=0)#
##########################################################

set outfile [open "$NAME.prmsd" w]
set ref [atomselect top "protein and name CA" frame 0]
set sel [atomselect top "protein and name CA"]
set n [molinfo top get numframes]
for {set i 1} {\$i < \$n} {incr i} {
  \$sel frame \$i
  set rmsd [measure rmsd \$sel \$ref]
  puts \$outfile "\$i \$rmsd"
  }
close \$outfile

##########################################################
#calculate & output epitope CA rmsd aligned w.r.t CA(t=0)#
##########################################################

set outfile [open "$NAME.ermsd" w]
set ref [atomselect top "protein and name CA and (resid $epitope)" frame 0]
set sel [atomselect top "protein and name CA and (resid $epitope)"]
for {set i 0} {\$i < \$n} {incr i} {
  \$sel frame \$i
  set rmsd [measure rmsd \$sel \$ref]
  puts \$outfile "\$i \$rmsd"
}
close \$outfile

##################################################
#calculate & output CA rmsf aligned w.r.t CA(t=0)#
##################################################

set outfile [open "$NAME.rmsf" w]
set sel [atomselect top "protein and name CA"]
set mol [\$sel molindex]
for {set i 0} {\$i < [\$sel num]} {incr i} {
    set rmsf [measure rmsf \$sel]
    puts \$outfile "[expr {\$i+1}] [lindex \$rmsf \$i]"
    }
close \$outfile

####################################
#calculate & output rgyr of protein#
####################################

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

##################################################################
#calculate & output CA native contact number series w/ cutoff 4.5#
##################################################################

proc calc_nc {old_list new_list} {
    set nc 0
    foreach elem \$old_list {
        if {\$elem in \$new_list} {
            set nc [expr \$nc+1]
        }
    }
    return \$nc
}

set outfile [open "$NAME.nc" w]
set nf [molinfo top get numframes]
set ca [atomselect top "name CA"]

for {set f 0} {\$f < \$nf} {incr f 25} {
    set ncf 0
    foreach index [\$ca get index] {
        set neighbor [atomselect top "(all within 4.5 of index \$index) and ({name \"C.*\"} or {name \"N.*\"} or {name \"O.*\"} and not {name \"OH.*\"})"]

        \$neighbor frame 0
        \$neighbor update
        set init_neighbor [\$neighbor get index]

        \$neighbor frame \$f
        \$neighbor update
        set current_neighbor [\$neighbor get index]

        set nci [calc_nc \$init_neighbor \$current_neighbor]
        set ncf [expr \$ncf+\$nci]
	if {\$f==0} {set ncf0 [expr \$ncf+0.]}
    }
    puts \$outfile  "\$f [expr \$ncf/\$ncf0]"
}
close \$outfile

exit
EOF

eval "$VMD_EXE -e vmd.tcl";
