LIST_PROTEIN='4e2h 3e2h 4h 4e1h'
LIST_REPLICA='r1 r2 r3 r4 r5'
LIST_MODIFICATION='g b c'
LIST_NUM_EACH='1 2 3'
VMD_EXE='vmd'
PYTHON_EXE='python'

for p in $LIST_PROTEIN; do
    for m in $LIST_MODIFICATION; do
        for n in $LIST_NUM_EACH; do
	    for r in $LIST_REPLICA; do

FIL_DIR="/data/slwang/result-replica/${p}/${m}${n}/${r}"
OUT_DIR="./result"
OUTPUT="log"
INPUT="script_${p}_${m}${n}_${r}"

ls $FIL_DIR > /dev/null 2>&1

if [ $? -eq 0 ]; then

    if [ "$p" = "4e1h" ]; then
        epitope='11 to 16'
    elif [ "$p" = "4e1h" ]; then
        epitope='12 to 17'
    elif [ "$p" = "4e2h" ]; then
        epitope='45 to 50'
    elif [ "$p" = "3e2h" ]; then
        epitope='11 to 39 57 to 71'
    fi

cat > $OUT_DIR/$INPUT << EOF

##########
#read trj#
##########

mol new $FIL_DIR/step3_charmm2gmx.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile $FIL_DIR/step3_charmm2gmx.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

mol addfile $FIL_DIR/step4.1_equilibration-noPBC.xtc type xtc first 0 last -1 step 20 filebonds 1 autobonds 1 waitfor all
mol addfile $FIL_DIR/step5_1-noPBC.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile $FIL_DIR/step5_1.part0002-noPBC.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile $FIL_DIR/step5_1.part0003-noPBC.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

##################
#set max duration#
##################

animate delete beg 2500 end 10000 skip -1

#########################
#alignment w.r.t CA(t=0)#
#########################
#align all atoms w.r.t protein's alpha carbon CA at frame 0

set ref [atomselect top "protein and name CA" frame 0]
set sel [atomselect top "protein and name CA"]
set all [atomselect top "all"]
set n [molinfo top get numframes]

for {set i 1} {\$i < \$n} {incr i} {
    \$sel frame \$i   
    \$all frame \$i
    \$all move [measure fit \$sel \$ref]
    }

############################
#calc & out protein CA rmsd#
############################

set outfile [open "$OUT_DIR/${p}_${m}${n}_${r}_protein_rmsd.txt" w]
set ref [atomselect top "protein and name CA" frame 0]
set sel [atomselect top "protein and name CA"]
set n [molinfo top get numframes]
for {set i 1} {\$i < \$n} {incr i} {
  \$sel frame \$i
  set rmsd [measure rmsd \$sel \$ref]
  puts \$outfile "\$i \$rmsd"
  }
close \$outfile

############################
#calc & out epitope CA rmsd#
############################

set outfile [open "$OUT_DIR/${p}_${m}${n}_${r}_epitope_rmsd.txt" w]
set ref [atomselect top "protein and name CA and (resid $epitope)" frame 0]
set sel [atomselect top "protein and name CA and (resid $epitope)"]
for {set i 0} {\$i < \$n} {incr i} {
  \$sel frame \$i
  set rmsd [measure rmsd \$sel \$ref]
  puts \$outfile "\$i \$rmsd"
}
close \$outfile

#################
#calc & out rmsf#
#################

set outfile [open "$OUT_DIR/${p}_${m}${n}_${r}_rmsf.txt" w]
set sel [atomselect top "protein and name CA"]
set mol [\$sel molindex]
for {set i 0} {\$i < [\$sel num]} {incr i} {
    set rmsf [measure rmsf \$sel]
    puts \$outfile "[expr {\$i+1}] [lindex \$rmsf \$i]"
    }
close \$outfile

#################
#calc & out rgyr#
#################

set mol [molinfo top]
set outfile [open "$OUT_DIR/${p}_${m}${n}_${r}_rgyr.txt" w]
set sel [atomselect top "protein"]
set frames [molinfo \$mol get numframes]

for {set i 0} {\$i < \$frames} {incr i} {
    \$sel frame \$i
    \$sel update
    set a [measure rgyr \$sel]
    puts \$outfile "\$i \$a"
    }
close \$outfile

###############
#calc & out cn#
###############

set outfile [open "$OUT_DIR/${p}_${m}${n}_${r}_cn.txt" w]
set nf [molinfo top get numframes]
set ca [atomselect top "name CA"]

set nca [\$ca num]
puts \$outfile "\$nca"
for {set f 0} {\$f < \$nf} {incr f 25} {
    foreach index [\$ca get index] {
        set a [atomselect top "(all within 4.5 of index \$index) and ({name \"C.*\"} or {name \"N.*\"} or {name \"O.*\"} and not {name \"OH.*\"})"]
        \$a frame \$f
        \$a update
        set cn [\$a get index]
        puts \$outfile "\$f \$index \$cn"
    }
}
close \$outfile
exit
EOF

  eval "$VMD_EXE -e $OUT_DIR/$INPUT"
  rm $OUT_DIR/$INPUT
 
else
  echo "PROTEIN ${p}_${m}${n}_${r} CANNOT FOUND" >> $OUT_DIR/$OUTPUT
fi
	    done
        done
    done
done

$PYTHON_EXE native_contact.py


