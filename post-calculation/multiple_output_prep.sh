if ! [ -d processed_data ]; then
	mkdir processed_data;
	cd processed_data;

	mkdir trj;
	#for r in $(seq 1 1 5); do
	for r in $(seq 1 1 2); do

		echo "processing sample $r..."
                mkdir trj/sample-$r;

		#for filename in $(ls /scratch/snx3000/swang/md-input); do
		for filename in "3E2H_b_1" "3E2H_g_1"; do

			echo "processing $filename ...";
			mkdir trj/sample-$r/$filename;
			cp /scratch/snx3000/swang/md-post/tool/single_output_prep.sh -t /scratch/snx3000/swang/md-run/sample-$r/$filename;

			cd /scratch/snx3000/swang/md-run/sample-$r/$filename;
			bash /scratch/snx3000/swang/md-run/sample-$r/$filename/single_output_prep.sh $filename;
			cd -;

			mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.nc r_${r}_$filename.nc;
                        mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.rgyr r_${r}_$filename.rgyr;
			mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.rmsf r_${r}_$filename.rmsf;
			mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.ermsd r_${r}_$filename.ermsd;
			mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.prmsd r_${r}_$filename.prmsd;
			cp /scratch/snx3000/swang/md-run/sample-$r/$filename/step3_charmm2gmx.pdb -t trj/sample-$r/$filename;
			mv /scratch/snx3000/swang/md-run/sample-$r/$filename/*noPBC.xtc -t trj/sample-$r/$filename;

			rm /scratch/snx3000/swang/md-run/sample-$r/$filename/view_state.vmd;
			rm /scratch/snx3000/swang/md-run/sample-$r/$filename/vmd.tcl;
			rm /scratch/snx3000/swang/md-run/sample-$r/$filename/single_output_prep.sh;
		done;
	done;
	cd ../;
fi
