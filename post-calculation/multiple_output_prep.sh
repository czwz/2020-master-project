if ! [ -d processed_output ]; then
	mkdir processed_output;
	cd processed_output;

	mkdir trj;
	for r in $(seq 1 1 5); do

		echo "processing sample $r..." >> process.log;
                mkdir trj/sample-$r;

		for filename in $(ls /scratch/snx3000/swang/md-input); do

			echo "processing $filename ..." >> process.log;
			mkdir trj/sample-$r/$filename;
			cp /scratch/snx3000/swang/md-post/tool/single_output_prep.sh -t /scratch/snx3000/swang/md-run/sample-$r/$filename;
			cp /scratch/snx3000/swang/md-post/tool/VMDextensions.tcl -t /scratch/snx3000/swang/md-run/sample-$r/$filename;

			cd /scratch/snx3000/swang/md-run/sample-$r/$filename;
			bash /scratch/snx3000/swang/md-run/sample-$r/$filename/single_output_prep.sh $filename;
			cd /scratch/snx3000/swang/md-post/processed_output;

                        mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.pnc r_${r}_$filename.pnc;
                        mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.enc r_${r}_$filename.enc;
                        mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.rgyr r_${r}_$filename.rgyr;
                        mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.rmsf r_${r}_$filename.rmsf;
                        mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.pprmsd r_${r}_$filename.pprmsd;
                        mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.eprmsd r_${r}_$filename.eprmsd;
                        mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.sprmsd r_${r}_$filename.sprmsd;
                        mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.eermsd r_${r}_$filename.eermsd;
                        mv /scratch/snx3000/swang/md-run/sample-$r/$filename/$filename.sermsd r_${r}_$filename.sermsd;

			cp /scratch/snx3000/swang/md-run/sample-$r/$filename/step3_charmm2gmx.pdb -t trj/sample-$r/$filename;
			mv /scratch/snx3000/swang/md-run/sample-$r/$filename/*noPBC.xtc -t trj/sample-$r/$filename;

			rm /scratch/snx3000/swang/md-run/sample-$r/$filename/view_state.vmd;
                        rm /scratch/snx3000/swang/md-run/sample-$r/$filename/vmd.log;
			rm /scratch/snx3000/swang/md-run/sample-$r/$filename/vmd.tcl;
			rm /scratch/snx3000/swang/md-run/sample-$r/$filename/single_output_prep.sh;
			rm /scratch/snx3000/swang/md-run/sample-$r/$filename/VMDextensions.tcl;

		done;
	done;
	echo "done." >> process.log;
	cd ../;
fi
