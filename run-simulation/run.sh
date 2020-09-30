for r in $(seq 1 1 3); do

	echo "processing sample $r..."
        if ! [ -d sample-$r ]; then

                mkdir sample-$r;
		cd sample-$r;

		for filename in $(ls /scratch/snx3000/swang/md-input); do

			mkdir $filename;
			cp -r /scratch/snx3000/swang/md-input/$filename/* -t $filename;
			cp /scratch/snx3000/swang/md-run/tool/* -t $filename;

			cd $filename;
			bash submit-single.sh;
			echo "sending $filename ...";
			cd ../;

		done;

		cd ../;
		squeue -u swang;
	fi;
done
