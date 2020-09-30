for r $(seq 1 1 3); do

	echo "processing sample $r ..."
	cd sample $r;

	for filename in $(ls /scratch/snx3000/swang/md-input); do
	
		cd $filename;
		bash submit-continue.sh;
		echo "continuing $filename ...";
		cd ../

	done;

done
