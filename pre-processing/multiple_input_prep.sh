if [ -d processed_input ]; then
	rm -rf processed_input;
	mkdir processed_input;
else
	mkdir processed_input;
fi

for i in $(ls pdb); do 

	filename=${i%.*};
	filetype=${i#*.}; 
	
	if [ $filetype == "pdb" ]; then
		echo "Process $filename ...";
		bash input_prep.sh $filename;
		mv $filename -t processed_input;
	fi;

done
