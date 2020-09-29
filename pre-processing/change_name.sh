for i in $(ls); do 

	new=$(echo $i | sed s/_enriched/_g/g \
		      | sed s/_enrich/_g/g \
		      | sed s/_noenriched/_b/g \
		      | sed s/_noenrich/_b/g \
		      | sed s/_best/_/g \
		      | sed s/_1st/_1/g \
		      | sed s/_worst/_/g);
	mv $i $new 
done
