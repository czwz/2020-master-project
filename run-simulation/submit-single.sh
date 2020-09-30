job1=$(sbatch --account=s959 gmx-min.job);
job1=$(echo "${job1//[!0-9]/}" );
echo "submit min: $job1 ..." >> submit_log;

job2=$(sbatch --account=s959 --dependency=afterany:$job1 gmx-equil.job);
job2=$(echo "${job2//[!0-9]/}" );
echo "submit eq : $job2 with dependency on $job1..." >> submit_log;

job3=$(sbatch --account=s959 --dependency=afterany:$job2 gmx-prod-1st.job);
old_job=$(echo "${job3//[!0-9]/}" );
echo "submit p1 : $old_job with dependency on $job2..." >> submit_log;

for i in $(seq 1 1 9); do
	old_id=$old_job
	new_job=$(sbatch --account=s959 --dependency=afterany:$old_job gmx-prod-continuation.job);
	old_job=$(echo "${new_job//[!0-9]/}" );
	echo "submit pc : $old_job with dependency on $old_id..." >> submit_log;
done
