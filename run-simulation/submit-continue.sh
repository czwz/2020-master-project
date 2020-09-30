job1=$(sbatch --account=s959 gmx-prod-continuation.job);
old_job=$(echo "${job1//[!0-9]/}" );
echo "submit pc: $old_job ..." >> submit_log;

for i in $(seq 1 1 5); do
        old_id=$old_job
        new_job=$(sbatch --account=s959 --dependency=afterany:$old_job gmx-prod-continuation.job);
        old_job=$(echo "${new_job//[!0-9]/}" );
        echo "submit pc : $old_job with dependency on $old_id..." >> submit_log;
done

