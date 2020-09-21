
job1=$(sbatch gmx-min.job);
echo "sbatch gmx-min.job";
job1=$(echo "${job1//[!0-9]/}" );

job2=$(sbatch --dependency=afterany:$job1 gmx-equil.job);
echo "sbatch --dependency=afterany:${job1} gmx-equil.job";
job2=$(echo "${job2//[!0-9]/}" );

job3=$(sbatch --dependency=afterany:$job2 gmx-prod-1st.job);
echo "sbatch --dependency=afterany:${job2} gmx-prod-1st.job";
job3=$(echo "${job3//[!0-9]/}" );

job4=$(sbatch --dependency=afterany:$job3 gmx-prod-continuation.job);
echo "sbatch --dependency=afterany:${job3} gmx-prod-continuation.job";
job4=$(echo "${job4//[!0-9]/}" );

squeue -u abriata
