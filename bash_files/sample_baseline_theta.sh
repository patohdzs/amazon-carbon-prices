

idarray=($(seq 1 1))



for id in "${idarray[@]}"; do
   
                            count=0
                                        
                            action_name="baseline_theta"

                            dataname="${action_name}"

                            mkdir -p ./job-outs/${action_name}/id_${id}/

                            if [ -f ./bash/${action_name}/id_${id}/run.sh ]; then
                                rm ./bash/${action_name}/id_${id}/run.sh
                            fi

                            mkdir -p ./bash/${action_name}/id_${id}/

                            touch ./bash/${action_name}/id_${id}/run.sh

                            tee -a ./bash/${action_name}/id_${id}/run.sh <<EOF
#!/bin/bash

#SBATCH --account=pi-lhansen
#SBATCH --job-name=id_${id}_${action_name}
#SBATCH --output=./job-outs/$job_name/${action_name}/id_${id}/run.out
#SBATCH --error=./job-outs/$job_name/${action_name}/id_${id}/run.err
#SBATCH --time=0-36:00:00
#SBATCH --account=pi-lhansen
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=21G

module load python/anaconda-2022.05  

echo "\$SLURM_JOB_NAME"

echo "Program starts \$(date)"
start_time=\$(date +%s)

python3 -u /project/lhansen/HMC_final/project-amazon/scripts/sample_baseline_theta.py


echo "Program ends \$(date)"
end_time=\$(date +%s)
elapsed=\$((end_time - start_time))

eval "echo Elapsed time: \$(date -ud "@\$elapsed" +'\$((%s/3600/24)) days %H hr %M min %S sec')"

EOF
    count=$(($count + 1))
    sbatch ./bash/${action_name}/id_${id}/run.sh

done
