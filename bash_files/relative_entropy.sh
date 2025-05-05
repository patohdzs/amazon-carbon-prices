
# xiarray=(5.0)
# pee=4.5
# sites=1043

xiarray=(1.0)
pee=2.2
sites=1043

for xi in "${xiarray[@]}"; do
    
   
                            count=0
                                        
                            action_name="relative_entropy"

                            dataname="${action_name}"

                            mkdir -p ./job-outs/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/

                            if [ -f ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/run.sh ]; then
                                rm ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/run.sh
                            fi

                            mkdir -p ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/

                            touch ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/run.sh

                            tee -a ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/run.sh <<EOF
#!/bin/bash

#SBATCH --account=pi-lhansen
#SBATCH --job-name=id_${id}_${action_name}
#SBATCH --output=./job-outs/$job_name/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/run.out
#SBATCH --error=./job-outs/$job_name/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/run.err
#SBATCH --time=1-11:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G

module load python/anaconda-2022.05  

echo "\$SLURM_JOB_NAME"

echo "Program starts \$(date)"
start_time=\$(date +%s)

python3 -u /project/lhansen/HMC_conju/amazon-carbon-prices/pysrc/bash/relative_entropy.py --xi ${xi} --sites ${sites} --pee ${pee}
echo "Program ends \$(date)"
end_time=\$(date +%s)
elapsed=\$((end_time - start_time))

eval "echo Elapsed time: \$(date -ud "@\$elapsed" +'\$((%s/3600/24)) days %H hr %M min %S sec')"

EOF
    count=$(($count + 1))
    sbatch ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/run.sh


done
