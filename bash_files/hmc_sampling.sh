
# xiarray=(1 5 10000)
# idarray=(0 10 15 20 25)
# peearray=(4.5 6.6)


# xiarray=(1)
# idarray=(0 10 15 20 25)
# peearray=(2.2)

# xiarray=(5)
# idarray=(0 10 15 20 25)
# peearray=(4.5)

xiarray=(10000)
idarray=(0 10 15 20 25)
peearray=(2.2 4.5 6.6)

sites=1043

for pee in "${peearray[@]}"; do
    for xi in "${xiarray[@]}"; do
        for id in "${idarray[@]}"; do
   
                            count=0
                                        
                            action_name="hmc_sampling"

                            dataname="${action_name}"

                            mkdir -p ./job-outs/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/id_${id}/

                            if [ -f ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/id_${id}/run.sh ]; then
                                rm ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/id_${id}/run.sh
                            fi

                            mkdir -p ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/id_${id}/

                            touch ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/id_${id}/run.sh

                            tee -a ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/id_${id}/run.sh <<EOF
#!/bin/bash

#SBATCH --account=pi-lhansen
#SBATCH --job-name=id_${id}_${action_name}
#SBATCH --output=./job-outs/$job_name/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/id_${id}/run.out
#SBATCH --error=./job-outs/$job_name/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/id_${id}/run.err
#SBATCH --time=1-11:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=12G

module load python/anaconda-2022.05  

echo "\$SLURM_JOB_NAME"

echo "Program starts \$(date)"
start_time=\$(date +%s)

python3 -u /project/lhansen/HMC_conju/amazon-carbon-prices/pysrc/bash/hmc_sampling.py --id ${id} --xi ${xi} --sites ${sites} --pee ${pee}
echo "Program ends \$(date)"
end_time=\$(date +%s)
elapsed=\$((end_time - start_time))

eval "echo Elapsed time: \$(date -ud "@\$elapsed" +'\$((%s/3600/24)) days %H hr %M min %S sec')"

EOF
    count=$(($count + 1))
    sbatch ./bash/${action_name}/site_${sites}/pee_${pee}/xi_${xi}/id_${id}/run.sh
        done
    done
done
