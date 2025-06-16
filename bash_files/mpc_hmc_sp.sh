

xiarray=(0.5 1 10000)
pearray=($(seq 5.0 0.1 6.3))
id=999

typearray=("unconstrained" "constrained")

for type in "${typearray[@]}"; do
    for xi in "${xiarray[@]}"; do
        for pe in "${pearray[@]}"; do
    
                            count=0
                                        
                            action_name="mpc_sp"

                            dataname="${action_name}"

                            mkdir -p ./job-outs/${action_name}/xi_${xi}/pe_${pe}/id_${id}/type_${type}

                            if [ -f ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/type_${type}/run.sh ]; then
                                rm ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/type_${type}/run.sh
                            fi

                            mkdir -p ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/type_${type}/

                            touch ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/type_${type}/run.sh

                            tee -a ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/type_${type}/run.sh <<EOF
#!/bin/bash

#SBATCH --account=pi-lhansen
#SBATCH --job-name=id_${id}_${action_name}
#SBATCH --output=./job-outs/$job_name/${action_name}/xi_${xi}/pe_${pe}/id_${id}/type_${type}/run.out
#SBATCH --error=./job-outs/$job_name/${action_name}/xi_${xi}/pe_${pe}/id_${id}/type_${type}/run.err
#SBATCH --time=1-11:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G

module load python/anaconda-2022.05  
module load gurobi/11.0

echo "\$SLURM_JOB_NAME"

echo "Program starts \$(date)"
start_time=\$(date +%s)

python3 -u /project/lhansen/HMC_0525/amazon-carbon-prices/pysrc/mpc/mpc_hmc_sp.py --pe ${pe} --xi ${xi} --type ${type}
echo "Program ends \$(date)"
end_time=\$(date +%s)
elapsed=\$((end_time - start_time))

eval "echo Elapsed time: \$(date -ud "@\$elapsed" +'\$((%s/3600/24)) days %H hr %M min %S sec')"

EOF
    count=$(($count + 1))
    sbatch ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/type_${type}/run.sh
        done
    done
done
