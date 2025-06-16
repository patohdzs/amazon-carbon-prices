


xiarray=(0.5 1 2 10000)


# xiarray=(0.5)


sites=1043

for xi in "${xiarray[@]}"; do
  if [[ "$xi" == "0.5" ]]; then
    idarray=($(seq 20 39))
  elif [[ "$xi" == "1" ]]; then
    idarray=($(seq 40 50))
  elif [[ "$xi" == "2" ]]; then
    idarray=($(seq 50 60))
  else
    idarray=($(seq 60 70))
  fi
    for id in "${idarray[@]}"; do
   
                            count=0
                                        
                            action_name="shadow_price"

                            dataname="${action_name}"

                            mkdir -p ./job-outs/${action_name}/xi_${xi}/id_${id}/sites_${sites}

                            if [ -f ./bash/${action_name}/xi_${xi}/id_${id}/run.sh ]; then
                                rm ./bash/${action_name}/xi_${xi}/id_${id}/run.sh
                            fi

                            mkdir -p ./bash/${action_name}/xi_${xi}/id_${id}/

                            touch ./bash/${action_name}/xi_${xi}/id_${id}/run.sh

                            tee -a ./bash/${action_name}/xi_${xi}/id_${id}/run.sh <<EOF
#!/bin/bash

#SBATCH --account=pi-lhansen
#SBATCH --job-name=id_${id}_${action_name}
#SBATCH --output=./job-outs/$job_name/${action_name}/xi_${xi}/id_${id}/sites_${sites}/run.out
#SBATCH --error=./job-outs/$job_name/${action_name}/xi_${xi}/id_${id}/sites_${sites}/run.err
#SBATCH --time=1-11:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=12G

module load python/anaconda-2022.05  
module load gurobi/11.0
source .venv/bin/activate

echo "\$SLURM_JOB_NAME"

echo "Program starts \$(date)"
start_time=\$(date +%s)

python -u /project/lhansen/HMC_0525/amazon-carbon-prices/pysrc/bash/shadow_price.py --id ${id} --xi ${xi} --sites ${sites} 
echo "Program ends \$(date)"
end_time=\$(date +%s)
elapsed=\$((end_time - start_time))

eval "echo Elapsed time: \$(date -ud "@\$elapsed" +'\$((%s/3600/24)) days %H hr %M min %S sec')"

EOF
    count=$(($count + 1))
    sbatch ./bash/${action_name}/xi_${xi}/id_${id}/run.sh

    done
done
