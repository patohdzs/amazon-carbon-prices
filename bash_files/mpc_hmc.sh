


pee=6.3
xi=10000
trig=0
type="unconstrained"
idarray=($(seq 1 50))
pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))



# pee=6.0
# xi=1
# trig=1
# type="unconstrained"
# idarray=($(seq 1 50))
# pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))


# pee=5.6
# xi=0.5
# trig=1
# type="unconstrained"
# idarray=($(seq 1 50))
# pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))


# pee=6.0
# xi=10000
# trig=0
# type="constrained"
# idarray=($(seq 1 50))
# pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))



# pee=5.7
# xi=1
# trig=1
# type="constrained"
# idarray=($(seq 1 50))
# pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))


# pee=5.2
# xi=0.5
# trig=1
# type="constrained"
# idarray=($(seq 1 50))
# pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))



for id in "${idarray[@]}"; do
    for pe in "${pearray[@]}"; do
   
                            count=0
                                        
                            action_name="mpc"

                            dataname="${action_name}"

                            mkdir -p ./job-outs/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}

                            if [ -f ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.sh ]; then
                                rm ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.sh
                            fi

                            mkdir -p ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}

                            touch ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.sh

                            tee -a ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.sh <<EOF
#!/bin/bash

#SBATCH --account=pi-lhansen
#SBATCH --job-name=id_${id}_${action_name}
#SBATCH --output=./job-outs/$job_name/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.out
#SBATCH --error=./job-outs/$job_name/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.err
#SBATCH --time=1-11:00:00
#SBATCH --partition=caslake	
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

module load python/anaconda-2022.05  
module load gurobi/11.0
source .venv/bin/activate

echo "\$SLURM_JOB_NAME"

echo "Program starts \$(date)"
start_time=\$(date +%s)

python3 -u /project/lhansen/HMC_0525/amazon-carbon-prices/pysrc/mpc/mpc_hmc.py --id ${id} --pe ${pe} --xi ${xi} --trig ${trig} --type ${type}
echo "Program ends \$(date)"
end_time=\$(date +%s)
elapsed=\$((end_time - start_time))

eval "echo Elapsed time: \$(date -ud "@\$elapsed" +'\$((%s/3600/24)) days %H hr %M min %S sec')"

EOF
    count=$(($count + 1))
    sbatch ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.sh

    done
done
