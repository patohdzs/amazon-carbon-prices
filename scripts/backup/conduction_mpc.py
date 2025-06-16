from pysrc.analysis.tables import value_decom_mpc
from pysrc.sampling.mpc_parallel import parallel_conduction
from pysrc.sampling.mpc_simulation import mpc_simulation


### simulation mpc samples
mpc_simulation()

### get mpc optimization results
for pe in [6.9,21.9,31.9]:
    for id in range(1,101,5):
        parallel_conduction(idstart=id,idend=id+4,pf=pe,model="unconstrained")
    print(f"MPC optimization finished! pe={pe}")

### then start plotting tables and figures

value_decom_mpc(pee=6.9,b=0,opt='gams')
value_decom_mpc(pee=6.9,b=15,opt='gams')
value_decom_mpc(pee=6.9,b=25,opt='gams')

print("MPC All done!")