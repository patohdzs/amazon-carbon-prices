from pysrc.services.get_opt import get_optimization
from pysrc.analysis.figures import trajectory_diff
from pysrc.analysis.tables import transfer_cost,ambiguity_decom
from pysrc.analysis.map import spatial_allocation
from pysrc.analysis.figures import density


## Section 7.3 Results with robustness to parameter uncertainty

pe_det_base = 6.8

## xi1
get_optimization(num_sites=1043,pee=pe_det_base,model="det",solver="gurobi")
get_optimization(num_sites=1043,pee=4.8,model="hmc",xi=1.0,solver="gurobi")
get_optimization(num_sites=1043,pee=pe_det_base,model="hmc",xi=1.0,solver="gurobi")

ambiguity_decom(num_sites=1043,pe_det=pe_det_base,pe_hmc=4.8,xi=1.0,solver="gurobi")
trajectory_diff(num_sites=1043,pe_hmc=pe_det_base,pe_det=pe_det_base,b=0,solver="gurobi",pa=41.11,xi=1.0) # Figure 11
trajectory_diff(num_sites=1043,pe_hmc=4.8,pe_det=pe_det_base,b=0,solver="gurobi",pa=41.11,xi=1.0) # Figure 14
trajectory_diff(num_sites=1043,pe_hmc=4.8,pe_det=pe_det_base,b=15,solver="gurobi",pa=41.11,xi=1.0) # Figure 14

transfer_cost(num_sites=1043,pee=4.8,xi=1.0,solver="gurobi",y=30,model="hmc")
transfer_cost(num_sites=1043,pee=4.8,xi=1.0,solver="gurobi",y=15,model="hmc")


spatial_allocation(num_sites=1043,pe_hmc=pe_det_base,pe_det=pe_det_base,xi=1.0,solver="gurobi",b=0) # Figure 12
spatial_allocation(num_sites=1043,pe_hmc=4.8,pe_det=pe_det_base,xi=1.0,solver="gurobi",b=0)
spatial_allocation(num_sites=1043,pe_hmc=4.8,pe_det=pe_det_base,xi=1.0,solver="gurobi",b=15) 





# xi0_5
get_optimization(num_sites=1043,pee=2.8,model="hmc",xi=0.5,solver="gurobi")
ambiguity_decom(num_sites=1043,pe_det=pe_det_base,pe_hmc=2.8,xi=0.5,solver="gurobi") 



# xi2
get_optimization(num_sites=1043,pee=5.6,model="hmc",xi=2.0,solver="gurobi")
ambiguity_decom(num_sites=1043,pe_det=pe_det_base,pe_hmc=5.6,xi=2.0,solver="gurobi") 


density(num_sites=1043,pee=4.8,xi=1.0,solver="gurobi")
density(num_sites=1043,pee=5.6,xi=2.0,solver="gurobi")
density(num_sites=1043,pee=2.8,xi=0.5,solver="gurobi")


print("hmc All done!")
