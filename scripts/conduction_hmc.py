from pysrc.services.get_opt import get_optimization
from pysrc.analysis.figures import trajectory_diff
from pysrc.analysis.tables import transfer_cost,ambiguity_decom
from pysrc.analysis.map import spatial_allocation



## Section 7.3 Results with robustness to parameter uncertainty


### xi1
get_optimization(num_sites=1043,pee=6.6,model="det",solver="gurobi")
get_optimization(num_sites=1043,pee=4.7,model="hmc",xi=1.0,solver="gurobi")
get_optimization(num_sites=1043,pee=6.6,model="hmc",xi=1.0,solver="gurobi")

ambiguity_decom(num_sites=1043,pe_det=6.6,pe_hmc=4.7,xi=1.0,solver="gurobi") 
trajectory_diff(num_sites=1043,pe_hmc=6.6,pe_det=6.6,b=0,solver="gurobi",pa=41.11,xi=1.0) # Figure 11
trajectory_diff(num_sites=1043,pe_hmc=4.7,pe_det=6.6,b=0,solver="gurobi",pa=41.11,xi=1.0) # Figure 14
trajectory_diff(num_sites=1043,pe_hmc=4.7,pe_det=6.6,b=15,solver="gurobi",pa=41.11,xi=1.0) # Figure 14

transfer_cost(num_sites=1043,pee=4.7,xi=1.0,solver="gurobi",y=30,model="hmc") 
transfer_cost(num_sites=1043,pee=4.7,xi=1.0,solver="gurobi",y=15,model="hmc") 


spatial_allocation(num_sites=1043,pe_hmc=6.6,pe_det=6.6,xi=1.0,solver="gurobi",b=0) # Figure 12
spatial_allocation(num_sites=1043,pe_hmc=4.7,pe_det=6.6,xi=1.0,solver="gurobi",b=0) 
spatial_allocation(num_sites=1043,pe_hmc=4.7,pe_det=6.6,xi=1.0,solver="gurobi",b=15) 


print("hmc All done!")


