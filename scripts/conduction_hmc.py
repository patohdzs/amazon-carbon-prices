from pysrc.services.get_opt import get_optimization
from pysrc.services.get_sample import get_sampling
from pysrc.analysis.figures import land_allocation,density,trajectory_diff
from pysrc.analysis.tables import value_decom,transfer_cost,ambiguity_decom
from pysrc.analysis.map import spatial_allocation
from pysrc.analysis.relative_entropy import relative_entropy



## Section 7.3 Results with robustness to parameter uncertainty


## xi5
get_optimization(num_sites=1043,pee=6.6,model="det",solver="gams")
get_optimization(num_sites=1043,pee=4.5,model="hmc",xi=5.0,solver="gams")
get_optimization(num_sites=1043,pee=6.6,model="hmc",xi=5.0,solver="gams")

ambiguity_decom(num_sites=1043,pe_det=6.6,pe_hmc=4.5,xi=5.0,solver="gams") 
trajectory_diff(num_sites=1043,pe_hmc=6.6,pe_det=6.6,b=0,solver="gams",pa=41.11,xi=5.0) # Figure 11
trajectory_diff(num_sites=1043,pe_hmc=4.5,pe_det=6.6,b=0,solver="gams",pa=41.11,xi=5.0) # Figure 14
trajectory_diff(num_sites=1043,pe_hmc=4.5,pe_det=6.6,b=15,solver="gams",pa=41.11,xi=5.0) # Figure 14

transfer_cost(num_sites=1043,pee=4.5,xi=5.0,solver="gams",y=30,model="hmc") 
transfer_cost(num_sites=1043,pee=4.5,xi=5.0,solver="gams",y=15,model="hmc") 

spatial_allocation(num_sites=1043,pe_hmc=6.6,pe_det=6.6,xi=5.0,solver="gams",b=0) # Figure 12
spatial_allocation(num_sites=1043,pe_hmc=4.5,pe_det=6.6,xi=5.0,solver="gams",b=0) 
spatial_allocation(num_sites=1043,pe_hmc=4.5,pe_det=6.6,xi=5.0,solver="gams",b=15) 

### xi1
get_optimization(num_sites=1043,pee=2.2,model="hmc",xi=1.0,solver="gams")

ambiguity_decom(num_sites=1043,pe_det=6.6,pe_hmc=2.2,xi=1.0,solver="gams") 


print("hmc All done!")

