from pysrc.services.get_opt import get_optimization
from pysrc.analysis.figures import land_allocation,plot_transfers
from pysrc.analysis.tables import value_decom,transfer_cost

### Section 7.2 Results for case without stochasticity or ambiguity aversion


get_optimization(num_sites=1043,pee=6.6,model="det",solver="gams")
land_allocation(num_sites=1043,solver="gams",pee=6.6) # Figure 5
value_decom(num_sites=1043,pee=6.6,solver="gams") # Table 2
transfer_cost(num_sites=1043,pee=6.6,y=30,solver="gams") # Table 3
transfer_cost(num_sites=1043,pee=6.6,y=15,solver="gams") # Table 3
plot_transfers(num_sites=1043,pee=6.6,solver="gams") # Figure 6


get_optimization(num_sites=78,pee=6.1,model="det",solver="gams")
value_decom(num_sites=78,pee=6.1,solver="gams") # Table 4



print("det All done!")

