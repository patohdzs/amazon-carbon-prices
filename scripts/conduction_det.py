from pysrc.services.get_opt import get_optimization
from pysrc.analysis.figures import land_allocation,plot_transfers
from pysrc.analysis.tables import value_decom,transfer_cost

### Section 7.2 Results for case without stochasticity or ambiguity aversion


get_optimization(num_sites=1043,pee=6.6,model="det",solver="gurobi")
land_allocation(num_sites=1043,solver="gurobi",pee=6.6) # Figure 5
value_decom(num_sites=1043,pee=6.6,solver="gurobi") # Table 2
transfer_cost(num_sites=1043,pee=6.6,y=30,solver="gurobi") # Table 3
transfer_cost(num_sites=1043,pee=6.6,y=15,solver="gurobi") # Table 3
plot_transfers(num_sites=1043,pee=6.6,solver="gurobi") # Figure 6

get_optimization(num_sites=78,pee=6.0,model="det")
value_decom(num_sites=78,pee=6.0) # Table 4

print("det All done!")

