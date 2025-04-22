from pysrc.services.get_opt import get_optimization
from pysrc.analysis.figures import land_allocation
from pysrc.analysis.tables import value_decom,transfer_cost




get_optimization(num_sites=78,pee=5.9,model="det",solver="gams",pa=35.71)

value_decom(num_sites=78,pee=5.9,solver="gams",pa=35.71)

get_optimization(num_sites=78,pee=6.3,model="det",solver="gams",pa= 44.26)

value_decom(num_sites=78,pee=6.3,solver="gams",pa= 44.26)

# get_optimization(num_sites=78,pee=6.1,model="det",solver="gams",pa= 41.11)

# value_decom(num_sites=78,pee=6.1,solver="gams",pa= 41.11)