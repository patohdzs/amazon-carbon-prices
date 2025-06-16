import numpy as np
from pysrc.services.file_service import get_path
import os
from pysrc.services.data_service import load_site_data
import matplotlib.pyplot as plt


solver="gurobi"
num_sites=78
xi1=10000.0
xi2=1.0
xi3=0.5
pe1=6.3
pe2=6.0
pe3=5.6


(zbar, _, _) = load_site_data(num_sites)

def read_file(result_directory):
    
    Z = np.loadtxt(os.path.join(result_directory, "Z.txt"), delimiter=",")
    X = np.sum(np.loadtxt(os.path.join(result_directory, "X.txt"), delimiter=","),axis=1)
    Xdot = np.diff(X, axis=0)
    U = np.loadtxt(os.path.join(result_directory, "U.txt"), delimiter=",")
    V = np.loadtxt(os.path.join(result_directory, "V.txt"), delimiter=",")

    return (Z / 1e2, Xdot / 1e2, U / 1e2, V / 1e2)

result_zper1=[]
result_zper2=[]
result_zper3=[]

for i in range(50):


    result_directory1 = (
        str(get_path("output"))
        + f"/optimization/mpc/{solver}/{num_sites}sites/xi_{xi1}/pa_41.1/"
        + f"pe_{pe1}/mc_{i+1}/unconstrained"
    )

    result_directory2 = (
        str(get_path("output"))
        + f"/optimization/mpc_worstcase/{solver}/{num_sites}sites/xi_{xi2}/pa_41.1/"
        + f"pe_{pe2}/mc_{i+1}/unconstrained"
    )

    result_directory3 = (
        str(get_path("output"))
        + f"/optimization/mpc_worstcase/{solver}/{num_sites}sites/xi_{xi3}/pa_41.1/"
        + f"pe_{pe3}/mc_{i+1}/unconstrained"
    )

    (dfz_np1,_,_,_) = read_file(
        result_directory1
    )
    (dfz_np2,_,_,_) = read_file(
        result_directory2
    )

    (dfz_np3,_,_,_) = read_file(
        result_directory3
    )


    dfz_np1=np.sum(dfz_np1,axis=1)[:50]
    dfz_np2=np.sum(dfz_np2,axis=1)[:50]
    dfz_np3=np.sum(dfz_np3,axis=1)[:50]

    result_zper1.append((dfz_np1/ np.sum(zbar)) * 100)
    result_zper2.append((dfz_np2/ np.sum(zbar)) * 100)
    result_zper3.append((dfz_np3/ np.sum(zbar)) * 100)



result_zper1=np.mean(np.array(result_zper1),axis=0)
result_zper2=np.mean(np.array(result_zper2),axis=0)
result_zper3=np.mean(np.array(result_zper3),axis=0)


time = list(range(len(dfz_np1)))
output_folder = str(get_path("output"))
plt.figure()

plt.plot(time, result_zper1*100,linewidth=4,color = 'red', label=r"$\widehat{\xi}=\infty$")
plt.plot(time, result_zper2*100,linewidth=4,color = 'green', label=r"$\widehat{\xi}=1$")
plt.plot(time, result_zper3*100,linewidth=4,color = 'blue', label=r"$\widehat{\xi}=0.5$")

# plt.title("Land allocation, transfer b=0")
plt.ylabel("Z(%)")
plt.xlabel("years")
plt.ylim(12,24)
plt.xlim(0,50)
plt.legend()
plt.savefig(os.path.join(output_folder, f"figures/mpc_landallocation_b_{0}_adjust"))
plt.show()