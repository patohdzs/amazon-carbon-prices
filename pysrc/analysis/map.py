
import geopandas as gpd
import os
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import numpy as np
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
import pandas as pd
from pysrc.services.file_service import get_path

def spatial_allocation(solver='gams',
                       pa=41.11,
                       pe_hmc=4.5,
                       pe_det=6.6,
                       b=0,
                       xi=5.0,
                       num_sites=1043):
    pe_hmc+=b
    pe_det+=b

    
    result_folder = os.path.join(str(get_path("output")), "optimization")
    output_folder = str(get_path("output"))
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    z_hmc= np.loadtxt(os.path.join(result_folder+ "/hmc/"+solver+f"/{num_sites}sites/xi_{xi}/pa_{pa}/pe_{pe_hmc}/", "Z.txt"), delimiter=",")
    z_hmc=z_hmc.T*1e9
    z_det= np.loadtxt(os.path.join(result_folder+ "/det/"+solver+f"/{num_sites}sites/pa_{pa}/pe_{pe_det}/", "Z.txt"), delimiter=",")
    z_det=z_det.T*1e9


    u_hmc= np.loadtxt(os.path.join(result_folder+ "/hmc/"+solver+f"/{num_sites}sites/xi_{xi}/pa_{pa}/pe_{pe_hmc}/", "U.txt"), delimiter=",")
    u_hmc=u_hmc.T[:,:50]*1e9
    v_hmc= np.loadtxt(os.path.join(result_folder+ "/hmc/"+solver+f"/{num_sites}sites/xi_{xi}/pa_{pa}/pe_{pe_hmc}/", "V.txt"), delimiter=",")
    v_hmc=v_hmc.T[:,:50]*1e9

    u_det= np.loadtxt(os.path.join(result_folder+ "/det/"+solver+f"/{num_sites}sites/pa_{pa}/pe_{pe_det}/", "U.txt"), delimiter=",")
    u_det=u_det.T[:,:50]*1e9
    v_det= np.loadtxt(os.path.join(result_folder+ "/det/"+solver+f"/{num_sites}sites/pa_{pa}/pe_{pe_det}/", "V.txt"), delimiter=",")
    v_det=v_det.T[:,:50]*1e9

    hmc_max = np.maximum(np.abs(u_hmc), np.abs(v_hmc))
    det_max = np.maximum(np.abs(u_det), np.abs(v_det))

    indicator_hmc = np.full(z_hmc.shape[0], 2)
    indicator_hmc[z_hmc[:, 49] > z_hmc[:, 0]] = 0
    indicator_hmc[z_hmc[:, 49] < z_hmc[:, 0]] = 1
    indicator_hmc[np.abs(z_hmc[:,49]-z_hmc[:,0])<0.05*z_hmc[:,0]]=2

    indicator_det = np.full(z_det.shape[0], 2)
    indicator_det[z_det[:, 49] > z_det[:, 0]] = 0
    indicator_det[z_det[:, 49] < z_det[:, 0]] = 1
    indicator_det[np.abs(z_det[:,49]-z_det[:,0])<0.05*z_det[:,0]]=2

    positions_hmc = []
    for i in range(hmc_max.shape[0]):
        if indicator_hmc[i] == 2:
            positions_hmc.append(111)
            continue
        series = hmc_max[i, :50]
        found = False  
        for j in range(len(series)):
            if series[j] == np.max(series) and series[j] > 1e-4:
                positions_hmc.append(j + 1)
                found = True
                break
        if not found:
            positions_hmc.append(111)

    positions_det = []
    for i in range(det_max.shape[0]):
        if indicator_det[i] == 2:
            positions_det.append(111)
            continue
        series = det_max[i, :50]
        found = False 
        for j in range(len(series)):
            if series[j] == np.max(series) and series[j] > 1e-4:
                positions_det.append(j + 1)
                found = True
                break
        if not found:
            positions_det.append(111)



    positions_hmc_array= np.array(positions_hmc)
    positions_det_array=np.array(positions_det)
    indicator_hmc_array=np.array(indicator_hmc)
    indicator_det_array=np.array(indicator_det)


    positions_hmc_array = np.where(indicator_hmc_array == 0, -positions_hmc_array, positions_hmc_array)
    positions_det_array = np.where(indicator_det_array == 0, -positions_det_array, positions_det_array)





    bin_edges = np.histogram_bin_edges(np.concatenate((positions_hmc_array[positions_hmc_array != 111], positions_det_array[positions_det_array != 111])), bins=40)

    plt.figure(figsize=(8, 6))
    plt.hist(positions_hmc_array[positions_hmc_array != 111],  bins=bin_edges, color='red', alpha=0.6,label='ambiguity averse')
    plt.hist(positions_det_array[positions_det_array != 111], bins=bin_edges, color='blue', alpha=0.6,label='ambiguity neutral')
    plt.xlabel("years (minus for deforestation)", fontsize=14)
    plt.ylabel("frequency", fontsize=14)
    if b==0:
        plt.xlim(-55,55)
    elif b==15:
        plt.xlim(0,40)
    plt.legend()
    plt.savefig(os.path.join(output_folder, f"figures/decision_histogram_pehmc_{pe_hmc}_pedet_{pe_det}_xi_{xi}.png"))
    plt.show()


    print(f"Pe_hmc is: {pe_hmc}, Pe_det is: {pe_det}, b is: {b}, xi is: {xi}")
    count_hmc = np.sum((positions_hmc_array == 111))
    count_det = np.sum((positions_det_array == 111))

    print(f"Number of observations missed for HMC: {count_hmc/1043}")
    print(f"Number of observations missed for DET: {count_det/1043}")
    

    
    return

