/// Data cleaning part
1. git clone the repo
2. mkdir data
3. download raw folder
4. put raw into data folder
5. run .Rprofile
6. (if missing required package) renv::restore()
7. run rsrc/masterfile_all.R



8. installation instruction 1 to 4

9. python pysrc/sampling/baseline.py 
10. python pysrc/sampling/baseline.py --sites 78 

output: output/tables/gamma_percentiles_78.csv 
	output/tables/theta_percentiles_78.csv 
	output/tables/sigma_percentiles_78.csv
	output/tables/gamma_percentiles_1043.csv (Table 21)
	output/tables/theta_percentiles_1043.csv (Table 22)
	output/tables/sigma_percentiles_1043.csv (Table 23)



//// Figure 2 and Figure 3, Figure 4
11. run rsrc/analysis/calibration_maps_78_sites.R
12. run rsrc/analysis/calibration_maps_1043_sites.R

output: plots/calibration/1043SitesModel/map_z2017_1043Sites.png
	plots/calibration/1043SitesModel/map_x2017_1043Sites.png
	plots/calibration/1043SitesModel/map_gamma_fit.png
	plots/calibration/78SitesModel/map_gamma_fit.png
	plots/calibration/1043SitesModel/map_theta_fit.png
	plots/calibration/78SitesModel/map_theta_fit.png



/// Table 1
13. python pysrc/bash/shadow_price.py  (xi=1,5,10,1000 id =20...70) or use bash bash_files/shadow_price.sh





/// Figure 5, Figure 6, Table 2, Table 3, Table 4, Table 13
14. python scripts/conduction_det.py

output: output/tables/present_value_site1043_pa41.11_det.tex
	output/tables/present_value_site78_pa41.11_det.tex
	output/tables/transfer_cost_1043site_41.11pa_15year_det.tex
	output/tables/transfer_cost_1043site_41.11pa_30year_det.tex
	output/figures/pred_zshare_1043_sites_det.png
	output/figures/plot_pred_x_1043_sites_det.png
	output/figures/net_transfers.png




/// Figure 7,8 

15. R rsrc/analysis/map_1043_det.R
output: plots/1043-det/map_z0z30GammaTheta_1043Sites_allPrices_det.png
	plots/1043-det/map_zDecades_1043Sites_pe21.6_det.png



/// hmc 
16. python pysrc/bash/hmc_sampling.py --id ${id} --xi ${xi} --sites ${sites} --pee ${pee} (or use bash bash_files/hmc_sampling.sh)


17. python pysrc/bash/relative_entropy.py --xi ${xi} --sites ${sites} --pee ${pee} (or use bash_files/relative_entropy.sh)



18. python scripts/conduction_hmc.py

output: output/tables/present_value_site_ambiguity_comparison_xi_5.0.tex (Table 4)
	output/tables/present_value_site_ambiguity_comparison_xi_1.0.tex (Table 16)
	output/figures/aggregate_percentage_Z_b0_pehmc_6.6_pedet_6.6.png (Figure 11)
	output/figures/aggregate_percentage_Z_b0_pehmc_4.5_pedet_6.6.png (Figure 13)
	output/figures/aggregate_percentage_Z_b15_pehmc_19.5_pedet_21.6.png (Figure 13)
	output/figures/decision_histogram_pehmc_4.5_pedet_6.6_xi_5.0.png (Figure 12)
	output/figures/decision_histogram_pehmc_6.6_pedet_6.6_xi_5.0.png (Figure 10)
	output/figures/decision_histogram_pehmc_19.5_pedet_21.6_xi_5.0.png (Figure 12)
	output/tables/transfer_cost_1043site_41.11pa_15year_hmc_xi_5.0.tex (Table 14)
	output/tables/transfer_cost_1043site_41.11pa_30year_hmc_xi_5.0.tex (Table 15)








///mpc



19. bash bash_files/mpc_prepare.sh	 

20. bash bash_files/mpc_hmc_sp.sh

21. python pysrc/mpc/mpc_compute_sp.py

22. bash bash_files/mpc_hmc.sh

23. python pysrc/mpc/mpc_compute.py

24. python scripts/mpc_trajectory.py




///Appendix


//Table 10,11 and Figure 15
25. python scripts/price_estimation.py 

output: output/tables/hmm_results_table.tex
	output/tables/hmm_information_criteria
	output/figures/smooth_prob_con.png
	output/figures/smooth_prob_uncon.png



26. R rsrc/analysis/map_1043_hmc_xi05.R
    R rsrc/analysis/map_1043_hmc_xi1.R
output: plots/1043-hmc_xi5/map_zDecades_1043Sites_pe19.5_hmc.png (Figure 18)
	plots/1043-hmc_xi1/map_zDecades_1043Sites_pe17.2_hmc.png (Figure 19)


27. python scripts/bayesian_R2.py





