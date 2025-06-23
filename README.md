# Project Amazon

## Requirements
- Python >= 3.9, <3.12
- Gurobi >= 10.0.3
## Data Requirements
To replicate, make sure to download the raw data into the directory structure below:
```
.
└── data
    └── raw
        ├── esa
        ├── fgv
        ├── global_forest_watch
        ├── ibge
        ├── ipea
        ├── mapbiomas
        ├── seabpr
        ├── seeg
        ├── worldbank
        └── worldclim
```

## Installation

0. Clone git repository and move into `project-amazon/`
1. Create and activate a new virtual environment
```
python -m venv .venv
source .venv/bin/activate
```
2. Install python dependencies
```
python -m pip install -e '.[all]'
```

3. Install CmdStan
```
install_cmdstan --overwrite
```

4. Install pre-commit hooks (required for contributors)
```
pre-commit install
```

## Replication


1. Clone the repository
```bash
git clone <repo-url>
```

2. Create the data folder
```bash
mkdir data
```

3. Download the `raw` data folder

4. Move `raw` into `data` folder

5. Run `.Rprofile`

6. Restore R packages if missing
```R
renv::restore()
```

7. Run masterfile
```R
rsrc/masterfile_all.R
```

---


8. Follow installation instructions (steps 1 to 4)

---


9. Run the baseline script:
```bash
python pysrc/sampling/baseline.py
```

10. Run with 78 sites:
```bash
python pysrc/sampling/baseline.py --sites 78
```

---

11. Run:
```R
rsrc/analysis/calibration_maps_78_sites.R
```

12. Run:
```R
rsrc/analysis/calibration_maps_1043_sites.R
```

---

13. Run:
```bash
python pysrc/bash/shadow_price.py
```
Or:
```bash
bash bash_files/shadow_price.sh
```

---


14. Run:
```bash
python scripts/conduction_det.py
```

---

15. Run:
```R
rsrc/analysis/map_1043_det.R
```

---

16. Run:
```bash
python pysrc/bash/hmc_sampling.py --id ${id} --xi ${xi} --sites ${sites} --pee ${pee}
```
Or:
```bash
bash bash_files/hmc_sampling.sh
```

17. Run:
```bash
python pysrc/bash/relative_entropy.py --xi ${xi} --sites ${sites} --pee ${pee}
```
Or:
```bash
bash bash_files/relative_entropy.sh
```

18. Run:
```bash
python scripts/conduction_hmc.py
```

---


19. Prepare MPC structure:
```bash
bash bash_files/mpc_prepare.sh
```

20. Run MPC HMC sampling (shadow price):
```bash
bash bash_files/mpc_hmc_sp.sh
```

21. Compute shadow prices:
```bash
python pysrc/mpc/mpc_compute_sp.py
```

22. Run MPC HMC sampling:
```bash
bash bash_files/mpc_hmc.sh
```

23. Compute results:
```bash
python pysrc/mpc/mpc_compute.py
```

24. Plot MPC trajectories:
```bash
python scripts/mpc_trajectory.py
```

---


25. Run:
```bash
python scripts/price_estimation.py
```

26. Run:
```R
rsrc/analysis/map_1043_hmc_xi05.R
rsrc/analysis/map_1043_hmc_xi1.R
```

27. Run:
```bash
python scripts/bayesian_R2.py



## Contributing
0. Open a new git branch
```
git checkout -b <new branch name>
```
1. Create new code changes

2. Stage changed files
```
git add <names of changed files>
```

3. Commit changed files
```
git commit
```

4. After several commits, push commits to remote (if it is the first time pushing this branch use the `--set-upstream` flag)
```
git push
```

5. Submit a pull request on GitHub
