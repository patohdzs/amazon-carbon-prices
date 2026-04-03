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


8. Follow installation instructions (steps 1 to 4)




9. Run the baseline script:
```bash
python pysrc/sampling/baseline.py
```

10. Run with 78 sites:
```bash
python pysrc/sampling/baseline.py --sites 78
```



11. Run:
```R
rsrc/analysis/calibration_maps_78_sites.R
```

12. Run:
```R
rsrc/analysis/calibration_maps_1043_sites.R
```



13. Run the deterministic model (shadow prices, optimization, maps):
```bash
./run.sh -d
```

14. Run time-consistency checks:
```bash
./run.sh -t
```

15. Run the HMC ambiguity model (adjusted sampling, shadow prices, relative entropy):
```bash
./run.sh -h
```

16. Run the MPC price-risk model:
```bash
./run.sh -m
```

17. Run remaining analysis (calibration maps, price estimation, Bayesian R2, HMC maps):
```bash
./run.sh -a
```


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
