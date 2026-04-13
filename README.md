# Project Amazon

## Requirements
- Python >= 3.9, <3.12
- Gurobi >= 10.0.3

## Fresh Machine Mise en Place (macOS)

Use this section when starting from a brand-new Mac and a fresh clone.

### 1) Install Apple Command Line Tools
```bash
xcode-select --install
```

### 2) Install Homebrew (if missing)
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

### 3) Clone and enter repo
```bash
git clone <repo-url>
cd amazon-carbon-prices
```

### 4) Install pinned system dependencies
```bash
brew bundle --file Brewfile
```

This enforces the system stack expected by the project (`python@3.11`, `r`, `gdal`, `geos`, `proj`, `udunits`, `cmake`, `make`, `gcc`, etc.).

### 5) Run preflight checks
```bash
./scripts/preflight_macos.sh
```

If this fails, fix the reported issue first (toolchain path, SDK headers, missing Brew deps, or conflicting env vars).

### 6) Install and license Gurobi
1. Install Gurobi (>= 10.0.3) from Gurobi.
2. Activate your license:
```bash
grbgetkey <YOUR-LICENSE-KEY>
```
3. Verify installation:
```bash
gurobi_cl --version
```

### 7) Put raw data in place
You need:
```
data/raw/{esa,fgv,global_forest_watch,ibge,ipea,mapbiomas,seabpr,seeg,worldbank,worldclim}
```

### 8) Run project bootstrap
```bash
chmod +x run.sh
./run.sh -s
```

This setup stage:
- runs `brew bundle` and macOS preflight checks
- creates `.venv` with a compatible Python (`>=3.9,<3.12`, preferring `python3.11`)
- installs Python dependencies
- installs CmdStan (default pin: `2.37.0`)
- restores R dependencies from `renv.lock`

You can override the CmdStan version if needed:
```bash
CMDSTAN_VERSION=2.37.0 ./run.sh -s
```

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

## Installation (Manual Alternative)

If you do not use `./run.sh -s`, install manually:

1. Create and activate a virtual environment
```bash
python -m venv .venv
source .venv/bin/activate
```
2. Install Python dependencies
```bash
python -m pip install -e '.[all]'
```
3. Install CmdStan
```bash
install_cmdstan --version 2.37.0 --overwrite
```

If you get `CmdStanInstallError: Command "make build" failed`, run `./scripts/preflight_macos.sh` and fix the reported CLT/SDK/toolchain issue first.
4. Restore R dependencies
```bash
Rscript -e "renv::restore()"
```
5. (Contributors) install pre-commit hooks
```bash
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

### Long Weekend Run (No Sleep)
```bash
mkdir -p logs
LOG="logs/weekend_pipeline_$(date +%Y%m%d_%H%M%S).log"
nohup caffeinate -imsu bash -lc './run.sh -spcrDdtHhrMma' > "$LOG" 2>&1 &
echo $! > logs/weekend_pipeline.pid
echo "PID: $(cat logs/weekend_pipeline.pid)"
echo "LOG: $LOG"
```

Monitor:
```bash
tail -f "$LOG"
```

Stop:
```bash
kill "$(cat logs/weekend_pipeline.pid)"
```

Notes:
- Keep machine plugged into power.
- Use lowercase `-r` in the pipeline flag string.
- Use `-T` instead of `-t` for terminal-only time-consistency checks.

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
