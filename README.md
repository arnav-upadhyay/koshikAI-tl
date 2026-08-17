# TXTLv2 — Translation model (background competition)

## 1. Installation
- Julia 1.10.3
- Clone / unzip this folder
- From the project root:

julia --project=. using Pkg Pkg.instantiate() # uses Project.toml + Manifest.toml


## 2. Load packages
Open `src/run_bg_sweep.ipynb` (or Julia REPL with cwd = `src/`) and run:
```julia
include("packages.jl")

## 3. Run call_model (increasing background)
In the same notebook, set N (or the list 5,10,15,20,25) and run the simulation cell. Outputs are written to data/output/.

## 4. Figures

Open src/figures.ipynb and run the plotting cells on the CSVs.