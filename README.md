# Data Sharing — Liquid Holdup and Mass Transfer in a Vertical Droplet Array

Data and code behind the manuscript "Predicting Mass Transfer in a Vertical Droplet Array" (AIChE Journal; revised version, November 2026).

The mass transfer and liquid holdup data are raw. The droplet data are given per imaging location (24 image stacks:
12 locations in the array × top and mid windows) and as the composite file used to compute the average array
characteristics. Model results cover the 18 experimental trials (2 flow rates × 3 CO2 levels × 3 trials).

`PD_gqmom_0_1.m` was shared by Alberto Passalacqua and should be cited as:
Fox, R. O., Laurent, F., & Passalacqua, A. (2023). The generalized quadrature method of moments.
*Journal of Aerosol Science*, 167, 106096. https://doi.org/10.1016/j.jaerosci.2022.106096

All paths in the scripts are relative. **Run every script with the working directory set to the folder the script is in.**

## Folder layout

```
Data Sharing/
├── DataProcessingCodeShare.R          Main analysis: experimental CO2 / liquid holdup results vs. model predictions,
│                                      uncertainties, and the manuscript figures (run from Data Sharing/)
├── Experimental Data/
│   ├── MassTransfer Experimental Data 10222025.xlsx    CO2 concentration measurements (one sheet per flow rate / CO2 level)
│   ├── Liquid Holdup Experimental Data 10222025.xlsx   Liquid holdup measurements
│   ├── LH_Experimental.csv                              Liquid holdup data used by DataProcessingCodeShare.R
│   ├── Droplet Data Checks/                             Manually checked per-stack droplet tables, 09/24/2024 (4.5 L/min)
│   └── Stack Analysis/
│       ├── 09242024/   StackAnalysis.R outputs, 4.5 L/min (droplet distribution used by the model)
│       └── 09232024/   StackAnalysis.R outputs, 3.0 L/min imaging rig (sensitivity check only)
├── Image Processing Files/
│   ├── createMask.m, dropletArray.m   MATLAB droplet detection on the high-speed image stacks
│   └── StackAnalysis.R                Combines the checked per-stack tables into the droplet dataset
│                                      (writes <date>_AllDropletData.csv, Frames.csv, Droplets.csv and summaries)
├── Droplet Tower Modeling Files/      MATLAB droplet tower model (run from this folder)
│   ├── runAll.m                Runs every trial; writes Model Results/<flow> CO2 <level> 10272026/T#/MT_Models.csv,
│   │                           <flow>LPM_LH_t_model_10272026.csv and PerTrialInputs_10272026.csv
│   ├── RunTrial.m              Model for one trial (GQMOM size distribution + ODE integration)
│   ├── trialInputs.m           Physical properties per trial (Supplement table)
│   ├── ODEsys_Amokrane.m, ODEsys_Angelo.m, ODEsys_Hsu.m, ODEsys_Ruckenstein.m   Mass transfer correlations
│   ├── PD_gqmom_0_1.m          Gaussian quadrature method of moments
│   ├── runUncertainty.m        Model uncertainty   -> Model Results/ModelUncertainty_10272026.csv
│   ├── runQuadratureCheck.m    Quadrature node check -> Model Results/QuadratureCheck_10272026.csv
│   ├── runSensitivityBatch.m   Sensitivity studies -> Model Results/*_10282026.csv
│   └── HistogramPlot.m         Droplet size histogram vs. GQMOM nodes
├── Model Results/                     Outputs of the scripts above (10272026 / 10282026 runs)
├── Supporting Analyses/
│   ├── SupportingCalculations.R (+ .csv)           Short calculations supporting the reviewer responses
│   ├── Velocity Analysis/VelocityAnalysis.R (+ CSV outputs)        Droplet velocity / acceleration from tracked drops
│   └── Aspect Ratio Analysis/AspectRatioAnalysis.R (+ CSV outputs) Droplet shape oscillation
```

## Workflow

1. **Image processing** — `dropletArray.m` detects droplets in each image stack; each stack's table is checked by hand
   (`Experimental Data/Droplet Data Checks`). `StackAnalysis.R` combines them into
   `Experimental Data/Stack Analysis/<date>/`. The raw `.tiff` image stacks are not included. Run on the shared files
   alone, `StackAnalysis.R` reproduces `09242024_AllDropletData.csv`, `Frames.csv`, `Droplets.csv` and the size
   statistics exactly; the columns that need the raw images (frames per stack, drops per 100 frames, and the
   deletion / manual-check tables) are left blank and are provided in the shared summary files instead.
2. **Model** — `runAll.m` reads `09242024_AllDropletData.csv` and the mass transfer data and writes `Model Results/`.
   `runUncertainty.m`, `runQuadratureCheck.m` and `runSensitivityBatch.m` add the uncertainty and sensitivity results.
3. **Analysis** — `DataProcessingCodeShare.R` compares experiments with `Model Results/`. The `Supporting Analyses`
   scripts read `Experimental Data/Stack Analysis/<date>/<date>_Frames.csv`.

## Changes from the June 2026 version

- Model code restructured: `NDFApproach_dropletTowerModel_09242024PDF.m` and `LiquidHoldup_Time_*_Model_09242024PDF.m`
  were replaced by `runAll.m` + `RunTrial.m` + `trialInputs.m`; `ODEsys_*.m` updated; uncertainty, quadrature and
  sensitivity scripts added.
- Model results: the `10222025` runs were replaced by the `10272026` runs; liquid holdup model files now
  `<flow>LPM_LH_t_model_10272026.csv` in `Model Results/`.
- Droplet data: `Experimental Data/AllDropletData.csv` and `dataCleaning.R` were replaced by `StackAnalysis.R` and
  `Experimental Data/Stack Analysis/09242024/09242024_AllDropletData.csv` (corrected distribution). Three checked
  stack files were updated.
- `DataProcessingCodeShare.R` updated to the revised analysis (instrument uncertainty, model uncertainty, 10272026 results).
- New: `Supporting Analyses/`, 3.0 L/min imaging-rig data (`Stack Analysis/09232024`).

The previous files are kept in `_Previous version (June 2026)/`.
