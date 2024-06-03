# Processing Eddy-Covariance Data

*This is a reorganization of scripts from my messy 'ibl_patch_snow' repository*

## Raw Data

1. import the raw logger data (not binary, but converted to ASCII using Campbell Scientific LoggerNet) using code/scripts/rawdata/read_raw_turb_data.jl (-> *.csv-file)
2. first preprocessing step using code/scripts/rawdata/offline_preproc.jl (-> *.csv-file or *.nc-file)
    
    Offline Preprocessing includes:
    1. Check according to physical limits (u,v,w,T,co2,h2o)
    2. Check quality control flags of device (CSAT & IRGASON)
    3. Make dataset time-continuos (adding rows with time stamps and 'missing' data, if needed)
    4. Despiking according to algorithm in Sigmund et al. (2022) (https://doi.org/10.1007/s10546-021-00653-x)
    5. save as *.csv or *.nc

## Evaluation

### load_data.jl

use /code/scripts/load_data.jl to load data in specific time period (evalstart - evalend) including following preprocessing
1. apply a NaN mask where data is not usable
2. show how much of data is missing or NaN
3. interpolate gaps if not too long (threshold can be set; default: 1s)
4. Double rotation (default: blocks of 30min length)

optional: load slow data

### mrd_script.jl
/code/scripts/mrd_script.jl: calculate Multiresolution Flux Decompositions. For references see Howell&Mahrt (1997) (https://doi.org/10.1023/A:1000210427798) and Vickers&Mahrt (2003) (https://doi.org/10.1175/1520-0426(2003)20<660:TCGATF>2.0.CO;2)

### turb_fluxes.jl

use /code/scripts/turb_fluxes.jl to calculate turbulent fluxes (need to know Reynolds averaging time before) and Obukhov lengths. Fluxes are calculated by subtracting a moving average from the 20Hz-data and smoothing afterwards. This script also includes multiple functionalities for plotting the data.

### Flux Footprints

All the scripts in this repository are adapted (or simply converted from Python to Julia) from Kljun et al. (2015) (https://doi.org/10.5194/gmd-8-3695-2015). Climatology aggregates single footprints over a longer period.

### 2D-MRDs
Moving-window multiresolution flux decomposition to see the time dependence of scales of turbulent variables. The general methodology is described in Haugeneder et al. (2024) (https://doi.org/10.1007/s10546-023-00856-4).

- 2dmrd_distributed.jl: Calculate the 2D-MRD on multiple cores (e.g. on HYPERION) using the Julia package 'Distributed'. Change variables 'path_to_data_to_use' (points to data source), 'cols_for_mrd' (which variables/flux should be decomposed), and 'savefilename' (according to the flux variable)
- stitch_2dmrds.jl: Stitch single 2D-MRD files (*.nc) in time domain.
- combine_2dmrds_xxx.jl: Combine 2D-MRDs to obtain a decomposition of a composed variable (as friction velocity or turbulent kinetic energy)
- cmp_nomrd_trad.jl: Compare non-orthogonal MRD with traditional (orthogonal) MRD (see Fig. 4 in Haugeneder et al. (2024))
- plot_2dmrd.jl: Collection of various scripts to plot the 2D-MRDs


### Other
- code/scripts/block_evaluation.jl: get starttime, endtime, block duration, winddir, mean wind speed, mean wT, mean T per block and plot

