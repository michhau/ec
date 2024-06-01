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


