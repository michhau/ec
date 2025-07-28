######################################################
###        OFFLINE PREPROCESSING OF EC DATA        ###
###            author: Michi Haugeneder            ###
######################################################
#=
Script to do preprocessing of the raw EC data:
- quality control with physical limits
- generate continuous data series
- spike removal (after Sigmund et al. (2022))
- export as .csv for further analysis
=#
using Dates, PyCall, DataFrames, Statistics, ProgressMeter, NCDatasets
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")

datapath = "/home/michi/Documents/slf/CONTRASTS25/data/processed/2a/raw/"
importdir = joinpath(@__DIR__, "..", "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .gen
PyPlot.pygui(true)

#########################################################################
#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
outfile_stam = joinpath(datapath, "..")
#tower_outfile_stam = joinpath(datapath, "tower")

writepath = joinpath(outfile_stam)
#########################################################################

#########################################################################
###                    I/O AND QUALITY CONTROL                        ###
#########################################################################
println()
println("-----------S-T-A-R-T-------------")
##########################################################################
irg = turb.csvtodataframe(joinpath(datapath, "2a_box1_irg.csv"))
csat = turb.csvtodataframe(joinpath(datapath, "2a_box1_csat.csv"))

irg = turb.sonicqualcontrol(irg)
csat = turb.sonicqualcontrol(csat)

####################################################
#which DataFrame should be preprocessed
to_preproc = csat#t2lcsatdb[60*20-1:end, :]
#preselect time period
#to_preproc = to_preproc[DateTime(2021, 05, 20, 13, 07, 00).<=to_preproc.time.<=DateTime(2021, 06, 11, 08, 17, 00, 00), :] #to_preproc[249:end, :] #to_preproc[249:end, :]
to_preproc = turb.makecontinuous(to_preproc)
dur = Minute(to_preproc.time[end] - to_preproc.time[1])
@show dur

#despiking after Sigmund et al. (2022)
after_desp = turb.despiking(to_preproc)

#interpolate missing values
#dataout = turb.interpolatemissing(after_desp)

#write to output array
CSV.write(joinpath(writepath, "2a_box1_csat_proc.csv"), after_desp)

####################################################
###            CONVERT TO NetCDF4                ###
####################################################

turb.saveturbasnetcdf(after_desp, joinpath(writepath, "2a_box1_csat_proc.nc"))#joinpath(writepath, "t2ucsatdb.nc"))