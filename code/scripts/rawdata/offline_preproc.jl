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

datapath = "/home/michi/Documents/slf/CONTRASTS25/data/processed/"
importdir = joinpath(@__DIR__, "..", "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .gen
PyPlot.pygui(true)

#########################################################################
#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
outfile_stam = joinpath(datapath, "preproc")
#tower_outfile_stam = joinpath(datapath, "tower")

writepath = joinpath(tower_outfile_stam, "preproc/")
#########################################################################

#########################################################################
###                    I/O AND QUALITY CONTROL                        ###
#########################################################################
println()
println("-----------S-T-A-R-T-------------")
##########################################################################
irg = turb.csvtodataframe(joinpath(datapath, "3b_t1_irg.csv"))
csat = turb.csvtodataframe(joinpath(datapath, "3b_t1_csat.csv"))

#fix a wrong naming in the CONTRASTS data
rename!(irg, :h20signalstrength => :h2osignalstrength)

irg = turb.qualcontrolflags(irg);
csat = turb.qualcontrolflagCSAT3(csat);

irg = turb.sonicqualcontrol(irg);
csat = turb.sonicqualcontrol(csat);

####################################################
#which DataFrame should be preprocessed
to_preproc_irg = irg;#t2lcsatdb[60*20-1:end, :]
to_preproc_csat = csat;
#preselect time period
#to_preproc = to_preproc[DateTime(2021, 05, 20, 13, 07, 00).<=to_preproc.time.<=DateTime(2021, 06, 11, 08, 17, 00, 00), :] #to_preproc[249:end, :] #to_preproc[249:end, :]
to_preproc_irg = turb.makecontinuous(to_preproc_irg);
to_preproc_csat = turb.makecontinuous(to_preproc_csat);
dur_irg = Minute(to_preproc_irg.time[end] - to_preproc_irg.time[1]);
@show dur_irg
dur_csat = Minute(to_preproc_csat.time[end] - to_preproc_csat.time[1]);
@show dur_csat

#despiking after Sigmund et al. (2022)
after_desp_irg = turb.despiking(to_preproc_irg);
after_desp_csat = turb.despiking(to_preproc_csat);

#interpolate missing values
#dataout = turb.interpolatemissing(after_desp)

#write to output array
#CSV.write(joinpath(writepath, "2a_box2_irg_proc.csv"), after_desp)

####################################################
###            CONVERT TO NetCDF4                ###
####################################################

turb.saveturbasnetcdf(after_desp_irg, joinpath(writepath, "3b_t1_irg_proc.nc"))#joinpath(writepath, "t2ucsatdb.nc"))
turb.saveturbasnetcdf(after_desp_csat, joinpath(writepath, "3b_t1_csat_proc.nc"))#joinpath(writepath, "t2ucsatdb.nc"))