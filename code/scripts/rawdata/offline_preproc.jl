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

datapath = "/home/haugened/Documents/data/"
importdir = joinpath(@__DIR__, "..", "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .gen
PyPlot.pygui(true)

#########################################################################
#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
kaijo_period_file = joinpath(datapath, "kaijo", "21_kaijo_periods.txt")
kaijo_outfile_stam = joinpath(datapath, "kaijo", "kaijo")
tower_outfile_stam = joinpath(datapath, "tower")
outfile_stam = joinpath(datapath, "WindSonic_Duerrboden_2005xx/tjk_sonic_200504-200514")
tower_outfile_stam = joinpath(datapath, "tower")

writepath = joinpath(tower_outfile_stam, "preproc/")
#########################################################################

#########################################################################
###                    I/O AND QUALITY CONTROL                        ###
#########################################################################
println()
println("-----------S-T-A-R-T-------------")
##########################################################################
#uncomment to read-in already created .csv files to DataFrame
#t1irg = turb.csvtodataframe(joinpath(tower_outfile_stam, "t1irg.csv"))
#t2irgdb = turb.csvtodataframe(joinpath(tower_outfile_stam, "t2irg_db.csv"))
t2lcsatdb = turb.csvtodataframe(joinpath(tower_outfile_stam, "t2lcsat_db.csv"))
#t2ucsatdb = turb.csvtodataframe(joinpath(tower_outfile_stam, "t2ucsat_db.csv"))
#t2irgpd = turb.csvtodataframe(joinpath(tower_outfile_stam, "t2irg_pd.csv"))
#t2lcsatpd = turb.csvtodataframe(joinpath(tower_outfile_stam, "t2lcsat_pd.csv"))
#t2ucsatpd = turb.csvtodataframe(joinpath(tower_outfile_stam, "t2ucsat_pd.csv"))

#kaijodf = turb.csvtodataframe(string(kaijo_outfile_stam, ".csv"))
tjkdf = turb.csvtodataframe("/home/haugened/Documents/data/tjk/sonic_complete/tjk_sonic.csv")
#sonicdf = turb.csvtodataframe(string(outfile_stam, ".csv"))

#kaijodf = turb.sonicqualcontrol(kaijodf,-5.001,5.001,-5.001,6.001,-5.501,5.501,5.001,27.001)
#t1irg = turb.sonicqualcontrol(t1irg)
#turb.qualcontrolflags!(t1irg)
#t1irgdb = t1irg[DateTime(2021,05,19,12,00,00).<t1irg.time, :]
#t2irgdb = turb.sonicqualcontrol(t2irgdb)
#turb.qualcontrolflags!(t2irgdb)
#t2irgdb = t2irgdb[DateTime(2021,05,19,12,00,00).<t2irgdb.time, :]
#t2lcsatdb = turb.sonicqualcontrol(t2lcsatdb)
#turb.qualcontrolflags!(t2lcsatdb)
#t2lcsatdb = t2lcsatdb[DateTime(2021,05,19,12,00,00).<t2lcsatdb.time, :]
#t2ucsatdb = turb.sonicqualcontrol(t2ucsatdb)
#turb.qualcontrolflags!(t2ucsatdb)
#t2ucsatdb = t2ucsatdb[DateTime(2021,05,19,12,00,00).<t2ucsatdb.time, :]
#t2irgpd = turb.sonicqualcontrol(t2irgpd)
#t2lcsatpd = turb.sonicqualcontrol(t2lcsatpd)
#t2ucsatpd = turb.sonicqualcontrol(t2ucsatpd)
#sonicdf = turb.sonicqualcontrol(sonicdf)
tjkdf = turb.sonicqualcontrol(tjkdf, -15.00001, 15.00001, -15.00001,
    15.00001, -2.500001, 2.500001, -21.00001, 21.00001, 0, 25)

####################################################
#which DataFrame should be preprocessed
to_preproc = t2lcsatdb[60*20-1:end, :]
#preselect time period
to_preproc = to_preproc[DateTime(2021, 05, 20, 13, 07, 00).<=to_preproc.time.<=DateTime(2021, 06, 11, 08, 17, 00, 00), :] #to_preproc[249:end, :] #to_preproc[249:end, :]
to_preproc = turb.makecontinuous(to_preproc)
dur = Minute(to_preproc.time[end] - to_preproc.time[1])
@show dur

#despiking after Sigmund et al. (2022)
after_desp = turb.despiking(to_preproc)

#interpolate missing values
#dataout = turb.interpolatemissing(after_desp)

#write to output array
CSV.write(joinpath(writepath, "tjkdf.csv"), after_desp)

####################################################
###            CONVERT TO NetCDF4                ###
####################################################

turb.saveturbasnetcdf(tjkdf, "/home/haugened/Documents/data/tower/preproc/tjkdf.nc")#joinpath(writepath, "t2ucsatdb.nc"))