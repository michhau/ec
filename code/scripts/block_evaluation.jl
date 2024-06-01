######################################################
###        GET INFORMATION ON BLOCKS OF DATA       ###
###            author: Michi Haugeneder            ###
######################################################
#=
NOTE: Read data with load_data.jl
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
#animation = pyimport("matplotlib.animation")
#GridSpec = pyimport("matplotlib.gridspec")
#mpwidgets = pyimport("matplotlib.widgets")

importdir = joinpath(@__DIR__, "..")
datapath = "/home/haugened/Documents/data/"
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .gen

PyPlot.pygui(true)

#########################################################################
# change variables here
#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)

println()
println("-----------S-T-A-R-T-------------")
#########################################################################
###                        DATA SELECTION                             ###
#########################################################################
#select data and measurement period to be evaluated
evaldf = copy(evaldf5)
#evaldf = t1irg

evalstart2 = DateTime(2021,04,28,10,57,07)
#evalend = DateTime(2021,04,23,12,30,00)
evalend2 = evalstart2 + Minute(10)
evaldf = evaldf[evalstart2 .<= evaldf[:,1] .<= evalend2,:]

#detrend
evaldf.u = turb.detrend(evaldf.u)
evaldf.v = turb.detrend(evaldf.v)
evaldf.w = turb.detrend(evaldf.w)

#=
#double rotation if needed (data was already double rotated in load_data.jl)
@info("Data should already have been double rotated in load_data.jl. Only do it again, if other period choosen")
turb.drdf!(evaldf)
=#

##
#=
#for Kaijo
#read-in information about the measurement periods
measureperiod = turb.readperiodfile(kaijo_period_file)
#select measurement period number for following evaluation (check DataFrame for times)
sel_per = 4

evaldf = evaldf[measureperiod[sel_per,1] .<= evaldf[:,1] .<= measureperiod[sel_per,2],:]
=#

#########################################################################
#get information about data blocks
(kaijo_block, datarot) = turb.blockevaluation(evaldf, Second(5), timestep)

#write block DataFrame to file
#CSV.write(string(outfile_stam, "_block.csv"), eval_block)

#########################################################################
#=
#read in the already exported files
kaijodf = turb.csvtodataframe(string(kaijo_outfile_stam, ".csv"))
eval_block = turb.csvtodataframe(string(outfile_stam, "_block.csv"))
=#
#########################################################################
###                           Plotting                                ###
#########################################################################
#plot block properties
turb.missing2nan!(evaldf)
fig2 = PyPlot.figure()
ax = fig2.add_subplot(111)
ax.set_title("Block data")
ax.set_xlabel("time")
ax.set_ylabel(L"\left(\overline{w'T_s'}\right)_{30s} \mathrm{[\frac{Km}{s}]}")
ax.grid(true)
#ax.set_ylim((-0.15,0.25))
#ax.plot(kaijo_block.time_middle, kaijo_block.wT, label="kaijo")
#ax.plot(evaldf.time, evaldf.w)
#ax.plot(t1irg_block.time_middle, t1irg_block.wT, label="T1IRG")
#ax.plot(t2irg_block.time_middle, t2irg_block.wT, label="T2IRG")
#ax.plot(t2lcsat_block.time_middle, t2lcsat_block.wT, label="lower CSAT", alpha=0.5)
#ax.plot(t2ucsat_block.time_middle, t2ucsat_block.wT, label="upper CSAT", alpha=0.5)
ax.legend()

#########################################################################
