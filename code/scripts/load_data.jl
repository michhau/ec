######################################################
###          LOADING DATA FROM EC SENSORS          ###
###            author: Michi Haugeneder            ###
######################################################
#=
Load the preprocessed (by offline_preproc.jl)
turbulence data from the Eddy-Covariance sensors.
Do further necessary preprocessing steps.
Further evaluation in other scripts.

evalstart and evalend to change period to import
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter, Distributed
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")

datapath = "/home/haugened/Documents/data/"
tjkpath = "/home/haugened/Documents/data/tjk/tjk_data.csv"
importdir = joinpath(@__DIR__, "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .gen
PyPlot.pygui(true)

######################################################
#variables

#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
kaijo_period_file = joinpath(datapath, "kaijo", "21_kaijo_periods.txt")
kaijo_outfile_stam = joinpath(datapath, "kaijo", "kaijo")
outfile_stam = joinpath(datapath, "WindSonic_Duerrboden_2005xx/tjk_sonic_200504-200514")
tower_outfile_stam = joinpath(datapath, "tower", "preproc")
ventair_stam = joinpath(datapath, "tower", "vent_air")
######################################################

println()
println("-----------S-T-A-R-T-------------")

######################################################
###            LOADING & PREPROCESSING             ###
######################################################
#select data and measurement period to be evaluated
evalstart = DateTime(2021, 05, 31, 13, 00, 00)
evalend   = DateTime(2021, 05, 31, 13, 30, 00)
#evalend = evalstart + Day(10)

evaldf1 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t1irgdb.nc"), evalstart, evalend)
evaldf2 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2irgdb.nc"), evalstart, evalend)
evaldf3 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2lcsatdb.nc"), evalstart, evalend)
evaldf4 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2ucsatdb.nc"), evalstart, evalend)
evaldf5 = turb.readturbasnetcdf(string(kaijo_outfile_stam, ".nc"), evalstart, evalend)
evaldf6 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "tjkdf.nc"), evalstart, evalend)

#rotate kaijo data due to mounting. Be careful and check for every new position!
if @isdefined evaldf5
    @warn("Careful!! Rotate Kaijo data due to measurement position. Check pictures of setup!")
    new_w = evaldf5.u
    new_u = -evaldf5.w
    evaldf5.u = new_u
    evaldf5.w = new_w
end

#apply NaN-mask to T1 & T2 when repositioned (for DR)
turb.repositionnanmask!(evaldf1)
turb.repositionnanmask!(evaldf2)
turb.repositionnanmask!(evaldf3)
turb.repositionnanmask!(evaldf4)

#show statistics about missing data
turb.printmissstats(evaldf1)
turb.printmissstats(evaldf2)
turb.printmissstats(evaldf3)
turb.printmissstats(evaldf4)
turb.printmissstats(evaldf5)
turb.printmissstats(evaldf6)

#interpolate missing
evaldf1 = turb.interpolatemissing(evaldf1)
evaldf2 = turb.interpolatemissing(evaldf2)
evaldf3 = turb.interpolatemissing(evaldf3)
evaldf4 = turb.interpolatemissing(evaldf4)
try
    evaldf5 = turb.interpolatemissing(evaldf5)
catch LoadError
    @warn("Kaijo DataFrame empty")
end
evaldf6 = turb.interpolatemissing(evaldf6)

#double rotation
turb.drdf!(evaldf1)
turb.drdf!(evaldf2)
turb.drdf!(evaldf3)
turb.drdf!(evaldf4)
turb.drdf!(evaldf5)
turb.drdf!(evaldf6, periodwise=false)

#turb.saveturbasnetcdf(evaldf5, "/home/haugened/Documents/openfoam/duerr_les/scripts/src/kaijotmp.nc")

######################################################
###               LOADING SLOW DATA                ###
######################################################

tjkmeteodata = turb.csvtodataframe(tjkpath)
tjkmeteodata = tjkmeteodata[evalstart.<=tjkmeteodata.time.<=evalend, :]

t2vent = turb.readturbasnetcdf(joinpath(ventair_stam, "t2.nc"), evalstart, evalend)
t3vent = turb.readturbasnetcdf(joinpath(ventair_stam, "t3.nc"), evalstart, evalend)

#remove unrealistic values
t2vent.vent_air_temp[.!(-30 .< replace!(t2vent.vent_air_temp, missing => NaN) .< 30)] .= NaN
t3vent.vent_air_temp[.!(-30 .< replace!(t3vent.vent_air_temp, missing => NaN) .< 30)] .= NaN
#=
######################################################
###                   PLOTTING                     ###
######################################################
#plot histograms to determine outliers

quantity1 = filter(!isnan, skipmissing(evaldf1.T))
quantity2 = filter(!isnan, skipmissing(evaldf2.T))
quantity3 = filter(!isnan, skipmissing(evaldf3.T))
quantity4 = filter(!isnan, skipmissing(evaldf4.T))
quantity5 = filter(!isnan, skipmissing(evaldf5.w))
binmin = -8
binmax = 8

##
comap = PyPlot.get_cmap("tab10");
hist = PyPlot.figure()
axh = hist.add_subplot(111)
axh.set_title("Kajio measurements (raw)")
axh.set_xlabel(L"w~\mathrm{[m~s^{-1}]}")
axh.grid()
PyPlot.yscale("log")
#axh.hist(vec(quantity1), bins=collect(binmin:binmax), density=true, label="T1IRG", alpha=0.5, color=comap(0))
#axh.hist(vec(quantity2), bins=collect(binmin:binmax), density=true, label="T2IRG", alpha=0.5, color=comap(1))
#axh.hist(vec(quantity3), bins=collect(binmin:binmax), density=true, label="T2LCSAT", alpha=0.5, color=comap(2))
#axh.hist(vec(quantity4), bins=collect(binmin:binmax), density=true, label="T2UCSAT", alpha=0.5, color=comap(3))
axh.hist(vec(quantity5), bins=collect(binmin:0.01:binmax), density=true, label="KAIJO", color=comap(4))
axh.legend()
##

##
#plot single parameter (e.g. u)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf4)
sifig = PyPlot.figure()
siax = sifig.add_subplot(111)
#siax.set_title(L"\mathrm{1~min~gliding~average}", fontsize=12)
siax.set_xlabel("03.06.2021")
siax.set_ylabel("winddir (w/o DR)")
#siax.set_ylabel(L"\overline{u}~\mathrm{[m~s^{-1}]}")
#siax.plot(evaldf1.time, gen.movingaverage(evaldf1.u,20*60), label="T1IRG")
siax.plot(evaldf2.time, gen.movingaverage(turb.simplewinddir.(evaldf2.u, evaldf2.v), 20*600), label="T2IRG, 10min avg")
#siax.plot(evaldf3.time, gen.movingaverage(evaldf3.u,20*60), label="T2LCSAT")
siax.plot(evaldf4.time, gen.movingaverage(turb.simplewinddir.(evaldf4.u, evaldf4.v), 20*600), label="T2UCSAT 10min avg")
#siax.plot(evaldf5.time, gen.movingaverage(evaldf5.u,20*60), label="KAIJO")
majorlocator = pydates.HourLocator(interval=1)
minorlocator = pydates.MinuteLocator([15,30,45])
siax.xaxis.set_major_locator(majorlocator)
siax.xaxis.set_minor_locator(minorlocator)
date_format = pydates.DateFormatter("%H:%M")
siax.xaxis.set_major_formatter(date_format)
#sifig.autofmt_xdate()
siax.grid()
siax.legend()
##
=#
#=
#plot time series with moving average
fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.plot(evaldf2.time[1:40:end], gen.movingaverage(evaldf2.u,40)[1:40:end], label="1m")
ax.plot(evaldf3.time[1:40:end], gen.movingaverage(evaldf3.u,40)[1:40:end], label="2m")
ax.plot(evaldf4.time[1:40:end], gen.movingaverage(evaldf4.u,40)[1:40:end], label="3m")
ax.plot(evaldf5.time[1:40:end], gen.movingaverage(evaldf5.u,40)[1:40:end], label="0.3m")
ax.plot(evaldf6.time[1:40:end], gen.movingaverage(evaldf6.u,40)[1:40:end], label="5m")
ax.set_ylabel(L"u~\mathrm{[m~s^{-1}]}")
ax.grid()
ax.legend()
=#