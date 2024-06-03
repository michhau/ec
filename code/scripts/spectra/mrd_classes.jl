######################################################
###MULTI RESOLUTION FLUX DECOMPOSITION FOR CLASSES ###
###            author: Michi Haugeneder            ###
######################################################
#=
Calculate MRD for data grouped in classes (for example)
according to wind direction
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
animation = pyimport("matplotlib.animation")

importdir = joinpath(@__DIR__, "..", "..")
datapath = "/home/haugened/Documents/data/"
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
import .turb
import .gen
import .MRD
PyPlot.pygui(true)

timestep = Millisecond(50)

evaldf = copy(evaldf6)

turb.missing2nan!(evaldf)

######################################################
###           SET CONDITION ON DATA                ###
######################################################
#wind direction class
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[evalstart.<=tjkmeteo[:, 1].<=evalend, :]
disallowmissing!(tjkdata)
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 1
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 2

ectime = collect(evaldf.time[1]:timestep:evaldf.time[end])
ecdata = fill(NaN, size(ectime, 1))
turb.extendtofinertimeseries!(ecdata, ectime, windclass, tjkdata.time)

#apply condition
condition = ecdata .== 1
evaldf[.!condition, [:u, :v, :w, :T]] .= NaN

#remove NaNs
nans = isnan.(evaldf.u)
data_tmp = evaldf[.!nans, :]


######################################################
###             MRD  VERCAUTEREN                   ###
######################################################
#uncomment to process the data
M = 15

D = MRD.mrd_std_block(data_tmp, evaldf.time[1], evaldf.time[end], 15, "u", "w", 2)

fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.set_title("TJK up valley wind, M=15")
mrd_std_plot = ax.pcolormesh(collect(0.5:14.5), collect(0.5:14.5), D)
ax.set_xlabel("w^2-scale [2^x data points]")
ax.set_ylabel("u-scale [2^x data points]")
cbar = fig.colorbar(mrd_std_plot, ax=ax)#, orientation = "horizontal")

######################################################
###                  NO-MRD                        ###
######################################################
M = 15
blocklength = ceil(2^M * timestep, Dates.Second)
(mrd_x, mrd_data, time_middle) = MRD.completemrd(evaldf, "w", "T", M, round(Int, 0.1 * 2^M))
(mrd_time, mrd_D_median, mrd_D_min, mrd_D_max, mrd_D_quant1, mrd_D_quant3) =
    MRD.mrdpp(mrd_x, mrd_data)

fig = PyPlot.figure()#figsize=(12,7))
ax = fig.add_subplot(111)
ax.set_title("T2LCSAT up valley")
ax.set_xlabel("time scale [s]")
ax.set_ylabel(L"$C_{wT_s} [10^{-3} \mathrm{Kms^{-1}}]$")
ax.grid(true)
PyPlot.xscale("log")
ax.plot(Dates.value.(mrd_time) ./ 1000, mrd_D_median .* 1000)
ax.fill_between(Dates.value.(mrd_time) ./ 1000, mrd_D_quant1 .* 1000, mrd_D_quant3 .* 1000, alpha=0.4)