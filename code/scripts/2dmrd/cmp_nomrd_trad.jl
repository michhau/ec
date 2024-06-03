######################################################
###          COMPARE TRADITIONAL AND NO-MRD        ###
###            author: Michi Haugeneder            ###
######################################################
#=
Compare 'traditional' MRD with non-orthogonal MRD.
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
animation = pyimport("matplotlib.animation")
gridspec = pyimport("matplotlib.gridspec")

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

turb.missing2nan!(evaldf1)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf3)
turb.missing2nan!(evaldf4)
turb.missing2nan!(evaldf5)
turb.missing2nan!(evaldf6)

evaldf = evaldf3

######################################################
###             TRADITIONAL MRD                    ###
######################################################
#uncomment to process the data
#Mulitresolution Flux Decomposition
##
M = 17
blocklength = ceil(2^M * timestep, Dates.Second)
(mrd_x, mrd_data, time_middle) = MRD.completemrd(evaldf, "w", "T", M, round(Int, 0.1 * 2^M))
(mrd_time, mrd_D_median, mrd_D_min, mrd_D_max, mrd_D_quant1, mrd_D_quant3) =
    MRD.mrdpp(mrd_x, mrd_data)
#------------------------------------------------------------------------
#=
#Plot result of MRD
PyPlot.pygui(true)
fig = PyPlot.figure()#figsize=(12,7))
ax = fig.add_subplot(111)
ax.set_title(string("MRD", evaldf1.time[1], " - ", evaldf1.time[end]))
ax.set_xlabel("avg. time [s]")
ax.set_ylabel(L"$C_{w\theta} [\cdot 10^{-3} \mathrm{Kms^{-1}}]$")
ax.grid(true)
PyPlot.xscale("log")
#ax.plot(Dates.value.(mrd_time) ./ 1000, mrd_D_median .* 1000, label="traditional")
#ax.fill_between(Dates.value.(mrd_time) ./ 1000, mrd_D_quant1 .* 1000, mrd_D_quant3 .* 1000, alpha=0.4)
ax.legend()
PyPlot.show()
=#
######################################################
###                    NO-MRD                      ###
######################################################
nrmodes1 = 50
nrmodes2 = M

(nomrd_x_1, nomrd_data_1, nomrd_time_middle_1) = MRD.completenomrd(evaldf, "w", "T", M, round(Int, 0.1 * 2^M), nrmodes1)
(nomrd_time_1, nomrd_D_median_1, nomrd_D_min_1, nomrd_D_max_1, nomrd_D_quant1_1, nomrd_D_quant3_1) =
    MRD.mrdpp(nomrd_x_1[:, 2], nomrd_data_1)
(nomrd_x_2, nomrd_data_2, nomrd_time_middle_2) = MRD.completenomrd(evaldf, "w", "T", M, round(Int, 0.1 * 2^M), nrmodes2)
(nomrd_time_2, nomrd_D_median_2, nomrd_D_min_2, nomrd_D_max_2, nomrd_D_quant1_2, nomrd_D_quant3_2) =
        MRD.mrdpp(nomrd_x_2[:, 2], nomrd_data_2)

#compare with 'traditional' MRD to check

to_plot_x1 = mrd_time
dmedian_1 = mrd_D_median
trad_q1 = mrd_D_quant1
trad_q3 = mrd_D_quant3
to_plot_x2 = nomrd_time_1
dmedian_2 = nomrd_D_median_1
new_q1 = nomrd_D_quant1_1
new_q3 = nomrd_D_quant3_1
to_plot_x3 = nomrd_time_2
dmedian_3 = nomrd_D_median_2
new_q1_3 = nomrd_D_quant1_2
new_q3_3 = nomrd_D_quant3_2

PyPlot.pygui(true)
fig = PyPlot.figure()#figsize=(8,5.5))
ax = fig.add_subplot(111)
ax.set_title(string("Comparison continuous - 'traditional' MRD"))
ax.set_xlabel("avg. time [s]")
ax.set_ylabel(L"$\mathrm{cospectral~density~w'T'}$")#\frac{w'T'}{\sum w'T'}}$")
ax.grid(true)
PyPlot.xscale("log")
ax.plot(Dates.value.(to_plot_x1) ./ 1000, dmedian_1, label="traditional")
ax.fill_between(Dates.value.(to_plot_x1) ./ 1000, trad_q1, trad_q3, alpha=0.4)
ax.plot(Dates.value.(to_plot_x2) ./ 1000, dmedian_2, label=string("continuous n = ", length(to_plot_x2)))
ax.fill_between(Dates.value.(to_plot_x2) ./ 1000, new_q1, new_q3, alpha=0.4)
ax.plot(Dates.value.(to_plot_x3) ./ 1000, dmedian_3, label=string("continuous n = ", length(to_plot_x3)))
ax.fill_between(Dates.value.(to_plot_x3) ./ 1000, new_q1_3, new_q3_3, alpha=0.4)
ax.legend()
PyPlot.show()

######################################################
###       COMPARE MEAN AND MEDIAN FOR NO-MRD       ###
######################################################
#=
compare the effect of mean and median of all the
windows in fastnomrd.
Therefore, need to change code in 
/src/mrd.jl/fastnomrd/second to last block
mean <-> median
=#

nrmodes3 = 60

@info("Code with 'median' (see comment)")
(nomrd_x_med, nomrd_data_med, nomrd_time_middle_med) = MRD.completenomrd(evaldf, "w", "T", M, round(Int, 0.1 * 2^M), nrmodes3)
(nomrd_time_med, nomrd_D_median_med, nomrd_D_min_med, nomrd_D_max_med, nomrd_D_quant1_med, nomrd_D_quant3_med) =
    MRD.mrdpp(nomrd_x_med[:, 2], nomrd_data_med)

@info("Code with 'mean' (see comment)")
(nomrd_x_mean, nomrd_data_mean, nomrd_time_middle_mean) = MRD.completenomrd(evaldf, "w", "T", M, round(Int, 0.1 * 2^M), nrmodes3)
(nomrd_time_mean, nomrd_D_median_mean, nomrd_D_min_mean, nomrd_D_max_mean, nomrd_D_quant1_mean, nomrd_D_quant3_mean) =
    MRD.mrdpp(nomrd_x_mean[:, 2], nomrd_data_mean)


to_plot_x1 = mrd_time
dmedian_1 = mrd_D_median
trad_q1 = mrd_D_quant1
trad_q3 = mrd_D_quant3
to_plot_x2 = nomrd_time_med
dmedian_2 = nomrd_D_median_med
new_q1 = nomrd_D_quant1_med
new_q3 = nomrd_D_quant3_med
to_plot_x3 = nomrd_time_mean
dmedian_3 = nomrd_D_median_mean
new_q1_3 = nomrd_D_quant1_mean
new_q3_3 = nomrd_D_quant3_mean

PyPlot.pygui(true)
fig = PyPlot.figure()#figsize=(8,5.5))
ax = fig.add_subplot(111)
ax.set_title(string("Comparison NO-MRD mean vs. median n = ", length(to_plot_x2)))
ax.set_xlabel("avg. time [s]")
ax.set_ylabel(L"$\mathrm{cospectral~density~w'T'}$")#\frac{w'T'}{\sum w'T'}}$")
ax.grid(true)
PyPlot.xscale("log")
ax.plot(Dates.value.(to_plot_x1) ./ 1000, dmedian_1, label="traditional")
ax.fill_between(Dates.value.(to_plot_x1) ./ 1000, trad_q1, trad_q3, alpha=0.4)
ax.plot(Dates.value.(to_plot_x2) ./ 1000, dmedian_2, label=string("cont. median"))
ax.fill_between(Dates.value.(to_plot_x2) ./ 1000, new_q1, new_q3, alpha=0.4)
ax.plot(Dates.value.(to_plot_x3) ./ 1000, dmedian_3, label=string("cont. mean"))
ax.fill_between(Dates.value.(to_plot_x3) ./ 1000, new_q1_3, new_q3_3, alpha=0.4)
ax.legend()
PyPlot.show()
