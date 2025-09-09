######################################################
###       MULTI RESOLUTION FLUX DECOMPOSITION      ###
###            author: Michi Haugeneder            ###
######################################################
#=
Calculate 'normal' and quasi-continuos MRD. Load data with
load_data.jl
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

turb.missing2nan!(evaldf1)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf3)
turb.missing2nan!(evaldf4)
#turb.missing2nan!(evaldf5)
#turb.missing2nan!(evaldf6)

######################################################
###                  MRD                           ###
######################################################
#Mulitresolution Flux Decomposition
##

M = 17
quantity1 = "u"
quantity2 = "w"

blocklength = ceil(2^M * timestep, Dates.Second)
(mrd_x_1, mrd_data_1, time_middle_1) = MRD.completemrd(evaldf1, quantity1, quantity2, M, round(Int, 0.1 * 2^M))
(mrd_time_1, mrd_D_median_1, mrd_D_min_1, mrd_D_max_1, mrd_D_quant1_1, mrd_D_quant3_1) =
    MRD.mrdpp(mrd_x_1, mrd_data_1)
(mrd_x_2, mrd_data_2, time_middle_2) = MRD.completemrd(evaldf2, quantity1, quantity2, M, round(Int, 0.1 * 2^M))
(mrd_time_2, mrd_D_median_2, mrd_D_min_2, mrd_D_max_2, mrd_D_quant1_2, mrd_D_quant3_2) =
    MRD.mrdpp(mrd_x_2, mrd_data_2)
(mrd_x_3, mrd_data_3, time_middle_3) = MRD.completemrd(evaldf3, quantity1, quantity2, M, round(Int, 0.1 * 2^M))
(mrd_time_3, mrd_D_median_3, mrd_D_min_3, mrd_D_max_3, mrd_D_quant1_3, mrd_D_quant3_3) =
    MRD.mrdpp(mrd_x_3, mrd_data_3)
(mrd_x_4, mrd_data_4, time_middle_4) = MRD.completemrd(evaldf4, quantity1, quantity2, M, round(Int, 0.1 * 2^M))
(mrd_time_4, mrd_D_median_4, mrd_D_min_4, mrd_D_max_4, mrd_D_quant1_4, mrd_D_quant3_4) =
    MRD.mrdpp(mrd_x_4, mrd_data_4)
#(mrd_x_5, mrd_data_5, time_middle_5) = MRD.completemrd(evaldf5, quantity1, quantity2, M, round(Int, 0.1 * 2^M))
#(mrd_time_5, mrd_D_median_5, mrd_D_min_5, mrd_D_max_5, mrd_D_quant1_5, mrd_D_quant3_5) =
#    MRD.mrdpp(mrd_x_5, mrd_data_5)
#(mrd_x_6, mrd_data_6, time_middle_6) = MRD.completemrd(evaldf6, quantity1, quantity2, M, round(Int, 0.1 * 2^M))
#(mrd_time_6, mrd_D_median_6, mrd_D_min_6, mrd_D_max_6, mrd_D_quant1_6, mrd_D_quant3_6) =
#    MRD.mrdpp(mrd_x_6, mrd_data_6)

#------------------------------------------------------------------------
#Plot result of MRD
PyPlot.pygui(true)
fig = PyPlot.figure()#figsize=(12,7))
ax = fig.add_subplot(111)
ax.set_title(string("2a MRD ", evaldf2.time[1], " - ", evaldf2.time[end]))
ax.set_xlabel("avg. time [s]")
ax.set_ylabel(L"$C_{wq} [\cdot 10^{-3} \mathrm{ms^{-1}}]$")
ax.grid(true)
PyPlot.xscale("log")
ax.plot(Dates.value.(mrd_time_1) ./ 1000, mrd_D_median_1 .* 1000, label="T1CSAT")
ax.fill_between(Dates.value.(mrd_time_1) ./ 1000, mrd_D_quant1_1 .* 1000, mrd_D_quant3_1 .* 1000, alpha=0.4)
ax.plot(Dates.value.(mrd_time_2) ./ 1000, mrd_D_median_2 .* 1000, label="IRG ice")
ax.fill_between(Dates.value.(mrd_time_2) ./ 1000, mrd_D_quant1_2 .* 1000, mrd_D_quant3_2 .* 1000, alpha=0.4)
ax.plot(Dates.value.(mrd_time_3) ./ 1000, mrd_D_median_3 .* 1000, label="T2CSAT")
ax.fill_between(Dates.value.(mrd_time_3) ./ 1000, mrd_D_quant1_3 .* 1000, mrd_D_quant3_3 .* 1000, alpha=0.4)
ax.plot(Dates.value.(mrd_time_4) ./ 1000, mrd_D_median_4 .* 1000, label="IRG lead")
ax.fill_between(Dates.value.(mrd_time_4) ./ 1000, mrd_D_quant1_4 .* 1000, mrd_D_quant3_4 .* 1000, alpha=0.4)
#ax.plot(Dates.value.(mrd_time_5) ./ 1000, mrd_D_median_5 .* 1000, label="Kaijo")
#ax.fill_between(Dates.value.(mrd_time_5) ./ 1000, mrd_D_quant1_5 .* 1000, mrd_D_quant3_5 .* 1000, alpha=0.4)
#ax.plot(Dates.value.(mrd_time_6) ./ 1000, mrd_D_median_6 .* 1000, label="TJK")
#ax.fill_between(Dates.value.(mrd_time_6) ./ 1000, mrd_D_quant1_6 .* 1000, mrd_D_quant3_6 .* 1000, alpha=0.4)
ax.legend()

######################################################
#MRD time series
#=
tjkdata = copy(evaldf6)

starts = DateTime(2021, 05, 20, 13, 00, 00)
ends = DateTime(2021, 06, 10, 17, 00, 00)

M = 16
mrd_time = zeros(Millisecond, M, length(starts))
mrd_D_median = zeros(Float64, M, length(starts))
mrd_D_min = similar(mrd_D_median)
mrd_D_max = similar(mrd_D_median)
mrd_D_quant1 = similar(mrd_D_median)
mrd_D_quant3 = similar(mrd_D_median)

for idx in 1:length(starts)
    datatouse = tjkdata[starts[idx].<=tjkdata.time.<=ends[idx], :]
    (mrd_x, mrd_data, time_middle) = MRD.completemrd(datatouse, "w", "T", M, round(Int, 0.1 * 2^M))
    (mrd_time[:, idx], mrd_D_median[:, idx], mrd_D_min[:, idx], mrd_D_max[:, idx], mrd_D_quant1[:, idx], mrd_D_quant3[:, idx]) =
        MRD.mrdpp(mrd_x, mrd_data)
end

#=
function animupdate(i::Int64)
    ax.set_title(string("MRD ", starts[i], " - ", timeend))
    pm.set_ydata(mrd_D_median[:, i] .* 1000)
    plq.set_ydata(mrd_D_quant1[:, i] .* 1000)
    puq.set_ydata(mrd_D_quant3[:, i] .* 1000)
    return pm, plq, puq
end
=#

#Plot
for i in 1:length(starts)
    fig = PyPlot.figure()#figsize=(12,7))
    ax = fig.add_subplot(111)
    ax.set_title(string("TJK MRD ", starts[i], " - ", timeend))
    ax.set_xlabel("avg. time [s]")
    ax.set_ylabel(L"$C_{w\theta} [\cdot 10^{-3} \mathrm{Kms^{-1}}]$")
    #ax.set_ylim(-10,10)
    ax.grid(true)
    PyPlot.xscale("log")
    pm, = ax.plot(Dates.value.(mrd_time[:, i]) ./ 1000, mrd_D_median[:, i] .* 1000, label="TJK")
    #puq, = ax.plot(Dates.value.(mrd_time[:, i]) ./ 1000, mrd_D_quant3[:, i] .* 1000, alpha=0.5, color=pm.get_color())
    #plq, = ax.plot(Dates.value.(mrd_time[:, i]) ./ 1000, mrd_D_quant1[:, i] .* 1000, alpha=0.5, color=pm.get_color())
    fl = ax.fill_between(Dates.value.(mrd_time[:, i]) ./ 1000, mrd_D_quant1[:, i] .* 1000, mrd_D_quant3[:, i] .* 1000, alpha=0.4)
    #fig.savefig(string("/home/haugened/Documents/plots/tjk_mrd/", string(starts[i])[1:10], ".png"))
end
#ani = animation.FuncAnimation(fig, animupdate, frames=collect(1:length(starts)), interval=5000, repeat_delay=500)
#ani.save(string("/home/haugened/Desktop/", "short_video_", matfileraw, ".mp4"), fps=framespersec)
##
=#
######################################################
#simple MRD
#=
M=15
evaldf = evaldf3
evaldf = evaldf[(1:2^M),:]
data_a = evaldf.w
data_b = evaldf.T

#orthogonal MRD
(D, ) = MRD.mrd(data_a, data_b, M, 0)
x = 2 .^ collect(1:M) .* Millisecond(50)

#non-orthogonal MRD
(nox, noD) = MRD.fastnomrd(data_a, data_b, 2 .^ collect(M:-1:0))

#plot
fig = PyPlot.figure()#figsize=(12,7))
ax = fig.add_subplot(111)
ax.set_xlabel("time scale [s]")
ax.set_ylabel(L"$C_{wT_s} [10^{-3} \mathrm{Kms^{-1}}]$")
ax.grid(true)
PyPlot.xscale("log")
ax.plot(Dates.value.(x) ./ 1000, D .* 1000, label="orthogonal")
ax.plot(Dates.value.(nox) ./ 1000, noD .* 1000, label="non-orthogonal")
ax.legend()
=#
######################################################