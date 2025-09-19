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
#datapath = "/home/haugened/Documents/data/"
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
quantity1 = "w"
quantity2 = "T"

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
ax.set_title(string("1a MRD ", evaldf2.time[1], " - ", evaldf2.time[end]))
ax.set_xlabel("avg. time [s]")
ax.set_ylabel(L"$C_{wT} [\cdot 10^{-3} \mathrm{Kms^{-1}}]$")
ax.grid(true)
PyPlot.xscale("log")
ax.plot(Dates.value.(mrd_time_1) ./ 1000, mrd_D_median_1 .* 1000, label="T1 IRG")
ax.fill_between(Dates.value.(mrd_time_1) ./ 1000, mrd_D_quant1_1 .* 1000, mrd_D_quant3_1 .* 1000, alpha=0.4)
ax.plot(Dates.value.(mrd_time_2) ./ 1000, mrd_D_median_2 .* 1000, label="T1 CSAT")
ax.fill_between(Dates.value.(mrd_time_2) ./ 1000, mrd_D_quant1_2 .* 1000, mrd_D_quant3_2 .* 1000, alpha=0.4)
ax.plot(Dates.value.(mrd_time_3) ./ 1000, mrd_D_median_3 .* 1000, label="T2 IRG")
ax.fill_between(Dates.value.(mrd_time_3) ./ 1000, mrd_D_quant1_3 .* 1000, mrd_D_quant3_3 .* 1000, alpha=0.4)
ax.plot(Dates.value.(mrd_time_4) ./ 1000, mrd_D_median_4 .* 1000, label="T2 CSAT")
ax.fill_between(Dates.value.(mrd_time_4) ./ 1000, mrd_D_quant1_4 .* 1000, mrd_D_quant3_4 .* 1000, alpha=0.4)
#ax.plot(Dates.value.(mrd_time_5) ./ 1000, mrd_D_median_5 .* 1000, label="Kaijo")
#ax.fill_between(Dates.value.(mrd_time_5) ./ 1000, mrd_D_quant1_5 .* 1000, mrd_D_quant3_5 .* 1000, alpha=0.4)
#ax.plot(Dates.value.(mrd_time_6) ./ 1000, mrd_D_median_6 .* 1000, label="TJK")
#ax.fill_between(Dates.value.(mrd_time_6) ./ 1000, mrd_D_quant1_6 .* 1000, mrd_D_quant3_6 .* 1000, alpha=0.4)
ax.legend()
######################################################
##
#Plot result in a 4-panel plot
PyPlot.pygui(true)
fig, axes = PyPlot.subplots(2, 2, figsize=(12, 8))
fig.suptitle(string("2a MRD ", evaldf2.time[1], " - ", evaldf2.time[end]))

# Common axis settings function
function setup_axis(ax, title_text)
    ax.set_title(title_text)
    ax.set_xlabel("avg. time [s]")
    ax.set_ylabel(L"$C_{wT} [\cdot 10^{-3} \mathrm{Kms^{-1}}]$")
    #ax.set_ylabel(L"$C_{uw} [\cdot 10^{-3} \mathrm{m^2s^{-2}}]$")
    #ax.set_ylabel(L"$C_{uv} [\cdot 10^{-3} \mathrm{m^2s^{-2}}]$")
    #ax.set_ylabel(L"$C_{wq} [\cdot 10^{-3} \mathrm{ms^{-1}gm^{-3}}]$")
    ax.grid(true)
    ax.set_xscale("log")
    ax.set_ylim(-5, 5)
end

# Top left - T1 IRG
setup_axis(axes[1,1], "Ice 1.1m")
axes[1,1].plot(Dates.value.(mrd_time_1) ./ 1000, mrd_D_median_1 .* 1000, label="T1 IRG")
axes[1,1].fill_between(Dates.value.(mrd_time_1) ./ 1000, mrd_D_quant1_1 .* 1000, mrd_D_quant3_1 .* 1000, alpha=0.4)

# Top right - T2 IRG
setup_axis(axes[1,2], "Lead 1.3m")
axes[1,2].plot(Dates.value.(mrd_time_3) ./ 1000, mrd_D_median_3 .* 1000, label="T2 IRG")
axes[1,2].fill_between(Dates.value.(mrd_time_3) ./ 1000, mrd_D_quant1_3 .* 1000, mrd_D_quant3_3 .* 1000, alpha=0.4)

# Bottom left - T1 CSAT
setup_axis(axes[2,1], "Ice 2.1m")
axes[2,1].plot(Dates.value.(mrd_time_2) ./ 1000, mrd_D_median_2 .* 1000, label="T1 CSAT")
axes[2,1].fill_between(Dates.value.(mrd_time_2) ./ 1000, mrd_D_quant1_2 .* 1000, mrd_D_quant3_2 .* 1000, alpha=0.4)

# Bottom right - T2 CSAT
setup_axis(axes[2,2], "Lead 2.3m")
axes[2,2].plot(Dates.value.(mrd_time_4) ./ 1000, mrd_D_median_4 .* 1000, label="T2 CSAT")
axes[2,2].fill_between(Dates.value.(mrd_time_4) ./ 1000, mrd_D_quant1_4 .* 1000, mrd_D_quant3_4 .* 1000, alpha=0.4)

# Adjust layout to prevent overlap
PyPlot.tight_layout()
##
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