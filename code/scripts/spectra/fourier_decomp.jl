######################################################
###               PERFORMING FOURIER               ###
###                 DECOMPOSITION                  ###
###            author: Michi Haugeneder            ###
######################################################
#=
Performing the necessary preprocessing steps and a Fourier trafo
according to R. Stull, "An Introduction to Boundary Layer Meteorology", 1988, p.307ff
load data with load_data.jl
=#
using LsqFit, FFTW, Statistics, Dates
using PyPlot, PyCall
animation = pyimport("matplotlib.animation")

importdir = joinpath(@__DIR__, "..", "..")
datapath = "/home/haugened/Documents/data/"

include(joinpath(importdir, "src", "ir_evaluation.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "fourier.jl"))
import .irev
import .gen
import .ft

######################################################
###              CHANGE VARIABLES HERE             ###
######################################################
#data to decompose
evaldftmp = copy(evaldf2[1:72000,:])

turb.missing2nan!(evaldftmp)
disallowmissing!(evaldftmp)
ftdata = evaldftmp.w .* evaldftmp.T

#duration in seconds
dur = Dates.value(evaldftmp.time[end] - evaldftmp.time[1])/1e3

######################################################

println()
println("-----------S-T-A-R-T-------------")


######################################################
###                   ALL-IN-ONE                   ###
######################################################

#@time (lp, IRdetrended) = ft.fouriercutoffmatrix(IRdatasmall, dur, fco)

######################################################
###              PREPROCESSING STEPS               ###
######################################################
println("Preprocessing...")
#select vector from dataset
rawvec = ftdata #IRdatasmallraw[500,200,:]
#detrend the vector (subtract linear fit)
vec1 = ft.detrend(rawvec)
#apply bell taper
vect = ft.belltaper(vec1)

######################################################
###               FOURIER TRANSFORM                ###
######################################################
println("Applying Fourier-Trafo...")
(FTraw, FTvec) = ft.fourierplus(vect)
(Soff, freq) = ft.spectralenergydensity(FTvec, dur)

#FT for window
#toprow = [500, 500, 325, 309]#[325, 190, 325, 190, 300, 335]
#leftcol = [230, 40, 780, 656]#[260, 260, 12, 800, 800, 900]

######################################################
###           INVERSE FOURIER TRANSFORM            ###
######################################################
#=
println("Applying inverse Fourier-Trafo...")
#find cutoff frequency
idxfco = ft.findcutofffreq(freq, fco)

#back-trafo
lpift = ft.lowpass(FTraw, idxfco)
#apply inverse bell tapering
lpf = ft.ibelltaper(lpift)
=#

######################################################
###                PLOT THE RESULTS                ###
######################################################

freqtimesSoff = freq.*Soff
(freq, freqtimesSoff) = ft.logavg(freq, freqtimesSoff, 0.05)
##
println("Plotting...")
PyPlot.pygui(true)
fig2 = PyPlot.figure()
ax = PyPlot.gca()
ax.set_title("TJK unstable 03.06. 08:14 - 11:20")
ax.set_xlabel("f [Hz]")
ax.set_ylabel("fS(f)")
tmp = gen.movingaverage(freqtimesSoff, 8)
ax.plot(freq, tmp, label=L"w'T_S'_{TJK}")
#=for idx in 1:length(toprow)
    (freq, freqtimesSoff) = ft.ftforwindow(IRdatasmallraw[(0:9).+toprow[idx],(0:9).+leftcol[idx],:], size(IRdatasmallraw,3)/framespersec, 0.05)
    tmp = gen.movingaverage(freqtimesSoff, 8)
    ax.plot(freq, tmp, label=idx)
end=#
x=collect(0.8:0.05:8)
y=exp.(-2*log.(x)/3)./1
ax.plot(x,y, color="black", label="exp(-2/3)")
PyPlot.xscale("log")
PyPlot.yscale("log")
#ax.set_xlim([0.01, 10])
ax.grid()
ax.legend()
##