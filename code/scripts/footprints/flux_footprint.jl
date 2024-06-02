######################################################
###           CALCULATE FLUX FOOTPRINTS            ###
###            author: Michi Haugeneder            ###
######################################################
#=
Calculate flux footprints according to
N. Kljun et al. (2015) A simple two-dimensional para-
meterisation for Flux Footprint Prediction (FFP), gmd
need to run turb_fluxes before for eg. Obukhov-length
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
LogNorm = pyimport("matplotlib.colors")
mpimg = pyimport("matplotlib.image")

importdir = joinpath(@__DIR__, "..", "..")
datapath = "/home/haugened/Documents/data/"
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "kljun_ffp.jl"))
import .turb
import .gen
import .kljun

#variables
names = [:evaldf1, :evaldf2, :evaldf3, :evaldf4, :evaldf5, :evaldf6]
meas_heights = [1.2, 0.9, 1.9, 2.8, 0.3, 5]
pbl_height = 1000
Ls = [:L1, :L2, :L3, :L4, :L5, :L6]
nrelemgrid = 1000
fluxes = [:fx1, :fx2, :fx3, :fx4, :fx5, :fx6]
rs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] #levels for plotting
winddir = mean(tjkmeteodata.wind_mean_vector_direction)#0 #rotating the footprint
rslayer = false #measurement within roughness sublayer (theory not working properly)
crop = true #crop output to maximum defined rs (max 0.9)
outnames = [:ffp1, :ffp2, :ffp3, :ffp4, :ffp5, :ffp6]

for ix in [2, 4]#1:size(names, 1)
    println("Calculating footprint for ", String(names[ix]))
    turb.missing2nan!(@eval $(names[ix]))
    (x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d, rs, fr, xr, yr, flag_err) = kljun.ffp(
        meas_heights[ix], nothing, 5, pbl_height, mean(filter(!isnan, @eval $(Ls[ix]).L)),
        std(filter(!isnan, @eval $(names[ix]).v)), mean(filter(!isnan, @eval $(fluxes[ix]).u_star)), winddir, rs, rslayer, nrelemgrid, crop)
    if flag_err == true
        @warn("Error flag set to true! There is an error!")
    end
    @eval $(outnames[ix]) = $(Dict("x_ci_max" => x_ci_max, "x_ci" => x_ci, "f_ci" => f_ci,
        "x_2d" => x_2d, "y_2d" => y_2d, "f_2d" => f_2d, "rs" => rs, "fr" => fr,
        "xr" => xr, "yr" => yr, "flag_err" => flag_err))
end

#PyPlot.pygui(true)
#kljun.plot_ffp(ffp6["x_2d"], ffp6["y_2d"], ffp6["f_2d"])
#PyPlot.figure()
#PyPlot.plot(ffp6["xr"], ffp6["yr"])


###############################################
#plotting the footprint on the ortho-mosaic

fileorthomosaic = joinpath(datapath, "pics", "Drone_orthomosaic_20210603_10cm_excerpt.png")
orthomosaic = mpimg.imread(fileorthomosaic)
#PyPlot.imshow(orthomosaic)
#location of flux measurements 1-6 in original image
#[row-location, col-location]
fluxloc = [1406 1507; 1136 1265; 1136 1265; 1136 1265; 1416 1387; 940 1474]

#extend of background [row, col]
#bgextend_m = [280, 280] #in m from measuring in GIS: 279.9
bgextend_pxl = [size(orthomosaic, 1), size(orthomosaic, 2)] #in pxl

#calculate m/pxl from it
meterperpxl_row = 0.1 #bgextend_m[1] / bgextend_pxl[1]
meterperpxl_col = 0.1 #bgextend_m[2] / bgextend_pxl[2]

#origin of figure
figorigin = [1136 1265] #tower 2

#calculate fluxloc in new coordinates [m]
fluxloc_final = Array{Float64}(undef, size(fluxloc, 1), size(fluxloc, 2))
fluxloc_final[:, 1] = (figorigin[1] .- fluxloc[:, 1]) .* meterperpxl_row
fluxloc_final[:, 2] = (fluxloc[:, 2] .- figorigin[2])  .* meterperpxl_col

#calculate extend in new coordinates
lft = (-figorigin[2]) * meterperpxl_col
rght = (bgextend_pxl[2]-1-figorigin[2]) *meterperpxl_col
tp = (bgextend_pxl[1]-(figorigin[1]- bgextend_pxl[1])-1) * meterperpxl_row
btm = (figorigin[1] - bgextend_pxl[1]) * meterperpxl_row

bgextend_final = (-figorigin[2], bgextend_pxl[2]-1-figorigin[2], -(bgextend_pxl[1]-figorigin[1]), bgextend_pxl[1]-(bgextend_pxl[1]-figorigin[1])-1).*meterperpxl_col
#(lft, rght, btm, tp)

##
ctab10 = PyPlot.cm.tab10
ffp_fig = PyPlot.figure()
ax1 = ffp_fig.add_subplot(111)
ax1.set_title("Flux footprints 80% - Kljun et al. (2015)")
bg = ax1.imshow(orthomosaic, extent=bgextend_final)
#bg = ax1.pcolormesh(orthomosaic)
ax1.set_xlabel("meter")
ax1.set_ylabel("meter")
#locfx1 = ax1.plot(fluxloc_final[1, 2], fluxloc_final[1, 1], ".", color=ctab10(0))#, label="T1IRG")
locfx2 = ax1.plot(fluxloc_final[2, 2], fluxloc_final[2, 1], ".", color=ctab10(1))#, label="T2IRG")
#locfx3 = ax1.plot(fluxloc_final[3, 2], fluxloc_final[3, 1], ".", color=ctab10(2))#, label="T2LCSAT")
locfx4 = ax1.plot(fluxloc_final[4, 2], fluxloc_final[4, 1], ".", color=ctab10(3))#, label="T2UCSAT")
#locfx5 = ax1.plot(fluxloc_final[5, 2], fluxloc_final[5, 1], ".", color=ctab10(4))#, label="Kaijo")
#locfx6 = ax1.plot(fluxloc_final[6, 2], fluxloc_final[6, 1], ".", color=ctab10(5))#, label="TJK")
#fp1 = ax1.plot(ffp1["xr"][:, end-1] .+ fluxloc_final[1, 2], ffp1["yr"][:, end-1] .+ fluxloc_final[1, 1], color=ctab10(0), label = "T1IRG")
fp2 = ax1.plot(ffp2["xr"][:, end-1] .+ fluxloc_final[2, 2], ffp2["yr"][:, end-1] .+ fluxloc_final[2, 1], color=ctab10(1), label = "T2IRG")
#fp3 = ax1.plot(ffp3["xr"][:, end-1] .+ fluxloc_final[3, 2], ffp3["yr"][:, end-1] .+ fluxloc_final[3, 1], color=ctab10(2), label = "T2LCSAT")
fp4 = ax1.plot(ffp4["xr"][:, end-1] .+ fluxloc_final[4, 2], ffp4["yr"][:, end-1] .+ fluxloc_final[4, 1], color=ctab10(3), label = "T2UCSAT")
#fp5 = ax1.plot(ffp5["xr"][:, end-1] .+ fluxloc_final[5, 2], ffp5["yr"][:, end-1] .+ fluxloc_final[5, 1], color=ctab10(4), label = "Kaijo")
#fp6 = ax1.plot(ffp6["xr"][:, end-1] .+ fluxloc_final[6, 2], ffp6["yr"][:, end-1] .+ fluxloc_final[6, 1], color=ctab10(5), label = "TJK")
ax1.legend()
##
println("------------D-O-N-E---------------")
println()