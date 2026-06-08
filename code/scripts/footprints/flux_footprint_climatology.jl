######################################################
### CALCULATE FLUX FOOTPRINTS CLIMATOLOGY (PYTHON) ###
###            author: Michi Haugeneder            ###
######################################################
#=
Calculate (climatological) flux footprints according to
N. Kljun et al. (2015) A simple two-dimensional para-
meterisation for Flux Footprint Prediction (FFP), gmd
need to run turb_fluxes before for eg. Obukhov-length
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter
using Images
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
LogNorm = pyimport("matplotlib.colors")
mpimg = pyimport("matplotlib.image")

importdir = joinpath(@__DIR__, "..", "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "kljun_ffp.jl"))
include(joinpath(importdir, "src", "footprint_plotting.jl"))
if !@isdefined stationcfg
    include(joinpath(importdir, "src", "station_config.jl"))
    import .stationcfg
end
import .turb
import .gen
import .kljun
@pyinclude(joinpath(importdir, "src", "kljun_ffp_climatology.py"))
#PyPlot.pygui(true)

if !@isdefined station_config
    error("Run load_data.jl before flux_footprint_climatology.jl so station_config is available.")
end
station_label = stationcfg.station_label(station_config)

#variables
names = [:evaldf1, :evaldf2, :evaldf3, :evaldf4]
meas_heights = Float64.(stationcfg.optional_key(station_config, heights, "footprint", "measurement_heights"))
pbl_height = 1000.0
fluxes = [:fx1, :fx2, :fx3, :fx4]
wd = [:wd1, :wd2, :wd3, :wd4] #wind directions
outnames = [:ffp1, :ffp2, :ffp3, :ffp4]
aggtime = Minute(30) #aggregation time
wind_direction_offsets = Float64.(stationcfg.optional_key(
    station_config,
    fill(0.0, length(names)),
    "footprint",
    "wind_direction_offsets",
))
contour_indices = Int.(stationcfg.optional_key(
    station_config,
    fill(-1, length(names)),
    "footprint",
    "contour_indices",
))

function footprint_contour(contours, end_offset::Integer)
    contour_index = lastindex(contours) + end_offset
    firstindex(contours) <= contour_index <= lastindex(contours) || error(
        "Footprint contour offset $end_offset is outside available contour range."
    )
    return contours[contour_index]
end

#optional input
domain = [nothing]#[-1e6, 1e6, -1e6, 1e6]
dx = nothing
dy = nothing
nx = nothing
ny = nothing
rs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]#, 0.9] #levels for plotting
rslayer = 0 #measurement within roughness sublayer (theory not working properly)
smooth_data = 1
crop = false #crop output to maximum defined rs (max 0.9)
pulse = nothing
verbosity = 2
fig = false

output = nothing

for ix in 1:size(names, 1)
    println("Calculating footprint for ", String(names[ix]))

    ecdata = @eval $(names[ix])
    turb.missing2nan!(ecdata)
    fluxdata = @eval $(fluxes[ix])
    turb.missing2nan!(fluxdata)
    obukl = fluxdata.L_highfreq
    wd_tmp = @eval $(wd[ix])

    #aggregate data to half-hour intervals
    duration = ecdata.time[end]-ecdata.time[1]
    nrblocks = div(duration, aggtime)
    aggidcs = round(Int, aggtime/Millisecond(50))
    println("Evaluating ", nrblocks, " blocks, last ", mod(duration, aggtime), " discarded.")
    umean = fill(NaN, nrblocks)
    h = fill(pbl_height, nrblocks)
    ol = fill(NaN, nrblocks)
    sigmav = fill(NaN, nrblocks)
    ustar = fill(NaN, nrblocks)
    wind_dir = fill(NaN, nrblocks)

    for j in 1:nrblocks
        six = 1 + (j - 1)*aggidcs
        eix = six + (aggidcs - 1)
        umean[j] = mean(filter(!isnan, sqrt.(ecdata.u[six:eix] .^2 .+ ecdata.v[six:eix] .^2 .+ ecdata.w[six:eix] .^2)))
        ol[j] = mean(filter(!isnan, obukl[six:eix]))
        sigmav[j] = std(filter(!isnan, ecdata.v[six:eix]))
        ustar[j] = mean(filter(!isnan, fluxdata.u_star[six:eix]))
        wind_dir_raw = filter(!isnan, wd_tmp[ecdata.time[six] .<= wd_tmp.time .< ecdata.time[eix], :α])
        wind_dir[j] = turb.mean_winddir(wind_dir_raw)
        wind_dir[j] = (wind_dir[j] + wind_direction_offsets[ix]) % 360
    end

    output = py"FFP_climatology"(meas_heights[ix], nothing, PyVector(umean), PyVector(h), PyVector(ol),
    PyVector(sigmav), PyVector(ustar), PyVector(wind_dir), PyVector(domain), dx, dy, nx, ny, PyVector(rs), rslayer, smooth_data, crop, pulse, verbosity, fig)
    if output["flag_err"] != 0
        @warn("Error flag set to true! There is an error!")
    end
    @eval $(outnames[ix]) = output
end

###############################################
#plotting the footprint on the ortho-mosaic

fileorthomosaic = String(stationcfg.require_key(station_config, "footprint", "orthomosaic"))
orthomosaic = mpimg.imread(fileorthomosaic)
orthomosaic_jl = load(fileorthomosaic)
#PyPlot.imshow(orthomosaic)
#location of flux measurements 1-6 in original image
#[row-location, col-location]
fluxloc = stationcfg.toml_matrix(stationcfg.require_key(station_config, "footprint", "fluxloc"); T=Float64)

#fallback extent of background [row, col], used only when no world file exists
bgextend_m_config = stationcfg.optional_key(station_config, nothing, "footprint", "bgextend_m")
bgextend_pxl_config = stationcfg.optional_key(station_config, nothing, "footprint", "bgextend_pxl")
bgextend_m = isnothing(bgextend_m_config) ? nothing : Float64.(bgextend_m_config)
bgextend_pxl = isnothing(bgextend_pxl_config) ? nothing : Float64.(bgextend_pxl_config)

#origin of figure
figorigin = Float64.(stationcfg.require_key(station_config, "footprint", "figorigin"))

fluxloc_final, bgextend_final = footprint_background_geometry(
    fluxloc,
    bgextend_m,
    bgextend_pxl,
    figorigin;
    image_file=fileorthomosaic,
    image_size=size(orthomosaic_jl),
)

##
ctab10 = PyPlot.cm.tab10
ffp_fig = PyPlot.figure(figsize=(10,8))
ax1 = ffp_fig.add_subplot(111)
ax1.set_title("Station $(station_label) Flux footprints")
bg = ax1.imshow(orthomosaic, extent=bgextend_final)
#bg = ax1.pcolormesh(orthomosaic)
ax1.set_xlabel("meter")
ax1.set_ylabel("meter")
locfx1 = ax1.plot(fluxloc_final[1, 2], fluxloc_final[1, 1], ".", color=ctab10(0), ms=15)
#locfx2 = ax1.plot(fluxloc_final[2, 2], fluxloc_final[2, 1], ".", color=ctab10(1), ms=15)
locfx3 = ax1.plot(fluxloc_final[3, 2], fluxloc_final[3, 1], ".", color=ctab10(2), ms=15)
#locfx4 = ax1.plot(fluxloc_final[4, 2], fluxloc_final[4, 1], ".", color=ctab10(3), ms=15)
fp1 = ax1.plot(footprint_contour(ffp1["xr"], contour_indices[1]) .+ fluxloc_final[1, 2], footprint_contour(ffp1["yr"], contour_indices[1]) .+ fluxloc_final[1, 1], color=ctab10(0), label = instr_labels[1])
fp2 = ax1.plot(footprint_contour(ffp2["xr"], contour_indices[2]) .+ fluxloc_final[2, 2], footprint_contour(ffp2["yr"], contour_indices[2]) .+ fluxloc_final[2, 1], color=ctab10(1), label = instr_labels[2])
fp3 = ax1.plot(footprint_contour(ffp3["xr"], contour_indices[3]) .+ fluxloc_final[3, 2], footprint_contour(ffp3["yr"], contour_indices[3]) .+ fluxloc_final[3, 1], color=ctab10(2), label = instr_labels[3])
fp4 = ax1.plot(footprint_contour(ffp4["xr"], contour_indices[4]) .+ fluxloc_final[4, 2], footprint_contour(ffp4["yr"], contour_indices[4]) .+ fluxloc_final[4, 1], color=ctab10(3), label = instr_labels[4])
ax1.legend()
PyPlot.tight_layout()
##
println("------------D-O-N-E---------------")
println()
