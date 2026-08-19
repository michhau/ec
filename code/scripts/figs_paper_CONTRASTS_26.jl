######################################################
###         FIGURES FOR PAPER CONTRASTS 2026       ###
###            author: Michi Haugeneder            ###
######################################################

###################################################################################################
#Figures for lead section (Figs. xxx)

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
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter, Distributed, Printf, Logging
using NCDatasets
using Images
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
cramericm = pyimport("cmcrameri.cm")
LogNorm = pyimport("matplotlib.colors")
mpimg = pyimport("matplotlib.image")
mpl_scatter_density = pyimport("mpl_scatter_density")
patches = pyimport("matplotlib.patches")
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
animation = pyimport("matplotlib.animation")

station_name = "1c" # Change this to "1a", "1b_1", ..., "3d" before including the script.

importdir = joinpath(@__DIR__, "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "station_config.jl"))
include(joinpath(importdir, "src", "kljun_ffp.jl"))
include(joinpath(importdir, "src", "footprint_plotting.jl"))
include(joinpath(importdir, "src", "ffp_block.jl"))
include(joinpath(importdir, "src", "block_analysis_src.jl"))
import .turb
import .gen
import .stationcfg
import .kljun
import .ffp_block
import .block_analyze
@pyinclude(joinpath(importdir, "src", "kljun_ffp_climatology.py"))

timestep = Millisecond(50)
ρ_air = 1.2 #kg m^{-3}
c_p = 1004 #J kg^{-1} K^{-1}
L_v = 2500e3 #J kg^{-1} (approx @0°C)

######################################################
#variables

station_config = stationcfg.load_station_config(station_name)
station_label = stationcfg.station_label(station_config)
station_file_stem = stationcfg.station_file_stem(station_config)

datapath = String(stationcfg.require_key(station_config, "data_root"))
#datapath = "/home/michi/Documents/slf/CONTRASTS25/data/processed/preproc/"

#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
nanfile_stam = String(stationcfg.optional_key(station_config, "/home/haugened/Documents/data/CONTRASTS/nan_periods", "manual_nanmask", "root"))
nanfiles = joinpath.(nanfile_stam, String.(stationcfg.optional_key(
    station_config,
    ["t1irg_nan.csv", "t1csat_nan.csv", "t2irg_nan.csv", "t2csat_nan.csv"],
    "manual_nanmask",
    "files",
)))
plot_output_root = stationcfg.plot_root(station_config)

######################################################

println()
println("-----------S-T-A-R-T-------------")
println("Loading station ", station_label, " from ", station_config["config_file"])

######################################################
###            LOADING & PREPROCESSING             ###
######################################################
#select data and measurement period to be evaluated
evalstart = stationcfg.station_datetime(station_config, "evalstart")
evalend   = stationcfg.station_datetime(station_config, "evalend")
#evalend = evalstart + Day(10)

data_files = String.(stationcfg.require_key(station_config, "data_files"))

evaldfs = [
    turb.readturbasnetcdf(joinpath(datapath, data_file), evalstart, evalend)
    for data_file in data_files
]
(evaldf1, evaldf2, evaldf3, evaldf4) = evaldfs

# Instrument metadata for labeling (corresponds to evaldf1-4 / fx1-4)
srf_type = String.(stationcfg.require_key(station_config, "surface_type"))
heights = Float64.(stationcfg.require_key(station_config, "heights"))
instr_type = String.(stationcfg.optional_key(station_config, ["IRG", "CSAT", "IRG", "CSAT"], "instrument_type"))
instr_labels = stationcfg.station_labels(station_config)

#show statistics about missing or NaN-data
for evaldf in evaldfs
    turb.printmissstats(evaldf)
end

#interpolate missing
evaldfs = [turb.interpolatemissing(evaldf) for evaldf in evaldfs]
(evaldf1, evaldf2, evaldf3, evaldf4) = evaldfs

wds = [turb.winddir(evaldf) for evaldf in evaldfs]

#double rotation
for evaldf in evaldfs
    turb.drdf!(evaldf, periodwise=false)
end

#CONTRASTS: apply manual nan-mask
if stationcfg.optional_key(station_config, true, "manual_nanmask", "enabled")
    nanmask_columns = String.(stationcfg.optional_key(station_config, ["u", "v", "w", "T"], "manual_nanmask", "columns"))
    manual_periods = [turb.read_nantimes_csv(nanfile) for nanfile in nanfiles]
    for (evaldf, manual_period) in zip(evaldfs, manual_periods)
        isbad_quality = turb.manual_nanmask(manual_period, evaldf.time)
        evaldf[isbad_quality, nanmask_columns] .= NaN
    end
end

#tower rotation correction for station 3d
stationcfg.apply_station_wind_direction_rotations!(wds, station_config)
(wd1, wd2, wd3, wd4) = wds

#Fluxes (from turb_fluxes.jl)

######################################################
###              TURBULENT FLUXES                  ###
######################################################
turb.missing2nan!(evaldf1)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf3)
turb.missing2nan!(evaldf4)

#Reynolds averaging times (determine usind MRD)
ra1 = Second(400) #Millisecond(2^11*timestep) #Second(330)
ra2 = Second(400)
ra3 = Second(400) #Millisecond(2^11*timestep) #Second(320)
ra4 = Second(400)

fx1_raw = turb.turbflux(evaldf1, ra1, 1.0)
fx2_raw = turb.turbflux(evaldf2, ra2, 1.0)
fx3_raw = turb.turbflux(evaldf3, ra3, 1.0)
fx4_raw = turb.turbflux(evaldf4, ra4, 1.0)

#averaging
fx1 = turb.avgflux(fx1_raw, ra1, true, 0.1)
fx2 = turb.avgflux(fx2_raw, ra2, true, 0.1)
fx3 = turb.avgflux(fx3_raw, ra3, true, 0.1)
fx4 = turb.avgflux(fx4_raw, ra4, true, 0.1)

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
##

#variables
names = [:evaldf1, :evaldf2, :evaldf3, :evaldf4]
meas_heights = Float64.(stationcfg.optional_key(station_config, heights, "footprint", "measurement_heights"))
pbl_height = 200.0
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

    #keep only periods for which every required footprint input is finite
    valid_blocks = map(eachindex(umean)) do j
        all(isfinite, (umean[j], ol[j], sigmav[j],ustar[j], wind_dir[j],))
    end

    nvalid = sum(valid_blocks)

    if nvalid == 0
        error("No valid periods available for $(names[ix]).")
    end

    println("Using $nvalid of $nrblocks periods; ", nrblocks-nvalid, " period(s) containing NaN/Inf were ignored.")

    output = py"FFP_climatology"(meas_heights[ix], nothing, PyVector(umean[valid_blocks]), PyVector(h[valid_blocks]), PyVector(ol[valid_blocks]),
    PyVector(sigmav[valid_blocks]), PyVector(ustar[valid_blocks]), PyVector(wind_dir[valid_blocks]), PyVector(domain), dx, dy, nx, ny, PyVector(rs), rslayer, smooth_data, crop, pulse, verbosity, fig)
    if output["flag_err"] != 0
        @warn("Error flag set to true! There is an error!")
    end
    @eval $(outnames[ix]) = output
end

###############################################
#plotting the footprint on the ortho-mosaic

fileorthomosaic = String(stationcfg.require_key(station_config, "footprint", "orthomosaic"))
orthomosaic = mpimg.imread(fileorthomosaic);
orthomosaic_jl = load(fileorthomosaic);

#=
PyPlot.pygui(true)
PyPlot.imshow(orthomosaic)
PyPlot.pygui(false)
=#
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

PyPlot.pygui(true)
ctab10 = PyPlot.cm.tab10
ffp_fig = PyPlot.figure()
ax1 = ffp_fig.add_subplot(111)
#ax1.set_title("Station $(station_label) Flux footprints")
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
ax1.legend()#loc="lower right")
ax1.set_xlim(-20,130)
ax1.set_ylim(-77,60)
PyPlot.tight_layout()
#PyPlot.savefig(joinpath("/home/haugened/Documents/data/CONTRASTS/plots/paper_CONTRASTS_26/", "1c_ffp.pdf"), bbox_inches="tight")

##

######################################################
# Flux time series with dual subplots for latent and sensible heat fluxes (CONTRASTS)
##
fig, (ax1, ax2) = PyPlot.subplots(2, 1, figsize=(10, 8), sharex=true)

# Define colors for consistent legend
colors = ["C0", "C1", "C2", "C3"]  # Default matplotlib color cycle

#define plotting step size
step = 20*60 #every 1min

#y-axis limits
wT_limits = Tuple(Float64.(stationcfg.optional_key(station_config, [-20.0, 20.0], "plot", "wT_limits")))
wq_limits = Tuple(Float64.(stationcfg.optional_key(station_config, [-5.0, 5.0], "plot", "wq_limits")))

# Upper subplot - Buoyancy fluxes (sensible heat)
wt1 = ax1.plot(fx1.time[1:step:end], fx1.wT[1:step:end] .* (ρ_air * c_p), color=colors[1])
wt2 = ax1.plot(fx2.time[1:step:end], fx2.wT[1:step:end] .* (ρ_air * c_p), color=colors[2])
wt3 = ax1.plot(fx3.time[1:step:end], fx3.wT[1:step:end] .* (ρ_air * c_p), color=colors[3])
wt4 = ax1.plot(fx4.time[1:step:end], fx4.wT[1:step:end] .* (ρ_air * c_p), color=colors[4])
ax1.set_ylabel(L"Q_H ~\mathrm{[W~m^{-2}]}")
ax1.grid()
ax1.set_ylim(wT_limits)

# Lower subplot - Latent heat fluxes
wq1 = ax2.plot(fx1.time[1:step:end], fx1.wq[1:step:end] .* (L_v * 1e-3), color=colors[1])
wq2 = ax2.plot(fx3.time[1:step:end], fx3.wq[1:step:end] .* (L_v * 1e-3), color=colors[3])
ax2.set_ylabel(L"Q_E ~\mathrm{[W~m^{-2}]}")
ax2.set_xlabel("Time")
ax2.xaxis_date()
ax2.grid()
ax2.set_ylim(wq_limits)

# Create a single legend for the entire figure
handles = [wt1[1], wt2[1], wt3[1], wt4[1]]  # Get line objects
labels = instr_labels
ax1.legend(handles, labels)#, loc="upper right", bbox_to_anchor=(1.0, 1))
ax2.legend([wq1[1], wq2[1]], [instr_labels[1], instr_labels[3]])

# Optional: Uncomment these lines if you want to set specific time limits
# ax2.set_xlim(DateTime(2021, 05, 31, 10, 30, 00), DateTime(2021, 05, 31, 20, 00))

# Format dates on x-axis if needed
# majorlocator = pydates.HourLocator(interval=1)
# minorlocator = pydates.MinuteLocator(interval=15)
# ax2.xaxis.set_major_locator(majorlocator)
# ax2.xaxis.set_minor_locator(minorlocator)
# date_format = pydates.DateFormatter("%H:%M")
# ax2.xaxis.set_major_formatter(date_format)
# fig.autofmt_xdate()

PyPlot.tight_layout()
output_folder = "/home/haugened/Documents/data/CONTRASTS/plots/paper_CONTRASTS_26/"
#PyPlot.savefig(joinpath(output_folder, "1c_heat_fluxes.pdf"), bbox_inches="tight")

#depending on the data length, the following block takes more than 10min!
#when no parameters are changed, the ready file is read afterwards automatically
#only uncomment, if you change parameters
#=
##
######################################################
###      CALCULATE FOOTPRINTS PER BLOCK FLUX       ###
###            author: Michi Haugeneder            ###
######################################################

#variables
reyavg_periods = [:ra1, :ra2, :ra3, :ra4]
block_flux_names = [:block_fluxes1, :block_fluxes2, :block_fluxes3, :block_fluxes4]
footprint_input_names = [:ffp_inputs1, :ffp_inputs2, :ffp_inputs3, :ffp_inputs4]
outnames_block = [:ffp_blocks1, :ffp_blocks2, :ffp_blocks3, :ffp_blocks4]
nrelemgrid = 1000

block_fluxes_all = Dict{Symbol, DataFrame}()
ffp_inputs_all = Dict{Symbol, DataFrame}()
ffp_blocks_all = Dict{Symbol, Vector{Dict{String, Any}}}()

#optional input
rs_block = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] #levels for later plotting/analysis
rslayer_block = false

println()
println("Calculating block fluxes and per-block footprints for station ", station_label)

PyPlot.pygui(false)
for ix in eachindex(names)
    println("Processing ", String(names[ix]))

    ecdata = @eval $(names[ix])
    fluxdata = @eval $(fluxes[ix])
    peri = @eval $(reyavg_periods[ix])
    wd_tmp = @eval $(wd[ix])

    turb.missing2nan!(ecdata)
    turb.missing2nan!(fluxdata)

    block_flux, center_indices, start_indices, end_indices = ffp_block.extract_block_fluxes(fluxdata, peri)
    println("Extracted ", nrow(block_flux), " full non-overlapping blocks from ", String(fluxes[ix]), ".")

    inputs = ffp_block.block_footprint_inputs(
        ecdata,
        block_flux,
        wd_tmp,
        start_indices,
        end_indices,
        meas_heights[ix],
        pbl_height,
        wind_direction_offsets[ix],
    )

    block_flux[!, :u] = [
        ffp_block.nanmean(
            sqrt.(ecdata.u[start_indices[j]:end_indices[j]] .^ 2 .+
                  ecdata.v[start_indices[j]:end_indices[j]] .^ 2)
        )
        for j in eachindex(start_indices)
    ]
    block_flux[!, :T] = [
        ffp_block.nanmean(ecdata.T[start_indices[j]:end_indices[j]])
        for j in eachindex(start_indices)
    ]
    block_flux[!, :wind_dir] = [
        turb.mean_winddir(wd_tmp[ecdata.time[start_indices[j]] .<= wd_tmp.time .<
            ecdata.time[end_indices[j]], Symbol("α")])
        for j in eachindex(start_indices)
    ]

    footprints = ffp_block.calculate_footprints(inputs, meas_heights[ix], rs_block, rslayer_block, nrelemgrid, crop)

    block_fluxes_all[fluxes[ix]] = block_flux
    ffp_inputs_all[fluxes[ix]] = inputs
    ffp_blocks_all[fluxes[ix]] = footprints

    @eval $(block_flux_names[ix]) = $block_flux
    @eval $(footprint_input_names[ix]) = $inputs
    @eval $(outnames_block[ix]) = $footprints
end
PyPlot.pygui(true)

##################################################
# Combine with a drone overview image to classify
using Images
classes_file = String(station_config["block_footprints"]["classes_file"])
classes_reference_file = String(stationcfg.require_key(station_config, "footprint", "orthomosaic"))
classes_scale_file = isfile(ffp_block.footprint_world_file(classes_file)) ? classes_file : classes_reference_file
classes_img = load(classes_file);
classes_img_size = size(classes_img);
#for plotting
mpimg = pyimport("matplotlib.image")
classes_img_pyplot = mpimg.imread(classes_file)

"""
    read_criterion_from_config(config)

Read RGB-criterions for ice and feature (e.g. lead, pond, ridge) from the config Dictionary
"""
function read_criterion_from_config(station_config::Dict)
    #read criterion for true classification = feature (e.g. lead, pond, ridge) from config file
    ice_red_crit_bigger_than = Int(station_config["block_footprints"]["crit_ice_red_bigger_eq"])
    ice_red_crit_smaller_eq = Int(station_config["block_footprints"]["crit_ice_red_smaller_eq"])
    ice_green_crit_bigger_than = Int(station_config["block_footprints"]["crit_ice_green_bigger_eq"])
    ice_green_crit_smaller_eq = Int(station_config["block_footprints"]["crit_ice_green_smaller_eq"])
    ice_blue_crit_bigger_than = Int(station_config["block_footprints"]["crit_ice_blue_bigger_eq"])
    ice_blue_crit_smaller_eq = Int(station_config["block_footprints"]["crit_ice_blue_smaller_eq"])

    ice_crit_rgb = [[ice_red_crit_bigger_than, ice_red_crit_smaller_eq],
                    [ice_green_crit_bigger_than, ice_green_crit_smaller_eq],
                    [ice_blue_crit_bigger_than, ice_blue_crit_smaller_eq]]

    feature_red_crit_bigger_than = Int(station_config["block_footprints"]["crit_feature_red_bigger_eq"])
    feature_red_crit_smaller_eq = Int(station_config["block_footprints"]["crit_feature_red_smaller_eq"])
    feature_green_crit_bigger_than = Int(station_config["block_footprints"]["crit_feature_green_bigger_eq"])
    feature_green_crit_smaller_eq = Int(station_config["block_footprints"]["crit_feature_green_smaller_eq"])
    feature_blue_crit_bigger_than = Int(station_config["block_footprints"]["crit_feature_blue_bigger_eq"])
    feature_blue_crit_smaller_eq = Int(station_config["block_footprints"]["crit_feature_blue_smaller_eq"])

    feature_crit_rgb = [[feature_red_crit_bigger_than, feature_red_crit_smaller_eq],
                            [feature_green_crit_bigger_than, feature_green_crit_smaller_eq],
                            [feature_blue_crit_bigger_than, feature_blue_crit_smaller_eq]]

    return ice_crit_rgb, feature_crit_rgb
end

#create classes matrix
"""
    create_classes_array(img::Matrix{RGB{N0f8}}, )

Create a class array with surface classes from an input image.
Criterions need to be set in the config file!
"""
function create_classes_array(img::Matrix{RGB{N0f8}},
    ice_crit::Vector, feature_crit::Vector)::Array{Union{Missing, Bool}}
    classes_array = Array{Union{Missing, Bool}}(undef, size(img))#fill(false, size(img))
    #missing: do not use
    #true: ice
    #false: feature such as lead, pond, ridge
    red_vals = round.(Int, float.(red.(img))*255)
    green_vals = round.(Int, float.(green.(img))*255)
    blue_vals = round.(Int, float.(blue.(img))*255)
    for i in eachindex(classes_array)
        if ice_crit[1][1] <= red_vals[i] <= ice_crit[1][2] &&
           ice_crit[2][1] <= green_vals[i] <= ice_crit[2][2] &&
           ice_crit[3][1] <= blue_vals[i] <= ice_crit[3][2]
            classes_array[i] = true
        elseif feature_crit[1][1] <= red_vals[i] <= feature_crit[1][2] &&
               feature_crit[2][1] <= green_vals[i] <= feature_crit[2][2] &&
               feature_crit[3][1] <= blue_vals[i] <= feature_crit[3][2]
            classes_array[i] = false
        else
            classes_array[i] = missing
        end
    end
    return classes_array
end

ice_crit_rgb, feature_crit_rgb = read_criterion_from_config(station_config)
surface_feature = String(station_config["block_footprints"]["feature"])

#read same as for footprints. Change, if pictures are not the same!
classes_array = create_classes_array(classes_img, ice_crit_rgb, feature_crit_rgb)

classes_fluxloc_pxl = stationcfg.toml_matrix(stationcfg.require_key(station_config, "footprint", "fluxloc"); T=Float64)
classes_extend_m_config = stationcfg.optional_key(station_config, nothing, "footprint", "bgextend_m")
classes_extend_pxl_config = stationcfg.optional_key(station_config, nothing, "footprint", "bgextend_pxl")
classes_extend_m = isnothing(classes_extend_m_config) ? nothing : Float64.(classes_extend_m_config)
classes_extend_pxl = isnothing(classes_extend_pxl_config) ? nothing : Float64.(classes_extend_pxl_config)
classes_figorigin_pxl = Float64.(stationcfg.require_key(station_config, "footprint", "figorigin"))

classes_scale = ffp_block.footprint_meterperpxl(classes_scale_file, classes_extend_m, classes_extend_pxl)
classes_meterperpxl_row = classes_scale.meterperpxl_row
classes_meterperpxl_col = classes_scale.meterperpxl_col
classes_coordinates_zero_based = true

"""
    footprint_weighted_class_fractions(footprint, classes_array, origin_pxl,
        meterperpxl_row, meterperpxl_col; origin_is_zero_based=true)

Calculate the footprint-weighted ice and feature (lead/pond/ridge) fractions for one footprint.
Footprint points outside `classes_array` or on `missing` class pixels are
excluded from the ice/feature denominator.
"""
function footprint_weighted_class_fractions(footprint::AbstractDict,
        classes_array::AbstractMatrix{Union{Missing, Bool}},
        origin_pxl::AbstractVector{<:Real},
        meterperpxl_row::Real,
        meterperpxl_col::Real;
        origin_is_zero_based::Bool=true)

    if get(footprint, "flag_err", false)
        return (
            ice_fraction=missing,
            feature_fraction=missing,
            classified_weight_fraction=0.0,
            outside_weight_fraction=missing,
            missing_class_weight_fraction=missing,
            ice_weight=0.0,
            feature_weight=0.0,
            classified_weight=0.0,
            total_weight=0.0,
        )
    end

    x = footprint["x_2d"]
    y = footprint["y_2d"]
    f = footprint["f_2d"]

    size(x) == size(y) == size(f) || error("Footprint x_2d, y_2d, and f_2d must have the same size.")

    row0 = Float64(origin_pxl[1])
    col0 = Float64(origin_pxl[2])
    index_offset = origin_is_zero_based ? 1 : 0
    nrows, ncols = size(classes_array)

    ice_weight = 0.0
    feature_weight = 0.0
    classified_weight = 0.0
    outside_weight = 0.0
    missing_class_weight = 0.0
    total_weight = 0.0

    @inbounds for k in eachindex(f)
        weight = f[k]
        isfinite(weight) && weight > 0 || continue
        total_weight += weight

        row = round(Int, row0 - y[k] / meterperpxl_row) + index_offset
        col = round(Int, col0 + x[k] / meterperpxl_col) + index_offset

        if !(1 <= row <= nrows && 1 <= col <= ncols)
            outside_weight += weight
            continue
        end

        class_value = classes_array[row, col]
        if ismissing(class_value)
            missing_class_weight += weight
            continue
        end

        classified_weight += weight
        if class_value
            ice_weight += weight
        else
            feature_weight += weight
        end
    end

    if classified_weight == 0
        return (
            ice_fraction=missing,
            feature_fraction=missing,
            classified_weight_fraction=total_weight > 0 ? classified_weight / total_weight : 0.0,
            outside_weight_fraction=total_weight > 0 ? outside_weight / total_weight : missing,
            missing_class_weight_fraction=total_weight > 0 ? missing_class_weight / total_weight : missing,
            ice_weight=ice_weight,
            feature_weight=feature_weight,
            classified_weight=classified_weight,
            total_weight=total_weight,
        )
    end

    return (
        ice_fraction=ice_weight / classified_weight,
        feature_fraction=feature_weight / classified_weight,
        classified_weight_fraction=total_weight > 0 ? classified_weight / total_weight : 0.0,
        outside_weight_fraction=total_weight > 0 ? outside_weight / total_weight : missing,
        missing_class_weight_fraction=total_weight > 0 ? missing_class_weight / total_weight : missing,
        ice_weight=ice_weight,
        feature_weight=feature_weight,
        classified_weight=classified_weight,
        total_weight=total_weight,
    )
end

function footprint_weighted_class_fractions_dataframe(footprints::AbstractVector,
        classes_array::AbstractMatrix{Union{Missing, Bool}},
        origin_pxl::AbstractVector{<:Real},
        meterperpxl_row::Real,
        meterperpxl_col::Real;
        origin_is_zero_based::Bool=true)

    fractions = [
        footprint_weighted_class_fractions(
            footprint,
            classes_array,
            origin_pxl,
            meterperpxl_row,
            meterperpxl_col;
            origin_is_zero_based=origin_is_zero_based,
        )
        for footprint in footprints
    ]
    return DataFrame(fractions)
end

footprint_class_fractions_all = Dict{Symbol, DataFrame}()
footprint_class_fraction_names = [
    :footprint_class_fractions1,
    :footprint_class_fractions2,
    :footprint_class_fractions3,
    :footprint_class_fractions4,
]

for ix in eachindex(fluxes)
    flux = fluxes[ix]
    fractions = footprint_weighted_class_fractions_dataframe(
        ffp_blocks_all[flux],
        classes_array,
        classes_fluxloc_pxl[ix, :],
        classes_meterperpxl_row,
        classes_meterperpxl_col;
        origin_is_zero_based=classes_coordinates_zero_based,
    )

    insertcols!(
        fractions,
        1,
        :time => block_fluxes_all[flux].time,
        :block_start => block_fluxes_all[flux].block_start,
        :block_end => block_fluxes_all[flux].block_end,
    )

    footprint_class_fractions_all[flux] = fractions

    block_fluxes_all[flux][!, :footprint_ice_fraction] = fractions.ice_fraction
    block_fluxes_all[flux][!, :footprint_feature_fraction] = fractions.feature_fraction
    block_fluxes_all[flux][!, :footprint_classified_weight_fraction] = fractions.classified_weight_fraction
    block_fluxes_all[flux][!, :footprint_outside_weight_fraction] = fractions.outside_weight_fraction
    block_fluxes_all[flux][!, :footprint_missing_class_weight_fraction] = fractions.missing_class_weight_fraction

    @eval $(footprint_class_fraction_names[ix]) = $fractions
end

block_fluxes_output_file = ffp_block.save_block_fluxes_netcdf(block_fluxes_all, station_id; output_dir=output_folder)
println("Saved block fluxes to ", block_fluxes_output_file)
##
=#
#######################################################
###                BLOCK ANALYSIS                   ###
#######################################################

save_csv = true
#######################################################
#read station specifics
station_id = String(stationcfg.require_key(station_config, "id"))
surface_type = String.(stationcfg.require_key(station_config, "surface_type"))
surface_feature = String(stationcfg.require_key(station_config, "block_footprints", "feature"))
feature_file_label = block_analyze.feature_file_label(surface_feature)
#######################################################
#load the block file
block_data = ffp_block.read_block_fluxes_netcdf(station_name; input_dir=output_folder)

#order of the subplots
#left top, right top, left bottom, right bottom
subplot_order = block_analyze.read_subplot_order(station_config, "subplot_order",
    [:fx4, :fx2, :fx3, :fx1])
latent_subplot_order = block_analyze.read_subplot_order(station_config, "latent_subplot_order",
    [:fx3, :fx1])
#######################################################
###                    PLOTS                        ###
#######################################################
wT_limits = Tuple(Float64.(stationcfg.optional_key(station_config, [-20.0, 20.0], "plot", "wT_limits")))
wq_limits = Tuple(Float64.(stationcfg.optional_key(station_config, [-5.0, 5.0], "plot", "wq_limits")))
##
#######################################################
#calculate correlations and export to .csv
correlation_specs = vcat(
    [(flux_name, :wT, "sensible_heat_flux", ρ_air * c_p) for flux_name in subplot_order],
    [(flux_name, :wq, "latent_heat_flux", L_v * 1e-3) for flux_name in latent_subplot_order],
)
flux_feature_correlations = block_analyze.feature_flux_correlations(
    correlation_specs, block_data, station_name, surface_type, heights;
    feature_label=surface_feature)
if save_csv
    CSV.write(joinpath(output_folder, "$(station_file_stem)_block_flux_$(feature_file_label)_correlations.csv"),
        flux_feature_correlations)
end

#######################################################
###     feature fraction and flux DIFFERENCES       ###
#######################################################
#plot time series
wT_limits_diff_timeseries = Tuple(Float64.(stationcfg.optional_key(station_config, [-10.0, 10.0], "block_footprints", "wT_limits_diff_timeseries")))
wq_limits_diff_timeseries = Tuple(Float64.(stationcfg.optional_key(station_config, [-10.0, 10.0], "block_footprints", "wq_limits_diff_timeseries")))

difference_time_series_specs = block_analyze.standard_difference_specs(
    subplot_order, latent_subplot_order, ρ_air * c_p, L_v * 1e-3)

flux_feature_difference_correlations = block_analyze.feature_flux_difference_correlations(
    difference_time_series_specs, block_data, station_name, surface_type, heights;
    feature_label=surface_feature)
if save_csv
    CSV.write(joinpath(output_folder, "$(station_file_stem)_block_flux_$(feature_file_label)_difference_correlations.csv"),
    flux_feature_difference_correlations)
end

#######################################################
#scatter plot
fig_diff_scatter, axs_diff_scatter = PyPlot.subplots(1, 3, figsize=(15, 4.2), sharex=false)
#fig_diff_scatter.suptitle("$(station_label) - Heat flux difference vs $(surface_feature) fraction difference")
for (ax, (flux_names, flux_column, conversion_factor, flux_kind, flux_ylabel)) in
        zip(vec(axs_diff_scatter), difference_time_series_specs)
    flux_ylim = flux_column == :wT ? wT_limits_diff_timeseries : wq_limits_diff_timeseries
    color = flux_column == :wT ? "black" : "C0"
    block_analyze.plot_feature_flux_difference_panel!(ax, block_data, flux_names,
        flux_column, conversion_factor, flux_kind, flux_ylabel, flux_ylim,
        surface_type, heights; color=color, feature_label=surface_feature)
end
axs_diff_scatter[1].set_xlim(0.2, 0.9)
axs_diff_scatter[2].set_xlim(0.55, 1.0)
axs_diff_scatter[3].set_xlim(0.35, 1.0)

subplot_labels = ["a)", "b)", "c)"]
for (ax, label) in zip(vec(axs_diff_scatter), subplot_labels)
    ax.text(-0.2, 1.1, label, transform=ax.transAxes, ha="left", va="top",
    fontsize=14, fontweight="normal", clip_on=false)
end

PyPlot.tight_layout(rect=[0, 0, 1, 0.91])

#fitting linear model to differences
linear_fit_params = []
predictions = []

for (ax, (flux_names, flux_column, conversion_factor, flux_kind, flux_ylabel)) in
        zip(vec(axs_diff_scatter), difference_time_series_specs)

    x, y = block_analyze.feature_flux_difference_xy(
        block_data, flux_names, flux_column, conversion_factor)
    fit = block_analyze.fit_feature_flux_difference(x, y)

    color = flux_column == :wT ? "black" : "C0"
    block_analyze.plot_feature_flux_difference_fit!(ax, fit; color=color)
    ax.legend(loc="best")

    difference_label = block_analyze.sensor_difference_label(flux_names, heights, surface_type)
    push!(linear_fit_params, block_analyze.feature_flux_difference_fit_summary(
        fit, surface_feature, flux_column, flux_kind, difference_label))
    push!(predictions, fit)
end

if save_csv
    CSV.write(joinpath(output_folder,
        "$(station_file_stem)_block_flux_$(feature_file_label)_difference_linear_fit.csv"),
        DataFrame(linear_fit_params))
end

#fig_diff_scatter.savefig(joinpath(output_folder, "$(station_file_stem)_block_flux_$(feature_file_label)_difference_scatter.pdf"), bbox_inches="tight")

######################################################
###                  IR PLOT                       ###
######################################################

"""
    struct ProfileSpec
        row_range :: UnitRange{Int}   # rows that form the vertical axis
        col_range :: UnitRange{Int}   # columns that will be collapsed (median)
        label     :: String           # legend entry (e.g. “Site A”)
    end
"""
struct ProfileSpec
    row_range :: UnitRange{Int}
    col_range :: UnitRange{Int}
    label     :: String
end

"""Convert a folder name like `250412_083015` → DateTime (UTC)."""
function foldername_to_dt(name::String)::DateTime
    parsed_date = DateTime(1900,01,01,00,00,00)
    try
        parsed_date = DateTime(name, dateformat"yymmdd_HHMMSS") + Year(2000)
    catch e
         @debug "Folder $(name) does not follow the expected naming scheme."
    end
    return parsed_date
end

"""
    file_index(fname::String) → Int

`fname` is something like `irdata_0012.nc`.  Returns the integer index
(`12`).  Throws an error if the pattern does not match.
"""
function fileindex(fname::String)
    m = match(r"irdata_(\d{2,})\.nc$", fname)
    m === nothing && error("Unexpected NetCDF filename: $fname")
    return parse(Int, m.captures[1])
end

"""
    find_files(root_dir::String,
               t_start::DateTime,
               t_end::DateTime) → Vector{String}

Return a **sorted** vector of absolute paths to NetCDF files that contain at
least one frame whose timestamp lies inside `[t_start, t_end]`.

Assumptions that the algorithm exploits
----------------------------------------
* Each folder is named `yyyymmdd_HHMMSS` – the start time of the *first*
  NetCDF file (`irdata_0000.nc`) inside that folder.
* Inside a folder the files are named `irdata_XXXX.nc` where `XXXX` is a
  zero‑padded counter that increases monotonically with time (≈ 16–17 s
  between successive files).
* The NetCDF variable `timestamp` stores “days since 1900‑01‑01”.
"""
function find_files(root_dir::String, t_start::DateTime, t_end::DateTime)::Vector{String}
    @info "Scanning $(root_dir) …"

    # ------------------------------------------------------------------
    # 1️⃣  Gather candidate folders (those that start early enough)
    # ------------------------------------------------------------------
    candidate_folders = String[]
    for (folder_root, dirs, _) in walkdir(root_dir; follow_symlinks=false)
        for d in dirs
            dt = foldername_to_dt(d)
            if t_start - Day(2) ≤ dt ≤ t_end
                push!(candidate_folders, joinpath(folder_root, d))
            end
        end
    end

    if isempty(candidate_folders)
        @warn "No candidate folders found under $(root_dir)."
        return String[]
    end

    # Sort folders chronologically (helps deterministic processing)
    sort!(candidate_folders)

    # ------------------------------------------------------------------
    # 2️⃣  Walk through the candidate folders in chronological order
    # ------------------------------------------------------------------
    matched_paths = String[]

    for folder in candidate_folders
        folder_start = foldername_to_dt(basename(folder))
        @debug "Processing folder $(folder) (starts $(folder_start))"

        # ----------------------------------------------------------------
        # 2a️⃣  List the NetCDF files inside the folder (already sorted)
        # ----------------------------------------------------------------
        nc_files = filter(f -> endswith(f, ".nc"),
                          sort(readdir(folder; join=true)))

        # Early‑out flag: once we have seen a file whose *first* timestamp >
        # t_end we can break out of the inner loop (no later file can be useful).
        stop_folder = false

        for nc_path in nc_files
            # ----------------------------------------------------------------
            # 2b️⃣  Quick numeric check using the file name index
            # ----------------------------------------------------------------
            file_index = fileindex(nc_path)

            # Approximate timestamp of the *first* frame in this file:
            #   folder_start + file_index * Δt_per_file
            # Δt_per_file ≈ 500 frames / 30 Hz ≈ 16.666… s
            Δt_per_file = Second(round(Int, 500 / 30))   # 16 s (rounded)
            approx_file_start = folder_start + file_index * Δt_per_file

            if approx_file_start > t_end
                @debug "Skipping remaining files in $(folder) – start $(approx_file_start) > t_end"
                stop_folder = true
                break
            end

            # ----------------------------------------------------------------
            # 2c️⃣  Read ONLY the first timestamp from the NetCDF file
            # ----------------------------------------------------------------
            ds = NCDatasets.Dataset(nc_path, "r")
            try
                ts = ds["timestamp"][:][1] 
                if t_start ≤ ts ≤ t_end
                    push!(matched_paths, nc_path)
                elseif ts > t_end
                    # This file already starts after the interval → stop scanning
                    stop_folder = true
                    break
                end
                # If ts_dt < t_start we simply continue – a later file may match.
            finally
                close(ds)
            end
        end   # ← file loop

        stop_folder && break   # no later folder can contain relevant data
    end   # ← folder loop

    sort!(matched_paths)
    @info "Found $(length(matched_paths)) matching NetCDF files."
    return matched_paths
end

# ----------------------------------------------------------------------
# Load a single frame (as Float64) and apply the conversion factor
# ----------------------------------------------------------------------
"""
    read_frame(path::String, frame_idx::Int) -> Matrix{Float64}

Loads `irdata[:,:,frame_idx]` from `path` and divides by 100 to obtain °C.
Only the requested slice is read, keeping memory usage minimal.
"""
function read_frame(path::String, frame_idx::Int)::Matrix{Float64}
    ds = NCDatasets.Dataset(path, "r")
    try
        # NetCDF stores data in C order (row‑major); NCDatasets returns an
        # Array with the same ordering, so slicing works directly.
        raw = ds["irdata"][:, :, frame_idx]   # 2‑D slice
        return raw ./ 100.0                    # convert to °C
    finally
        close(ds)
    end
end

# ----------------------------------------------------------------------
# Plot the first matching frame for interactive selection
# ----------------------------------------------------------------------
function plot_first_frame(filelist::Vector{String}; tmin=-99.9, tmax=99.9)
    @info "Loading first frame of the earliest file …"
    first_file = filelist[1]
    frame = read_frame(first_file, 1)   # frame index 1 (NetCDF is 1‑based)
    (fig, ax) = PyPlot.subplots()
    if tmin != -99.9 || tmax != 99.9
        im = ax.imshow(frame, origin="upper", cmap=cramericm.batlow, vmin=tmin, vmax=tmax)#; cmap=cmap, origin="upper")
    else
        im = ax.imshow(frame, origin="upper", cmap=cramericm.batlow)#; cmap=cmap, origin="upper")
    end
    ax.set_title("First frame (file: $(basename(first_file)))")
    fig.colorbar(im, ax=ax, label="Temperature (°C)")
    return fig, ax, frame
end

"""
    count_matching_frames_opt(filelist::Vector{String},
                             t_start::DateTime,
                             t_end::DateTime) → Int

Counts how many frames lie inside `[t_start, t_end]` **without opening
most files**.  It assumes that every file except the *last file in its
folder* (and the *last file of the whole selection*) contains exactly
`FRAMES_PER_FILE` frames.  Those two exceptional files are opened to read
their actual `timestamp` length.

The function returns the total number of frames that intersect the interval.
"""
function count_matching_frames(filelist::Vector{String},
                                  t_start::DateTime,
                                  t_end::DateTime,
                                  secs_per_file::Real,
                                  frames_per_file::Int)::Int
    # ------------------------------------------------------------------
    # 1️⃣  Group files by folder (they are already sorted alphabetically)
    # ------------------------------------------------------------------
    sort!(filelist)                                 # ensure deterministic order
    folder_groups = Dict{String,Vector{String}}()

    for p in filelist
        folder = dirname(p)                         # parent directory
        push!(get!(folder_groups, folder, String[]), p)
    end

    total_frames = 0

    # ------------------------------------------------------------------
    # 2️⃣  Process each folder
    # ------------------------------------------------------------------
    for (folder, files) in sort(collect(folder_groups); by=first)
        folder_dt = foldername_to_dt(folder)
        folder_dt === nothing && continue           # skip malformed folders

        # Files inside a folder are already sorted because we sorted `filelist`
        # The *last* file in this folder may be truncated.
        last_file_in_folder = files[end]

        for (i, fpath) in enumerate(files)
            idx = fileindex(basename(fpath))

            # --------------------------------------------------------------
            # Compute the *theoretical* start‑time of this file
            # --------------------------------------------------------------
            file_start = folder_dt + Millisecond(round(Int,
                                 idx * secs_per_file * 1000))

            # By default we assume a full‑size file
            n_frames_this_file = frames_per_file

            # --------------------------------------------------------------
            # Detect the two “exceptional” files that need real inspection
            # --------------------------------------------------------------
            is_last_in_folder = (fpath == last_file_in_folder)
            is_last_overall   = (fpath == filelist[end])

            if (is_last_in_folder) .&& .!(is_last_overall)
                # Open the file and read the actual length of the timestamp vector
                ds = NCDatasets.Dataset(fpath, "r")
                try
                    n_frames_this_file = length(ds["timestamp"][:])
                finally
                    close(ds)
                end
            end
            if is_last_overall
                ds = NCDatasets.Dataset(fpath, "r")
                try
                    timestamps = ds["timestamp"][:]

                    # Identify indices within the interval
                    idx_in = findall(t -> t ≤ t_end, timestamps)
                    n_frames_this_file = length(idx_in)
                finally
                    close(ds)
                end
            end

            frames_here = max(0, n_frames_this_file)
            total_frames += frames_here
        end
    end

    @info "Counted $total_frames frames that intersect the requested interval."
    return total_frames
end

"""
    read_and_collapse!(
        profile_data::Matrix{Float64}
        filelist::Vector{String},
        t_start::DateTime,
        t_end::DateTime,
        row_range::UnitRange{Int},
        col_range::UnitRange{Int},
    ) → Nothing

Fills `profile_data` column‑wise.  Each column corresponds to one frame.
If `length(col_range) > 1` the function first computes the **row‑wise median**
across the selected columns; otherwise the single column is copied unchanged.
Only the rows in `row_range` are ever read from disk.
"""
function read_and_collapse!(profile_data::AbstractArray{Float64},
                           filelist::Vector{String},
                           t_start::DateTime,
                           t_end::DateTime,
                           row_range::UnitRange{Int},
                           col_range::UnitRange{Int})::Nothing
    col_len = length(col_range)
    cur_col = 1                               # which column of profile_data we are filling

    total_frames = size(profile_data, 2)       # number of columns we allocated
    prog = Progress(total_frames; desc="Reading frames", dt=0.5)

    for path in filelist
        ds = NCDatasets.Dataset(path, "r")
        try
            ts_vec = ds["timestamp"][:]

            # Find the indices of frames that belong to the interval
            idxs = findall(t -> t_start .≤ t .≤ t_end, ts_vec)

            # If there are no matching frames in this file, skip it
            isempty(idxs) && continue

            # ----------------------------------------------------------------
            # Load **only** the rows/columns we need, for *all* matching frames
            # ----------------------------------------------------------------
            # NetCDF slicing: var[:, :, idx]  →  (row, col, frame)
            # We request a view of shape (nRows, nCols, nFrames)
            # NCDatasets allows us to pass a tuple of ranges.
            # NOTE: Julia uses 1‑based indexing, same as NetCDF.
            # ----------------------------------------------------------------
            ir = ds["irdata"]
            sub = Float32.(ir[row_range, col_range, idxs])   # 3‑D Array{Float32,3} (or Float64)

            # Convert to °C once (divide by 100)
            sub ./= 100.0

            n_frames = size(sub, 3)

            if col_len == 1
                # No column collapsing – just copy the single column per frame
                # `sub` has dimensions (nRows, 1, nFrames)
                for f = 1:n_frames
                    profile_data[:, cur_col] = view(sub, :, 1, f)
                    cur_col += 1
                    next!(prog)
                end
            else
                # Collapse columns by **row‑wise median** (fast, vectorised)
                # `mapslices(median, sub; dims=2)` returns an (nRows, 1, nFrames) array
                collapsed = mapslices(median, sub; dims=2)   # median across cols
                for f = 1:n_frames
                    profile_data[:, cur_col] = view(collapsed, :, 1, f)
                    cur_col += 1
                    next!(prog)
                end
            end
        finally
            close(ds)
        end
    end

    @assert cur_col - 1 == size(profile_data, 2) "Mismatch in counted frames vs. filled columns."
    @info "Finished loading and collapsing $((cur_col-1)) frames into profile_data."
    return nothing
end

"""
    compute_profile_statistics(profile_data::Matrix{Float64}) →
        NamedTuple{(:median, :q25, :q75, :min, :max), Tuple{Vector{Float64},...}}

`profile_data` is (n_rows × n_total_frames).  Each row contains all the
collapsed temperature values for a given height.  The function returns the
per‑row median, the 25‑th and 75‑th quantiles (IQR), and the min / max.
"""
function compute_profile_statistics(profile_data::AbstractArray{Float64})
    # Quantiles are cheap when applied row‑wise via `mapslices`
    median_vec    = mapslices(x -> mean(x), profile_data; dims=2)[:]
    qlow_vec    = mapslices(x -> quantile(x, 0.25), profile_data; dims=2)[:]
    qhigh_vec    = mapslices(x -> quantile(x, 0.75), profile_data; dims=2)[:]
    min_vec    = mapslices(minimum, profile_data; dims=2)[:]
    max_vec    = mapslices(maximum, profile_data; dims=2)[:]

    return (median = median_vec,
            qlow   = qlow_vec,
            qhigh  = qhigh_vec,
            min    = min_vec,
            max    = max_vec)
end

"""
    prepare_multiple_profiles(root_dir::String,
                              t_start::DateTime,
                              t_end::DateTime,
                              specs::Vector{ProfileSpec}) → Vector{NamedTuple}

For each `ProfileSpec` in `specs` the function:
   1. finds the matching NetCDF files,
   2. counts how many frames belong to the interval,
   3. allocates a reduced matrix,
   4. fills it with the collapsed data,
   5. computes median/IQR/min/max per row.

The returned vector has the same order as `specs`; each element is a
`NamedTuple` with fields `median, q25, q75, min, max, row_range, label`.
"""
function prepare_multiple_profiles(files::Vector{String},
                                   t_start::DateTime,
                                   t_end::DateTime,
                                   specs::Vector{ProfileSpec},
                                   frames_per_file::Int,
                                   sample_rate::Float64)

    secs_per_frame = 1.0 / sample_rate
    secs_per_file   = frames_per_file * secs_per_frame   # ≈ 16.666 s

    # Count total frames (shared by all profiles)
    n_frames = count_matching_frames(files, t_start, t_end, secs_per_file, frames_per_file)

    # Allocate a *single* buffer that will be reused for every profile.
    n_rows_max = maximum(length.(getfield.(specs, :row_range)))   # biggest row span
    profile_buffer = Matrix{Float64}(undef, n_rows_max, n_frames)

    results = NamedTuple[]   # will hold one NamedTuple per spec

    @info "Loading $(size(specs,1)) profiles"

    for spec in specs
        # -----------------------------------------------------------------
        # Resize the buffer to the exact number of rows needed for this spec
        # (the buffer is larger than necessary for some specs – that’s fine).
        # -----------------------------------------------------------------
        rows = spec.row_range
        n_rows = length(rows)
        view_buf = @view profile_buffer[1:n_rows, :]   # cheap view, no copy

        # -----------------------------------------------------------------
        # Fill the view with the collapsed data for this spec
        # -----------------------------------------------------------------
        read_and_collapse!(view_buf, files, t_start, t_end,
                           spec.row_range, spec.col_range)

        # -----------------------------------------------------------------
        # Compute statistics on the *filled* view
        # -----------------------------------------------------------------
        stats = compute_profile_statistics(view_buf)

        # -----------------------------------------------------------------
        # Pack everything we need for plotting
        # -----------------------------------------------------------------
        push!(results, (median = stats.median,
                        qlow    = stats.qlow,
                        qhigh    = stats.qhigh,
                        min    = stats.min,
                        max    = stats.max,
                        row_range = rows,
                        label  = spec.label))
    end

    return results
end

function plot_multi_profiles(first_frame::Matrix{Float64},
                             results::Vector{NamedTuple},
                             specs::Vector{ProfileSpec},
                             t_start::DateTime,
                             t_end::DateTime,
                             custom_title::String = "";
                             tmin::Real = -2.0,
                             tmax::Real = 3.0,
                             pxl2meter_row::Float64 = 0.005,
                             pxl2meter_col::Float64 = 0.005,
                             figsize = (10,5),
                             width_ratios = (1.3, 1.0),
                             colorbar_fraction::Real = 0.06,
                             colorbar_shrink::Real = 0.9,
                             colorbar_pad::Real = 0.2,
                             colorbar_aspect::Real = 30)

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12
    # --------------------------------------------------------------
    # Create the figure with two sub‑plots (profile | reference frame)
    # --------------------------------------------------------------
    fig = PyPlot.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 2, width_ratios=width_ratios, wspace=0.30,)
    ax_ref = fig.add_subplot(gs[1, 1])
    ax_prof = fig.add_subplot(gs[1, 2])
    # --------------------------------------------------------------
    # Colour handling – use Matplotlib's default cycle unless you supply yours
    # --------------------------------------------------------------
    colors = PyPlot.get_cmap("tab10")

    # --------------------------------------------------------------
    # 1️⃣  Plot each profile (left panel)
    # --------------------------------------------------------------
    for (i, res) in enumerate(results)
        col = colors(i-1)

        y = (collect(res.row_range).-res.row_range[1])[end:-1:1] .* pxl2meter_row

        # Median line
        ax_prof.plot(res.median, y; color=col, lw=2, label="$(res.label) – mean")

        # IQR shaded region
        ax_prof.fill_betweenx(y, res.qlow, res.qhigh; color=col, alpha=0.25, label="$(res.label) – IQR")

        # Min / Max dashed lines
        #ax_prof.plot(res.min, y; color=col, ls="--", lw=1, label="$(res.label) – min/max")
        #ax_prof.plot(res.max, y; color=col, ls="--", lw=1)   # same label not needed
    end

    #ax_prof.invert_yaxis()
    ax_prof.set_xlabel(L"T~\mathrm{[^\circ C]}")
    ax_prof.set_ylabel("Profile Height [m]")
    #ax_prof.set_title("Vertical temperature profiles")
    #ax_prof.legend(framealpha=0.8)
    ax_prof.grid(true)
    ax_prof.set_ylim(0,3.5)
    ax_prof.set_xlim(-1.9,-0.9)

    # --------------------------------------------------------------
    # 2️⃣  Plot the reference frame (right panel) and draw coloured rectangles
    # --------------------------------------------------------------
    nrows, ncols = size(first_frame)

    extent = (
        0,                     # left  (metres)
        (ncols-1) * pxl2meter_col,            # right (metres)
        0,                      # bottom
        (nrows-1) * pxl2meter_row            #  top
    )

    im = ax_ref.imshow(first_frame;
                       cmap=cramericm.batlow,
             # origin="upper",
                       vmin=tmin, vmax=tmax,
                       extent=extent)   

    #ax_ref.set_title("Sample frame with profile locations")
    ax_ref.set_xlabel(L"x~\mathrm{[m]}")
    ax_ref.set_ylabel("Frame height [m]")

    # One rectangle per profile – colour matches the curve
    for (i, spec) in enumerate(specs)
        col = colors(i-1)

        # Rectangle geometry (Matplotlib expects lower‑left corner)
        rect = patches.Rectangle(
            (first(spec.col_range) * pxl2meter_col,
            (nrows - last(spec.row_range)) * pxl2meter_row),
            length(spec.col_range) * pxl2meter_col,
            length(spec.row_range) * pxl2meter_row,
            linewidth = 2,
            edgecolor = col,
            facecolor = "none",
            zorder = 10,
            alpha = 0.4)
        ax_ref.add_patch(rect)
    end

    # Colour‑bar for the thermal image
    cbar = fig.colorbar(im, ax=ax_ref, orientation="horizontal", fraction=colorbar_fraction,
        shrink=colorbar_shrink, pad=colorbar_pad, aspect=colorbar_aspect,)
    cbar.set_label(L"T~\mathrm{[^\circ C]}")

    # --------------------------------------------------------------
    #  Add wind direction arrow
    # --------------------------------------------------------------
    # Calculate arrow position (bottom-left area)
    arrow_x_start = extent[1] + 0.05 * (extent[2] - extent[1])  
    arrow_y = extent[3] + 0.88 * (extent[4] - extent[3])         
    arrow_length = 0.15 * (extent[2] - extent[1])                # 15% of width
    
    ax_ref.annotate("", 
                   xy=(arrow_x_start + arrow_length, arrow_y),
                   xytext=(arrow_x_start, arrow_y),
                   arrowprops=Dict("arrowstyle"=>"->", 
                                  "lw"=>2, 
                                  "color"=>"white",
                                  "mutation_scale"=>20))
    
    # Add "Wind" label
    ax_ref.text(arrow_x_start + arrow_length/2, 
               arrow_y + 0.03 * (extent[4] - extent[3]),
               "wind",
               horizontalalignment="center",
               verticalalignment="bottom",
               color="white",
               #fontsize=10,
               fontweight="bold")

    # --------------------------------------------------------------
    # 3️⃣  Overall title (interval + optional custom part)
    # --------------------------------------------------------------
    interval_str = @sprintf("%s – %s",
                            Dates.format(t_start, "yyyy‑mm‑dd HH:MM:SS"),
                            Dates.format(t_end,   "yyyy‑mm‑dd HH:MM:SS"))
    suptitle = "IR Temperature Profiles  [$interval_str]"
    if !isempty(custom_title)
        suptitle *= "  –  $custom_title"
    end
    #fig.suptitle(suptitle, fontsize=12)

    for (ax, label) in zip((ax_ref, ax_prof), ("a)", "b)"))
        ax.text(-0.2, 1.03, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=14, clip_on=false, zorder=20)
    end

    PyPlot.tight_layout()#rect=[0, 0.03, 1, 0.95])
    return fig
end

# ------------------------------------------------------------------
# 1️⃣  Input variables
# ------------------------------------------------------------------
root_dir = "/media/haugened/8f2c07c5-c77d-4280-99f1-e932f1d72d33/contrasts/converted/" #dir containing the converted folders
t_start = DateTime(2025,07,15,10,15,00) #start time for analysis
t_end   = t_start + Minute(30)          #end time for analysis

pxl2meter_row = 0.0069124424   # m per pixel vertically
pxl2meter_col = pxl2meter_row #0.006   # m per pixel horizontally
#enter col_range and row_range below!!!
#=
#pxl2meter values
1a (1000-1030) 0.0064655172
2a (1000-1045) 0.0069124424
3a (0900-0930) row: 0.0083, col: 0.0080
3a 250721_111824 row: 0.00638, col 0.00647
=#
# ----------------------------------------------------------------------
# Configuration (adjust once if you like)
# ----------------------------------------------------------------------
const LOG_LEVEL = Info          # change to Debug for more output
global_logger(ConsoleLogger(stderr, LOG_LEVEL))

const FRAMES_PER_FILE = 500                # normal size
const SAMPLE_RATE     = 30.0                # Hz

# ----------------------------------------------------------------------
# Main 
# ----------------------------------------------------------------------
println("\n=== IR Vertical‑Profile Plotter ===\n")

# ------------------------------------------------------------------
# 2️⃣  Locate matching files
# ------------------------------------------------------------------
files = find_files(root_dir, t_start, t_end)
isempty(files) && error("No files found for the given interval.")

# ------------------------------------------------------------------
# 3️⃣  Show first frame for visual selection
# ------------------------------------------------------------------
#fig, ax, first_frame = plot_first_frame(files; tmin=-1.8, tmax=0.8)
#PyPlot.show(fig)

# ------------------------------------------------------------------
# 4️⃣  Ask user for column & row ranges (inclusive)
# ------------------------------------------------------------------
#row, col
#=
#1a
specs = [
    irev.ProfileSpec(400:690, 250:260, "A"),
    irev.ProfileSpec(110:715, 800:810, "B")   # col_range length == 1 → no median collapse
 
]
=#
#2a
specs = [
    ProfileSpec(460:730, 290:300, "A"),
#    irev.ProfileSpec(460:740, 450:460, "B"),
#    irev.ProfileSpec(200:740, 650:660, "C"),
    ProfileSpec(200:730, 770:780, "B")   # col_range length == 1 → no median collapse
]

#=
#3a (250721_085440)
specs = [
    irev.ProfileSpec(375:605, 295:305, "A"),
    irev.ProfileSpec(130:610, 645:655, "B")   # col_range length == 1 → no median collapse
]
=#
#=
#3a (250721_111824)
specs = [
    irev.ProfileSpec(380:670, 840:860, "A"),
    irev.ProfileSpec(390:685, 640:660, "B"),   # col_range length == 1 → no median collapse
    irev.ProfileSpec(100:685, 240:260, "C")   # col_range length == 1 → no median collapse
]=#
# ------------------------------------------------------------------
# 6️⃣  Compute profile statistics
# ------------------------------------------------------------------
results = prepare_multiple_profiles(files, t_start, t_end, specs, FRAMES_PER_FILE, SAMPLE_RATE)

# ------------------------------------------------------------------
# 7  Plot
# ------------------------------------------------------------------
fig = plot_multi_profiles(first_frame, results, specs, t_start, t_end, "2a"; tmin=-1.8, tmax=0.8, pxl2meter_row, pxl2meter_col);
fig.savefig(joinpath(output_folder, "2a_t_profiles_1015_1045.pdf"), bbox_inches="tight")


###################################################################################################
#Figures for ridge section (Figs. xxx)
#Heat fluxes station 3a2
#=
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
using Images
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
LogNorm = pyimport("matplotlib.colors")
mpimg = pyimport("matplotlib.image")
mpl_scatter_density = pyimport("mpl_scatter_density")

station_name = "3a" # Change this to "1a", "1b_1", ..., "3d" before including the script.

importdir = joinpath(@__DIR__, "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "station_config.jl"))
include(joinpath(importdir, "src", "kljun_ffp.jl"))
include(joinpath(importdir, "src", "footprint_plotting.jl"))
import .turb
import .gen
import .stationcfg
import .kljun
@pyinclude(joinpath(importdir, "src", "kljun_ffp_climatology.py"))

timestep = Millisecond(50)
ρ_air = 1.2 #kg m^{-3}
c_p = 1004 #J kg^{-1} K^{-1}
L_v = 2500e3 #J kg^{-1} (approx @0°C)

######################################################
#variables

station_config = stationcfg.load_station_config(station_name)
station_label = stationcfg.station_label(station_config)
station_file_stem = stationcfg.station_file_stem(station_config)

datapath = String(stationcfg.require_key(station_config, "data_root"))
#datapath = "/home/michi/Documents/slf/CONTRASTS25/data/processed/preproc/"

#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
nanfile_stam = String(stationcfg.optional_key(station_config, "/home/haugened/Documents/data/CONTRASTS/nan_periods", "manual_nanmask", "root"))
nanfiles = joinpath.(nanfile_stam, String.(stationcfg.optional_key(
    station_config,
    ["t1irg_nan.csv", "t1csat_nan.csv", "t2irg_nan.csv", "t2csat_nan.csv"],
    "manual_nanmask",
    "files",
)))
plot_output_root = stationcfg.plot_root(station_config)

######################################################

println()
println("-----------S-T-A-R-T-------------")
println("Loading station ", station_label, " from ", station_config["config_file"])

######################################################
###            LOADING & PREPROCESSING             ###
######################################################
#select data and measurement period to be evaluated
evalstart = DateTime(2025,07,21,09,00,00)
evalend   = stationcfg.station_datetime(station_config, "evalend")
#evalend = evalstart + Day(10)

data_files = String.(stationcfg.require_key(station_config, "data_files"))

evaldfs = [
    turb.readturbasnetcdf(joinpath(datapath, data_file), evalstart, evalend)
    for data_file in data_files
]
(evaldf1, evaldf2, evaldf3, evaldf4) = evaldfs

# Instrument metadata for labeling (corresponds to evaldf1-4 / fx1-4)
srf_type = String.(stationcfg.require_key(station_config, "surface_type"))
heights = Float64.(stationcfg.require_key(station_config, "heights"))
instr_type = String.(stationcfg.optional_key(station_config, ["IRG", "CSAT", "IRG", "CSAT"], "instrument_type"))
instr_labels = ["ridge 0.9m", "ridge 2.0m", "upwind 1.1m", "upwind 2.0m"]

#show statistics about missing or NaN-data
for evaldf in evaldfs
    turb.printmissstats(evaldf)
end

#interpolate missing
evaldfs = [turb.interpolatemissing(evaldf) for evaldf in evaldfs]
(evaldf1, evaldf2, evaldf3, evaldf4) = evaldfs

wds = [turb.winddir(evaldf) for evaldf in evaldfs]

#double rotation
for evaldf in evaldfs
    turb.drdf!(evaldf, periodwise=false)
end

#CONTRASTS: apply manual nan-mask
if stationcfg.optional_key(station_config, true, "manual_nanmask", "enabled")
    nanmask_columns = String.(stationcfg.optional_key(station_config, ["u", "v", "w", "T"], "manual_nanmask", "columns"))
    manual_periods = [turb.read_nantimes_csv(nanfile) for nanfile in nanfiles]
    for (evaldf, manual_period) in zip(evaldfs, manual_periods)
        isbad_quality = turb.manual_nanmask(manual_period, evaldf.time)
        evaldf[isbad_quality, nanmask_columns] .= NaN
    end
end

#tower rotation correction for station 3d
stationcfg.apply_station_wind_direction_rotations!(wds, station_config)
(wd1, wd2, wd3, wd4) = wds

#Fluxes (from turb_fluxes.jl)

######################################################
###              TURBULENT FLUXES                  ###
######################################################
turb.missing2nan!(evaldf1)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf3)
turb.missing2nan!(evaldf4)

#Reynolds averaging times (determine usind MRD)
ra1 = Second(400) #Millisecond(2^11*timestep) #Second(330)
ra2 = Second(400)
ra3 = Second(400) #Millisecond(2^11*timestep) #Second(320)
ra4 = Second(400)

fx1_raw = turb.turbflux(evaldf1, ra1, 1.0)
fx2_raw = turb.turbflux(evaldf2, ra2, 1.0)
fx3_raw = turb.turbflux(evaldf3, ra3, 1.0)
fx4_raw = turb.turbflux(evaldf4, ra4, 1.0)

#averaging
fx1 = turb.avgflux(fx1_raw, ra1, true, 0.1)
fx2 = turb.avgflux(fx2_raw, ra2, true, 0.1)
fx3 = turb.avgflux(fx3_raw, ra3, true, 0.1)
fx4 = turb.avgflux(fx4_raw, ra4, true, 0.1)

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
##

#variables
names = [:evaldf1, :evaldf2, :evaldf3, :evaldf4]
meas_heights = Float64.(stationcfg.optional_key(station_config, heights, "footprint", "measurement_heights"))
pbl_height = 200.0
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

    #keep only periods for which every required footprint input is finite
    valid_blocks = map(eachindex(umean)) do j
        all(isfinite, (umean[j], ol[j], sigmav[j],ustar[j], wind_dir[j],))
    end

    nvalid = sum(valid_blocks)

    if nvalid == 0
        error("No valid periods available for $(names[ix]).")
    end

    println("Using $nvalid of $nrblocks periods; ", nrblocks-nvalid, " period(s) containing NaN/Inf were ignored.")

    output = py"FFP_climatology"(meas_heights[ix], nothing, PyVector(umean[valid_blocks]), PyVector(h[valid_blocks]), PyVector(ol[valid_blocks]),
    PyVector(sigmav[valid_blocks]), PyVector(ustar[valid_blocks]), PyVector(wind_dir[valid_blocks]), PyVector(domain), dx, dy, nx, ny, PyVector(rs), rslayer, smooth_data, crop, pulse, verbosity, fig)
    if output["flag_err"] != 0
        @warn("Error flag set to true! There is an error!")
    end
    @eval $(outnames[ix]) = output
end

###############################################
#plotting the footprint on the ortho-mosaic

fileorthomosaic = String(stationcfg.require_key(station_config, "footprint", "orthomosaic"))
orthomosaic = mpimg.imread(fileorthomosaic);
orthomosaic_jl = load(fileorthomosaic);

#=
PyPlot.pygui(true)
PyPlot.imshow(orthomosaic)
PyPlot.pygui(false)
=#
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

PyPlot.pygui(true)
ctab10 = PyPlot.cm.tab10
ffp_fig = PyPlot.figure()
ax1 = ffp_fig.add_subplot(111)
#ax1.set_title("Station $(station_label) Flux footprints")
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
ax1.legend(loc="lower right")
ax1.set_xlim(-60,50)
ax1.set_ylim(-90,30)
PyPlot.tight_layout()
#PyPlot.savefig(joinpath("/home/haugened/Documents/data/CONTRASTS/plots/paper_CONTRASTS_26/", "3a2_ffp.pdf"), bbox_inches="tight")

##

######################################################
# Flux time series with dual subplots for latent and sensible heat fluxes (CONTRASTS)
##
fig, (ax1, ax2) = PyPlot.subplots(2, 1, figsize=(10, 8), sharex=true)

# Define colors for consistent legend
colors = ["C0", "C1", "C2", "C3"]  # Default matplotlib color cycle

#define plotting step size
step = 20*60 #every 1min

#y-axis limits
wT_limits = Tuple(Float64.(stationcfg.optional_key(station_config, [-20.0, 20.0], "plot", "wT_limits")))
wq_limits = Tuple(Float64.(stationcfg.optional_key(station_config, [-5.0, 5.0], "plot", "wq_limits")))

# Upper subplot - Buoyancy fluxes (sensible heat)
wt1 = ax1.plot(fx1.time[1:step:end], fx1.wT[1:step:end] .* (ρ_air * c_p), color=colors[1])
wt2 = ax1.plot(fx2.time[1:step:end], fx2.wT[1:step:end] .* (ρ_air * c_p), color=colors[2])
wt3 = ax1.plot(fx3.time[1:step:end], fx3.wT[1:step:end] .* (ρ_air * c_p), color=colors[3])
wt4 = ax1.plot(fx4.time[1:step:end], fx4.wT[1:step:end] .* (ρ_air * c_p), color=colors[4])
ax1.set_ylabel(L"Q_H ~\mathrm{[W~m^{-2}]}")
ax1.grid()
ax1.set_ylim(wT_limits)

# Lower subplot - Latent heat fluxes
wq1 = ax2.plot(fx1.time[1:step:end], fx1.wq[1:step:end] .* (L_v * 1e-3), color=colors[1])
wq2 = ax2.plot(fx3.time[1:step:end], fx3.wq[1:step:end] .* (L_v * 1e-3), color=colors[3])
ax2.set_ylabel(L"Q_E ~\mathrm{[W~m^{-2}]}")
ax2.set_xlabel("Time")
ax2.xaxis_date()
ax2.grid()
ax2.set_ylim(wq_limits)

# Create a single legend for the entire figure
handles = [wt1[1], wt2[1], wt3[1], wt4[1]]  # Get line objects
labels = instr_labels
ax1.legend(handles, labels)#, loc="upper right", bbox_to_anchor=(1.0, 1))
ax2.legend([wq1[1], wq2[1]], [instr_labels[1], instr_labels[3]])

# Optional: Uncomment these lines if you want to set specific time limits
# ax2.set_xlim(DateTime(2021, 05, 31, 10, 30, 00), DateTime(2021, 05, 31, 20, 00))

# Format dates on x-axis if needed
# majorlocator = pydates.HourLocator(interval=1)
# minorlocator = pydates.MinuteLocator(interval=15)
# ax2.xaxis.set_major_locator(majorlocator)
# ax2.xaxis.set_minor_locator(minorlocator)
# date_format = pydates.DateFormatter("%H:%M")
# ax2.xaxis.set_major_formatter(date_format)
# fig.autofmt_xdate()

PyPlot.tight_layout()
output_folder = "/home/haugened/Documents/data/CONTRASTS/plots/paper_CONTRASTS_26/"
#PyPlot.savefig(joinpath(output_folder, "3a2_heat_fluxes.pdf"), bbox_inches="tight")

##############################
#Scatter plot height comparison heat fluxes CONTRASTS
# Comparing sensible heat fluxes at different heights for each tower
##
scatter1 = fx1.wT[.!isnan.(fx1.wT) .&& .!isnan.(fx2.wT)]
scatter2 = fx2.wT[.!isnan.(fx1.wT) .&& .!isnan.(fx2.wT)]
scatter3 = fx3.wT[.!isnan.(fx3.wT) .&& .!isnan.(fx4.wT)]
scatter4 = fx4.wT[.!isnan.(fx3.wT) .&& .!isnan.(fx4.wT)]

scatter1 .*= ρ_air .* c_p
scatter2 .*= ρ_air .* c_p
scatter3 .*= ρ_air .* c_p
scatter4 .*= ρ_air .* c_p

ratio_ht1 = scatter2 ./ scatter1
ratio_ht2 = scatter4 ./ scatter3

fig_ht = PyPlot.figure(figsize=(7.5, 3.5))
gs_ht = gridspec.GridSpec(1, 2, width_ratios=(1, 1))

# --- Row 1: Scatter density ---
# Column 1: Tower 1 (fx1 vs fx2)
ax_ht1 = fig_ht.add_subplot(gs_ht[1, 1], projection="scatter_density")
ax_ht1.set_title("ridge")
ax_ht1.axhline(0, color="grey", alpha=0.6)
ax_ht1.axvline(0, color="grey", alpha=0.6)
ax_ht1.plot([-25,25], [-25,25], color="grey", alpha=0.3)
ax_ht1.scatter_density(scatter1, scatter2, color="black", vmin=0, vmax=1500)
lim_ht1 = (-25, 25)
ax_ht1.set_xlim(lim_ht1)
ax_ht1.set_ylim(lim_ht1)
ax_ht1.set_xlabel("\$Q_{H,$(heights[1])m}~\\mathrm{[W~m^{-2}]}\$")
ax_ht1.set_ylabel("\$Q_{H,$(heights[2])m}~\\mathrm{[W~m^{-2}]}\$")
ax_ht1.grid()
ax_ht1.set_aspect("equal", adjustable="box")

# Column 2: Tower 2 (fx3 vs fx4)
ax_ht2 = fig_ht.add_subplot(gs_ht[1, 2], projection="scatter_density")
ax_ht2.set_title("upwind")
ax_ht2.axhline(0, color="grey", alpha=0.6)
ax_ht2.axvline(0, color="grey", alpha=0.6)
ax_ht2.plot([-25,25], [-25,25], color="grey", alpha=0.3)
ax_ht2.scatter_density(scatter3, scatter4, color="black", vmin=0, vmax=1500)
lim_ht2 = (-6, 6)
ax_ht2.set_xlim(lim_ht2)
ax_ht2.set_ylim(lim_ht2)
ax_ht2.set_xlabel("\$Q_{H,$(heights[3])m}~\\mathrm{[W~m^{-2}]}\$")
ax_ht2.set_ylabel("\$Q_{H,$(heights[4])m}~\\mathrm{[W~m^{-2}]}\$")
ax_ht2.grid()
ax_ht2.set_aspect("equal", adjustable="box")

PyPlot.tight_layout()
#PyPlot.savefig(joinpath("/home/haugened/Documents/data/CONTRASTS/plots/paper_CONTRASTS_26/", "3a2_corr_per_tower.pdf"), bbox_inches="tight")

##
#PyPlot.savefig(joinpath(output_folder_ht, "$(station_file_stem)_per_tower.pdf"), bbox_inches="tight")

#calculating correlations for height comparison
cor_ht1 = nancor(fx1.wT, fx2.wT)
cor_ht2 = nancor(fx3.wT, fx4.wT)

#calculating mean bias for height comparison
mad_ht1 = nanmean(scatter2) - nanmean(scatter1)
mad_ht2 = nanmean(scatter4) - nanmean(scatter3)
##

##############################

##
# Momentum fluxes and friction velocity plot
fig, ax1 = PyPlot.subplots(1, 1, figsize=(8.5,3.5), sharex=true)

# Define colors for consistent legend
colors = ["C0", "C1", "C2", "C3"]

# Define plotting step size
step = 20*60  # every 1min

# Upper subplot - Vertical momentum fluxes
uw1 = ax1.plot(fx1.time[1:step:end], fx1.uw[1:step:end].*100, color=colors[1])
uw2 = ax1.plot(fx2.time[1:step:end], fx2.uw[1:step:end].*100, color=colors[2])
uw3 = ax1.plot(fx3.time[1:step:end], fx3.uw[1:step:end].*100, color=colors[3])
uw4 = ax1.plot(fx4.time[1:step:end], fx4.uw[1:step:end].*100, color=colors[4])
ax1.set_ylabel(L"\overline{u'w'} ~\mathrm{[10^{-2}~m^2~s^{-2}]}")
ax1.set_xlabel("Time")
ax1.xaxis_date()
ax1.grid()

# Create a single legend for the entire figure
handles = [uw1[1], uw2[1], uw3[1], uw4[1]]
#labels same as above
labels = instr_labels
ax1.legend(handles, labels)

PyPlot.tight_layout()
#PyPlot.savefig(joinpath("/home/haugened/Documents/data/CONTRASTS/plots/paper_CONTRASTS_26/", "3a2_momentum_fluxes.pdf"), bbox_inches="tight")
##
#############################
#=
##
#just, if additionally needed
# Momentum fluxes and friction velocity plot with uw and vw
fig, (ax1, ax2, ax3) = PyPlot.subplots(3, 1, figsize=(10, 8), sharex=true)

# Define colors for consistent legend
colors = ["C0", "C1", "C2", "C3"]

# Define plotting step size
step = 20*60  # every 1min

# Upper subplot - Vertical momentum fluxes u'w'
ax1.set_title("Vertical Momentum Fluxes and Friction Velocity")
uw1 = ax1.plot(fx1.time[1:step:end], fx1.uw[1:step:end], color=colors[1])
uw2 = ax1.plot(fx2.time[1:step:end], fx2.uw[1:step:end], color=colors[2])
uw3 = ax1.plot(fx3.time[1:step:end], fx3.uw[1:step:end], color=colors[3])
uw4 = ax1.plot(fx4.time[1:step:end], fx4.uw[1:step:end], color=colors[4])
ax1.set_ylabel(L"\overline{u'w'} ~\mathrm{[m^2~s^{-2}]}")
ax1.grid()

# Middle subplot - Vertical momentum fluxes v'w'
vw1 = ax2.plot(fx1.time[1:step:end], fx1.vw[1:step:end], color=colors[1])
vw2 = ax2.plot(fx2.time[1:step:end], fx2.vw[1:step:end], color=colors[2])
vw3 = ax2.plot(fx3.time[1:step:end], fx3.vw[1:step:end], color=colors[3])
vw4 = ax2.plot(fx4.time[1:step:end], fx4.vw[1:step:end], color=colors[4])
ax2.set_ylabel(L"\overline{v'w'} ~\mathrm{[m^2~s^{-2}]}")
ax2.grid()

# Lower subplot - Friction velocity
us1 = ax3.plot(fx1.time[1:step:end], fx1.u_star[1:step:end], color=colors[1])
us2 = ax3.plot(fx2.time[1:step:end], fx2.u_star[1:step:end], color=colors[2])
us3 = ax3.plot(fx3.time[1:step:end], fx3.u_star[1:step:end], color=colors[3])
us4 = ax3.plot(fx4.time[1:step:end], fx4.u_star[1:step:end], color=colors[4])
ax3.set_ylabel(L"u_*~\mathrm{[m~s^{-1}]}")
ax3.set_xlabel("Time")
ax3.xaxis_date()
ax3.grid()

# Create a single legend for the entire figure
handles = [uw1[1], uw2[1], uw3[1], uw4[1]]
#labels same as above
ax1.legend(handles, labels)

PyPlot.tight_layout()
##
=#
=#