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
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter, Distributed
using Images
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
LogNorm = pyimport("matplotlib.colors")
mpimg = pyimport("matplotlib.image")
mpl_scatter_density = pyimport("mpl_scatter_density")

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
output_folder = stationcfg.plot_dir(station_config, "wT_wq")
#PyPlot.savefig(joinpath("/home/haugened/Documents/data/CONTRASTS/plots/paper_CONTRASTS_26/", "1c_heat_fluxes.pdf"), bbox_inches="tight")


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
block_data = ffp_block.read_block_fluxes_netcdf(station_name)

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
output_folder = "/home/haugened/Documents/data/CONTRASTS/plots/paper_CONTRASTS_26/"
mkpath(output_folder)
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
axs_diff_scatter[2].set_xlim(0.5, 1.0)
axs_diff_scatter[3].set_xlim(0.3, 1.0)
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

if save_figs
    fig_diff_scatter.savefig(joinpath(output_folder,
    "$(station_file_stem)_block_flux_$(feature_file_label)_difference_scatter.pdf"), bbox_inches="tight")
end

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
output_folder = stationcfg.plot_dir(station_config, "wT_wq")
#PyPlot.savefig(joinpath("/home/haugened/Documents/data/CONTRASTS/plots/paper_CONTRASTS_26/", "3a2_heat_fluxes.pdf"), bbox_inches="tight")

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