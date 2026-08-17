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

station_name = "3a" # Change this to "1a", "1b_1", ..., "3d" before including the script.

importdir = joinpath(@__DIR__, "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "station_config.jl"))
import .turb
import .gen
import .stationcfg
#PyPlot.pygui(true)

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

######################################################
###               LOADING SLOW DATA                ###
######################################################
#=
"""
    export_1s_mean_winddirs(output_file, wd1, wd2, wd3, wd4)

Calculate circular one-second means for the four wind-direction DataFrames and
write them to a CSV with columns `timestamp`, `wd1`, `wd2`, `wd3`, and `wd4`.
"""
function export_1s_mean_winddirs(output_file, wd1, wd2, wd3, wd4)
    function one_second_means(wd, column_name)
        per_second = DataFrame(timestamp=floor.(wd.time, Second), direction=wd.α)
        return combine(
            groupby(per_second, :timestamp),
            :direction => (α -> turb.mean_winddir(collect(skipmissing(α)))) => column_name,
        )
    end

    means = map(
        one_second_means,
        (wd1, wd2, wd3, wd4),
        (:wd1, :wd2, :wd3, :wd4),
    )
    output = reduce((left, right) -> outerjoin(left, right, on=:timestamp), means)
    sort!(output, :timestamp)
    CSV.write(output_file, output)
    return output
end

# Example:
export_1s_mean_winddirs("/home/haugened/Documents/data/CONTRASTS/plots/wind_roses/3a_wind_dir_1s.csv", wds[1], wds[2], wds[3], wds[4])
=#

######################################################
###                   PLOTTING                     ###
######################################################

#plot wind roses for CONTRASTS

windrose = pyimport("windrose")
np = pyimport("numpy")


"""
    plot_windrose_panel(wd1, ws1, wd2, ws2, wd3, ws3, wd4, ws4; 
                        bins=nothing, cmap=nothing, figsize=(10, 10))

Create a 2x2 panel of wind roses.

Layout:
    [upper left: tower 1 CSAT]  [upper right: tower 2 CSAT]
    [lower left: tower 1 IRG]   [lower right: tower 2 IRG]

# Arguments
- `wd1, ws1`: Wind direction/speed for tower 1 IRG (lower left)
- `wd2, ws2`: Wind direction/speed for tower 1 CSAT (upper left)
- `wd3, ws3`: Wind direction/speed for tower 2 IRG (lower right)
- `wd4, ws4`: Wind direction/speed for tower 2 CSAT (upper right)
- `bins`: Wind speed bins (default: automatic)
- `cmap`: Colormap (default: "viridis")
- `figsize`: Figure size tuple

# Returns
- `fig`: Matplotlib figure object
"""
function plot_windrose_panel(wd1, ws1, wd2, ws2, wd3, ws3, wd4, ws4;
                             bins=nothing, cmap=nothing, figsize=(10, 10))
    
    # Set defaults
    isnothing(cmap) && (cmap = PyPlot.cm.viridis)
    
    # Data and labels for each panel: (row, col, wd, ws, title)
    panels = [
        (1, 1, wd2, ws2, "tower 1 CSAT"),  # upper left
        (1, 2, wd4, ws4, "tower 2 CSAT"),  # upper right
        (2, 1, wd1, ws1, "tower 1 IRG"),   # lower left
        (2, 2, wd3, ws3, "tower 2 IRG"),   # lower right
    ]
    
    fig = PyPlot.figure(figsize=figsize)
    
    for (row, col, wd, ws, title) in panels
        # Calculate subplot index (1-based, row-major)
        idx = (row - 1) * 2 + col
        
        # Filter out NaN values
        valid_mask = .!isnan.(wd) .& .!isnan.(ws)
        wd_clean = np.array(wd[valid_mask])
        ws_clean = np.array(ws[valid_mask])
        
        # Create WindroseAxes
        ax = fig.add_subplot(2, 2, idx, projection="windrose")
        
        # Plot stacked histogram with percentage normalization
        if isnothing(bins)
            ax.bar(wd_clean, ws_clean, normed=true, opening=0.8, edgecolor="white", cmap=cmap)
        else
            ax.bar(wd_clean, ws_clean, normed=true, opening=0.8, edgecolor="white", 
                   bins=bins, cmap=cmap)
        end
        
        ax.set_title(title, fontweight="bold")
        ax.set_legend(loc="lower right", fontsize=8)
    end

    fig.suptitle("Station $(station_label)", fontsize=14, fontweight="bold")
    PyPlot.tight_layout()
    
    return fig
end

##
@info("No double rotation before wind rose plot!")

#calculate wind speeds
ws1 = sqrt.(evaldf1.u .^2 + evaldf1.v .^2 + evaldf1.w .^2);
ws2 = sqrt.(evaldf2.u .^2 + evaldf2.v .^2 + evaldf2.w .^2);
ws3 = sqrt.(evaldf3.u .^2 + evaldf3.v .^2 + evaldf3.w .^2);
ws4 = sqrt.(evaldf4.u .^2 + evaldf4.v .^2 + evaldf4.w .^2);

#wind directions
wds = [turb.winddir(evaldf) for evaldf in evaldfs]
apply_station_wind_direction_rotations!(wds, station_config)
(wd1, wd2, wd3, wd4) = wds

#block average
block_length = Minute(1)
(ws1_block_time, ws1_block) = gen.block_average(evaldf1.time, ws1, block_length)
(ws2_block_time, ws2_block) = gen.block_average(evaldf2.time, ws2, block_length)
(ws3_block_time, ws3_block) = gen.block_average(evaldf3.time, ws3, block_length)
(ws4_block_time, ws4_block) = gen.block_average(evaldf4.time, ws4, block_length)

(wd1_block_time, wd1_block) = gen.block_average(wd1.time, wd1.α, block_length)
(wd2_block_time, wd2_block) = gen.block_average(wd2.time, wd2.α, block_length)
(wd3_block_time, wd3_block) = gen.block_average(wd3.time, wd3.α, block_length)
(wd4_block_time, wd4_block) = gen.block_average(wd4.time, wd4.α, block_length)

#plot
# Assuming you have your block-averaged data ready
fig = plot_windrose_panel(wd1_block, ws1_block, wd2_block, ws2_block, wd3_block, ws3_block, wd4_block, ws4_block;
                          bins=collect(0:2:10),  # optional custom bins
                          figsize=(12, 10))

output_folder = stationcfg.plot_dir(station_config, "wind_roses")
PyPlot.savefig(joinpath(output_folder, "$(station_file_stem).pdf"), bbox_inches="tight")
##
#########################################################
#Plot wind speed, direction, sonic temperature, and water vapor density
##
block_length = Minute(10)
percentiles = (5, 95) #for max/min shading

# Calculate wind speeds
ws1 = sqrt.(evaldf1.u .^ 2 .+ evaldf1.v .^ 2 .+ evaldf1.w .^ 2)
ws2 = sqrt.(evaldf2.u .^ 2 .+ evaldf2.v .^ 2 .+ evaldf2.w .^ 2)
ws3 = sqrt.(evaldf3.u .^ 2 .+ evaldf3.v .^ 2 .+ evaldf3.w .^ 2)
ws4 = sqrt.(evaldf4.u .^ 2 .+ evaldf4.v .^ 2 .+ evaldf4.w .^ 2)

# Compute block statistics for wind speed
ws1_time, ws1_mean, ws1_lo, ws1_hi = gen.block_stats(evaldf1.time, ws1, block_length; percentiles=percentiles)
ws2_time, ws2_mean, ws2_lo, ws2_hi = gen.block_stats(evaldf2.time, ws2, block_length; percentiles=percentiles)
ws3_time, ws3_mean, ws3_lo, ws3_hi = gen.block_stats(evaldf3.time, ws3, block_length; percentiles=percentiles)
ws4_time, ws4_mean, ws4_lo, ws4_hi = gen.block_stats(evaldf4.time, ws4, block_length; percentiles=percentiles)

# Compute block statistics for sonic temperature
T1_time, T1_mean, T1_lo, T1_hi = gen.block_stats(evaldf1.time, evaldf1.T, block_length; percentiles=percentiles)
T2_time, T2_mean, T2_lo, T2_hi = gen.block_stats(evaldf2.time, evaldf2.T, block_length; percentiles=percentiles)
T3_time, T3_mean, T3_lo, T3_hi = gen.block_stats(evaldf3.time, evaldf3.T, block_length; percentiles=percentiles)
T4_time, T4_mean, T4_lo, T4_hi = gen.block_stats(evaldf4.time, evaldf4.T, block_length; percentiles=percentiles)

# Compute block statistics for water vapor (IRG sensors only)
h2o1_time, h2o1_mean, h2o1_lo, h2o1_hi = gen.block_stats(evaldf1.time, evaldf1.h2o, block_length; percentiles=percentiles)
h2o3_time, h2o3_mean, h2o3_lo, h2o3_hi = gen.block_stats(evaldf3.time, evaldf3.h2o, block_length; percentiles=percentiles)

# Block average wind direction (CSAT sensors only)
wd1_time, wd1_avg = gen.block_average_winddir(wd1.time, wd1.α, block_length)
wd2_time, wd2_avg = gen.block_average_winddir(wd2.time, wd2.α, block_length)
wd3_time, wd3_avg = gen.block_average_winddir(wd3.time, wd3.α, block_length)
wd4_time, wd4_avg = gen.block_average_winddir(wd4.time, wd4.α, block_length)

# Create figure with 5x2 subplots
fig, axes = PyPlot.subplots(5, 2, figsize=(16, 12), sharex=true)

# Color scheme
cmap = PyPlot.get_cmap("tab10")
colors = [cmap(0), cmap(1), cmap(2), cmap(3)]

# Panel data for wind speed (left column): (ax, time, mean, lower, upper, title, color)
ws_panels = [
    (axes[1, 1], ws1_time, ws1_mean, ws1_lo, ws1_hi, "Tower 1 IRG", colors[1]),
    (axes[2, 1], ws2_time, ws2_mean, ws2_lo, ws2_hi, "Tower 1 CSAT", colors[2]),
    (axes[3, 1], ws3_time, ws3_mean, ws3_lo, ws3_hi, "Tower 2 IRG", colors[3]),
    (axes[4, 1], ws4_time, ws4_mean, ws4_lo, ws4_hi, "Tower 2 CSAT", colors[4]),
]

# Panel data for sonic temperature (right column)
T_panels = [
    (axes[1, 2], T1_time, T1_mean, T1_lo, T1_hi, "Tower 1 IRG", colors[1]),
    (axes[2, 2], T2_time, T2_mean, T2_lo, T2_hi, "Tower 1 CSAT", colors[2]),
    (axes[3, 2], T3_time, T3_mean, T3_lo, T3_hi, "Tower 2 IRG", colors[3]),
    (axes[4, 2], T4_time, T4_mean, T4_lo, T4_hi, "Tower 2 CSAT", colors[4]),
]

# Plot wind speed panels (left column)
for (ax, t, m, lo, hi, title, col) in ws_panels
    ax.fill_between(t, lo, hi, alpha=0.3, color=col, 
                    label="$(percentiles[1])th-$(percentiles[2])th pctl.")
    ax.plot(t, m, color=col, linewidth=1.0, label="mean")
    ax.set_ylabel(L"U~\mathrm{[m~s^{-1}]}")
    ax.set_title(title, fontsize=10, loc="left")
    ax.grid(true, alpha=0.3)
    ax.legend(fontsize=7, loc="upper right")
end

# Plot sonic temperature panels (right column)
for (ax, t, m, lo, hi, title, col) in T_panels
    ax.fill_between(t, lo, hi, alpha=0.3, color=col,
                    label="$(percentiles[1])th-$(percentiles[2])th pctl.")
    ax.plot(t, m, color=col, linewidth=1.0, label="mean")
    ax.set_ylabel(L"T_{sonic}~\mathrm{[^\circ C]}")
    ax.set_title(title, fontsize=10, loc="left")
    ax.grid(true, alpha=0.3)
    ax.legend(fontsize=7, loc="upper right")
end

# Plot wind direction panel (bottom left)
ax_wd = axes[5, 1]
ax_wd.scatter(wd1_time, wd1_avg, s=8, color=colors[1], label="Tower 1 IRG", alpha=0.7)
ax_wd.scatter(wd2_time, wd2_avg, s=8, color=colors[2], label="Tower 1 CSAT", alpha=0.7)
ax_wd.scatter(wd3_time, wd3_avg, s=8, color=colors[3], label="Tower 2 IRG", alpha=0.7)
ax_wd.scatter(wd4_time, wd4_avg, s=8, color=colors[4], label="Tower 2 CSAT", alpha=0.7)
ax_wd.set_ylabel(L"\alpha~\mathrm{[°]}")
ax_wd.set_ylim(0, 360)
ax_wd.set_yticks([0, 90, 180, 270, 360])
ax_wd.set_title("Wind Direction (in sensor coordinates!)", fontsize=10, loc="left")
ax_wd.grid(true, alpha=0.3)
ax_wd.legend(fontsize=7, loc="upper right")
ax_wd.set_xlabel("Time")

# Plot water vapor panel (bottom right)
ax_h2o = axes[5, 2]
ax_h2o.fill_between(h2o1_time, h2o1_lo, h2o1_hi, alpha=0.3, color=colors[1])
ax_h2o.plot(h2o1_time, h2o1_mean, color=colors[1], linewidth=1.0, label="Tower 1 IRG")
ax_h2o.fill_between(h2o3_time, h2o3_lo, h2o3_hi, alpha=0.3, color=colors[3])
ax_h2o.plot(h2o3_time, h2o3_mean, color=colors[3], linewidth=1.0, label="Tower 2 IRG")
ax_h2o.set_ylabel(L"\rho_{H_2O}~\mathrm{[g~m^{-3}]}")
ax_h2o.set_title("Water Vapor Density", fontsize=10, loc="left")
ax_h2o.grid(true, alpha=0.3)
ax_h2o.legend(fontsize=7, loc="upper right")
ax_h2o.set_xlabel("Time")

# Auto-format dates
fig.autofmt_xdate()

# Add overall title
#fig.suptitle("Wind Speed, Temperature, and Humidity ($(block_length) averages)", 
#             fontsize=12, fontweight="bold")

PyPlot.tight_layout()

##
# Save the figure
output_folder = stationcfg.plot_dir(station_config, "wind_temperature")
#PyPlot.savefig(joinpath(output_folder, "$(station_file_stem).pdf"), bbox_inches="tight")
