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

datapath = "/home/haugened/Documents/data/CONTRASTS/EC_offline_preproc/cut/"
#datapath = "/home/michi/Documents/slf/CONTRASTS25/data/processed/preproc/"
#tjkpath = "/home/haugened/Documents/data/tjk/tjk_data.csv"
importdir = joinpath(@__DIR__, "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .gen
PyPlot.pygui(true)

######################################################
#variables

#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
#kaijo_period_file = joinpath(datapath, "kaijo", "21_kaijo_periods.txt")
#kaijo_outfile_stam = joinpath(datapath, "kaijo", "kaijo")
#outfile_stam = joinpath(datapath, "WindSonic_Duerrboden_2005xx/tjk_sonic_200504-200514")
#tower_outfile_stam = joinpath(datapath, "tower", "preproc")
#ventair_stam = joinpath(datapath, "tower", "vent_air")
nanfile_stam = "/home/haugened/Documents/data/CONTRASTS/nan_periods"
nanfiles = joinpath.(nanfile_stam, ["t1irg_nan.csv", "t1csat_nan.csv", "t2irg_nan.csv", "t2csat_nan.csv"])

######################################################

println()
println("-----------S-T-A-R-T-------------")

######################################################
###            LOADING & PREPROCESSING             ###
######################################################
#select data and measurement period to be evaluated
evalstart = DateTime(2024, 07, 01, 10, 15, 00)
evalend   = DateTime(2027, 07, 10, 16, 00, 00)
#evalend = evalstart + Day(10)

evaldf1 = turb.readturbasnetcdf(joinpath(datapath, "1a_t1_irg_proc_cut.nc"), evalstart, evalend)
evaldf2 = turb.readturbasnetcdf(joinpath(datapath, "1a_t1_csat_proc_cut.nc"), evalstart, evalend)
evaldf3 = turb.readturbasnetcdf(joinpath(datapath, "1a_t2_irg_proc_cut.nc"), evalstart, evalend)
evaldf4 = turb.readturbasnetcdf(joinpath(datapath, "1a_t2_csat_proc_cut.nc"), evalstart, evalend)
#evaldf5 = turb.readturbasnetcdf(string(kaijo_outfile_stam, ".nc"), evalstart, evalend)
#evaldf6 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "tjkdf.nc"), evalstart, evalend)

#apply NaN-mask to T1 & T2 when repositioned (for DR)
#turb.repositionnanmask!(evaldf1)
#turb.repositionnanmask!(evaldf2)
#turb.repositionnanmask!(evaldf3)
#turb.repositionnanmask!(evaldf4)

#show statistics about missing or NaN-data
turb.printmissstats(evaldf1)
turb.printmissstats(evaldf2)
turb.printmissstats(evaldf3)
turb.printmissstats(evaldf4)
#turb.printmissstats(evaldf5)
#turb.printmissstats(evaldf6)

#interpolate missing
evaldf1 = turb.interpolatemissing(evaldf1);
evaldf2 = turb.interpolatemissing(evaldf2);
evaldf3 = turb.interpolatemissing(evaldf3);
evaldf4 = turb.interpolatemissing(evaldf4);
#try
#    evaldf5 = turb.interpolatemissing(evaldf5)
#catch LoadError
#    @warn("Kaijo DataFrame empty")
#end
#evaldf6 = turb.interpolatemissing(evaldf6)

wd1 = turb.winddir(evaldf1)
wd2 = turb.winddir(evaldf2)
wd3 = turb.winddir(evaldf3)
wd4 = turb.winddir(evaldf4)

#double rotation
turb.drdf!(evaldf1, periodwise=false)
turb.drdf!(evaldf2, periodwise=false)
turb.drdf!(evaldf3, periodwise=false)
turb.drdf!(evaldf4, periodwise=false)
#turb.drdf!(evaldf5)
#turb.drdf!(evaldf6, periodwise=false)

#CONTRASTS: apply manual nan-mask
manual_period_1 = turb.read_nantimes_csv(nanfiles[1]);
manual_period_2 = turb.read_nantimes_csv(nanfiles[2]);
manual_period_3 = turb.read_nantimes_csv(nanfiles[3]);
manual_period_4 = turb.read_nantimes_csv(nanfiles[4]);
isbad_quality_1 = turb.manual_nanmask(manual_period_1, evaldf1.time);
isbad_quality_2 = turb.manual_nanmask(manual_period_2, evaldf2.time);
isbad_quality_3 = turb.manual_nanmask(manual_period_3, evaldf3.time);
isbad_quality_4 = turb.manual_nanmask(manual_period_4, evaldf4.time);
evaldf1[isbad_quality_1, ["u", "v", "w", "T"]] .= NaN;
evaldf2[isbad_quality_2, ["u", "v", "w", "T"]] .= NaN;
evaldf3[isbad_quality_3, ["u", "v", "w", "T"]] .= NaN;
evaldf4[isbad_quality_4, ["u", "v", "w", "T"]] .= NaN;

#turb.saveturbasnetcdf(evaldf5, "/home/haugened/Documents/openfoam/duerr_les/scripts/src/kaijotmp.nc")

######################################################
###               LOADING SLOW DATA                ###
######################################################

#tjkmeteodata = turb.csvtodataframe(tjkpath)
#tjkmeteodata = tjkmeteodata[evalstart.<=tjkmeteodata.time.<=evalend, :]

#t2vent = turb.readturbasnetcdf(joinpath(ventair_stam, "t2.nc"), evalstart, evalend)
#t3vent = turb.readturbasnetcdf(joinpath(ventair_stam, "t3.nc"), evalstart, evalend)

#remove unrealistic values
#t2vent.vent_air_temp[.!(-30 .< replace!(t2vent.vent_air_temp, missing => NaN) .< 30)] .= NaN
#t3vent.vent_air_temp[.!(-30 .< replace!(t3vent.vent_air_temp, missing => NaN) .< 30)] .= NaN

######################################################
###                   PLOTTING                     ###
######################################################
#=
##
#plot u,v,w,T components to visually inspect data
evaldf = evaldf1
skip_n = 200  # Plot every nth data point

fig, (ax1, ax2, ax3, ax4) = PyPlot.subplots(4, 1, figsize=(10, 12), sharex=true)

# Plot u component
ax1.plot(evaldf.time[1:skip_n:end], evaldf.u[1:skip_n:end])
ax1.set_ylabel(L"u~\mathrm{[m~s^{-1}]}")
ax1.grid(true)

# Plot v component
ax2.plot(evaldf.time[1:skip_n:end], evaldf.v[1:skip_n:end])
ax2.set_ylabel(L"v~\mathrm{[m~s^{-1}]}")
ax2.grid(true)

# Plot w component
ax3.plot(evaldf.time[1:skip_n:end], evaldf.w[1:skip_n:end])
ax3.set_ylabel(L"w~\mathrm{[m~s^{-1}]}")
ax3.grid(true)

# Plot T component
ax4.plot(evaldf.time[1:skip_n:end], evaldf.T[1:skip_n:end])
ax4.set_ylabel(L"T~\mathrm{[K]}")
ax4.grid(true)

# Adjust layout to prevent overlapping
fig.tight_layout()
##
=#
#######################################################
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

    fig.suptitle("Station 3d", fontsize=14, fontweight="bold")
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
wd1 = turb.winddir(evaldf1);
wd2 = turb.winddir(evaldf2);
wd3 = turb.winddir(evaldf3);
wd4 = turb.winddir(evaldf4);

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

output_folder = "/home/haugened/Documents/data/CONTRASTS/plots/wind_roses"
PyPlot.savefig(joinpath(output_folder, "3d.pdf"), bbox_inches="tight")
##
#########################################################
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
output_folder = "/home/haugened/Documents/data/CONTRASTS/plots/wind_temperature/"
PyPlot.savefig(joinpath(output_folder, "3d.pdf"), bbox_inches="tight")
#PyPlot.savefig(joinpath(output_folder, "3a.png"), dpi=150, bbox_inches="tight")

##

#######################################################
#=
#######################################################
#plot histograms to determine outliers

quantity1 = filter(!isnan, skipmissing(evaldf1.T))
quantity2 = filter(!isnan, skipmissing(evaldf2.T))
quantity3 = filter(!isnan, skipmissing(evaldf3.T))
quantity4 = filter(!isnan, skipmissing(evaldf4.T))
quantity5 = filter(!isnan, skipmissing(evaldf5.w))
binmin = -8
binmax = 8

##
comap = PyPlot.get_cmap("tab10");
hist = PyPlot.figure()
axh = hist.add_subplot(111)
axh.set_title("Kajio measurements (raw)")
axh.set_xlabel(L"w~\mathrm{[m~s^{-1}]}")
axh.grid()
PyPlot.yscale("log")
#axh.hist(vec(quantity1), bins=collect(binmin:binmax), density=true, label="T1IRG", alpha=0.5, color=comap(0))
#axh.hist(vec(quantity2), bins=collect(binmin:binmax), density=true, label="T2IRG", alpha=0.5, color=comap(1))
#axh.hist(vec(quantity3), bins=collect(binmin:binmax), density=true, label="T2LCSAT", alpha=0.5, color=comap(2))
#axh.hist(vec(quantity4), bins=collect(binmin:binmax), density=true, label="T2UCSAT", alpha=0.5, color=comap(3))
axh.hist(vec(quantity5), bins=collect(binmin:0.01:binmax), density=true, label="KAIJO", color=comap(4))
axh.legend()
##

##
#plot single parameter (e.g. u)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf4)
sifig = PyPlot.figure()
siax = sifig.add_subplot(111)
#siax.set_title(L"\mathrm{1~min~gliding~average}", fontsize=12)
siax.set_xlabel("03.06.2021")
siax.set_ylabel("winddir (w/o DR)")
#siax.set_ylabel(L"\overline{u}~\mathrm{[m~s^{-1}]}")
#siax.plot(evaldf1.time, gen.movingaverage(evaldf1.u,20*60), label="T1IRG")
siax.plot(evaldf2.time, gen.movingaverage(turb.simplewinddir.(evaldf2.u, evaldf2.v), 20*600), label="T2IRG, 10min avg")
#siax.plot(evaldf3.time, gen.movingaverage(evaldf3.u,20*60), label="T2LCSAT")
siax.plot(evaldf4.time, gen.movingaverage(turb.simplewinddir.(evaldf4.u, evaldf4.v), 20*600), label="T2UCSAT 10min avg")
#siax.plot(evaldf5.time, gen.movingaverage(evaldf5.u,20*60), label="KAIJO")
majorlocator = pydates.HourLocator(interval=1)
minorlocator = pydates.MinuteLocator([15,30,45])
siax.xaxis.set_major_locator(majorlocator)
siax.xaxis.set_minor_locator(minorlocator)
date_format = pydates.DateFormatter("%H:%M")
siax.xaxis.set_major_formatter(date_format)
#sifig.autofmt_xdate()
siax.grid()
siax.legend()
##
=#
#=
#plot time series with moving average
fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.plot(evaldf2.time[1:40:end], gen.movingaverage(evaldf2.u,40)[1:40:end], label="1m")
ax.plot(evaldf3.time[1:40:end], gen.movingaverage(evaldf3.u,40)[1:40:end], label="2m")
ax.plot(evaldf4.time[1:40:end], gen.movingaverage(evaldf4.u,40)[1:40:end], label="3m")
ax.plot(evaldf5.time[1:40:end], gen.movingaverage(evaldf5.u,40)[1:40:end], label="0.3m")
ax.plot(evaldf6.time[1:40:end], gen.movingaverage(evaldf6.u,40)[1:40:end], label="5m")
ax.set_ylabel(L"u~\mathrm{[m~s^{-1}]}")
ax.grid()
ax.legend()
=#