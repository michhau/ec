######################################################
###              ANALYSE THE BLOCK DATA            ###
###  (FLUXES, FOOTPRINTS, SURFACE FRACTIONS)       ###
###            author: Michi Haugeneder            ###
######################################################
#=
Run load_data.jl, turb_fluxes.jl, ffp_per_flux_value.jl
in this order before this script to generate
netcdf_files containing the block data
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")

importdir = joinpath(@__DIR__, "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "kljun_ffp.jl"))
include(joinpath(importdir, "src", "ffp_block.jl"))
include(joinpath(importdir, "src", "block_analysis_src.jl"))
if !@isdefined stationcfg
    include(joinpath(importdir, "src", "station_config.jl"))
    import .stationcfg
end
import .ffp_block
import .block_analyze

PyPlot.pygui(true)

#variables to convert fluxes to W/m^2
ρ_air = 1.2 #kg m^{-3}
c_p = 1004 #J kg^{-1} K^{-1}
L_v = 2500e3 #J kg^{-1} (approx @0°C)

#######################################################
station_name = "2a"
#######################################################
#read station specifics
station_config = stationcfg.load_station_config(station_name)
station_label = stationcfg.station_label(station_config)
station_file_stem = stationcfg.station_file_stem(station_config)
station_id = String(stationcfg.require_key(station_config, "id"))
surface_type = String.(stationcfg.require_key(station_config, "surface_type"))
heights = Float64.(stationcfg.require_key(station_config, "heights"))
#######################################################
#load the block file
block_data = ffp_block.read_block_fluxes_netcdf(station_name)

#order of the subplots
subplot_order = [:fx4, :fx2, :fx3, :fx1] #left top, right top, left bottom, right bottom
latent_subplot_order = [:fx3, :fx1]
#######################################################
###                    PLOTS                        ###
#######################################################
wT_limits = Tuple(Float64.(stationcfg.optional_key(station_config, [-20.0, 20.0], "plot", "wT_limits")))
wq_limits = Tuple(Float64.(stationcfg.optional_key(station_config, [-5.0, 5.0], "plot", "wq_limits")))
output_folder = stationcfg.plot_dir(station_config, "footprints", "blocks", station_name)
mkpath(output_folder)

#######################################################
#calculate correlations and export to .csv
correlation_specs = vcat(
    [(flux_name, :wT, "sensible_heat_flux", ρ_air * c_p) for flux_name in subplot_order],
    [(flux_name, :wq, "latent_heat_flux", L_v * 1e-3) for flux_name in latent_subplot_order],
)
flux_lead_correlations = block_analyze.lead_flux_correlations(
    correlation_specs, block_data, station_name, surface_type, heights)
CSV.write(joinpath(output_folder, "$(station_file_stem)_block_flux_lead_correlations.csv"),
    flux_lead_correlations)

#######################################################
#plot time series of fluxes and lead fraction
time_series_columns = [
    (subplot_order[1], subplot_order[3]),
    (subplot_order[2], subplot_order[4]),
]

fig_ts, axs_ts = PyPlot.subplots(2, 2, figsize=(12, 7.5), sharex=true)
fig_ts.suptitle("$(station_label) - Block fluxes and lead fraction")
for (col, flux_names) in enumerate(time_series_columns)
    axs_ts[1, col].set_title(block_analyze.column_title(flux_names, surface_type))

    block_analyze.plot_block_timeseries_panel!(axs_ts[1, col], block_data,
        flux_names, :wT, ρ_air * c_p, L"\overline{w'T_s'}~\mathrm{[W~m^{-2}]}",
        heights;
        show_ylabel=col == 1, show_right_ylabel=col == 2, show_legend=true)
    axs_ts[1, col].set_ylim(wT_limits)

    block_analyze.plot_block_timeseries_panel!(axs_ts[2, col], block_data,
        flux_names, :wq, L_v * 1e-3, L"\overline{w'q'}~\mathrm{[W~m^{-2}]}",
        heights;
        show_xlabel=true, show_ylabel=col == 1, show_right_ylabel=col == 2,
        show_legend=false)
    axs_ts[2, col].set_ylim(wq_limits)
    axs_ts[2, col].xaxis_date()
end
fig_ts.autofmt_xdate()
PyPlot.tight_layout(rect=[0, 0, 1, 0.95])
fig_ts.savefig(joinpath(output_folder, "$(station_file_stem)_block_flux_lead_timeseries.pdf"),
    bbox_inches="tight")

#######################################################
#plot correlation fluxes vs. lead fraction
fig_wT, axs_wT = PyPlot.subplots(2, 2, figsize=(9, 7), sharex=true, sharey=true)
fig_wT.suptitle("$(station_label) - Sensible heat flux vs lead fraction")
wT_panel_specs = [
    (axs_wT[1, 1], subplot_order[1], false, true),
    (axs_wT[1, 2], subplot_order[2], false, false),
    (axs_wT[2, 1], subplot_order[3], true, true),
    (axs_wT[2, 2], subplot_order[4], true, false),
]
for (ax, flux_name, show_xlabel, show_ylabel) in wT_panel_specs
    block_analyze.plot_lead_flux_panel!(ax, block_data, flux_name, :wT, ρ_air * c_p,
        L"\overline{w'T_s'}~\mathrm{[W~m^{-2}]}", wT_limits, surface_type, heights;
        color="black", show_xlabel=show_xlabel, show_ylabel=show_ylabel)
end
PyPlot.tight_layout(rect=[0, 0, 1, 0.95])
fig_wT.savefig(joinpath(output_folder, "$(station_file_stem)_block_wT_lead_fraction.pdf"),
    bbox_inches="tight")

fig_wq, axs_wq = PyPlot.subplots(1, 2, figsize=(9, 3.8), sharex=true, sharey=true)

for (ix, (ax, flux_name)) in enumerate(zip(vec(axs_wq), latent_subplot_order))
    block_analyze.plot_lead_flux_panel!(ax, block_data, flux_name, :wq, L_v * 1e-3,
        L"\overline{w'q'}~\mathrm{[W~m^{-2}]}", wq_limits, surface_type, heights;
        color="C0", show_ylabel=ix == 1)
end
fig_wq.suptitle("$(station_label) - Latent heat flux vs lead fraction")
PyPlot.tight_layout(rect=[0, 0, 1, 0.92])
fig_wq.savefig(joinpath(output_folder, "$(station_file_stem)_block_wq_lead_fraction.pdf"),
    bbox_inches="tight")

#######################################################
#######################################################
#plot flux DIFFERENCES and lead fraction DIFFERENCES

#######################################################
#time series
wT_limits_diff_timeseries = Tuple(Float64.(stationcfg.optional_key(station_config, [-10.0, 10.0], "block_footprints", "wT_limits_diff_timeseries")))
wq_limits_diff_timeseries = Tuple(Float64.(stationcfg.optional_key(station_config, [-10.0, 10.0], "block_footprints", "wq_limits_diff_timeseries")))

difference_time_series_specs = [
    ((subplot_order[1], subplot_order[2]), :wT, ρ_air * c_p, "sensible",
        L"\Delta H~\mathrm{[W~m^{-2}]}"),
    ((subplot_order[3], subplot_order[4]), :wT, ρ_air * c_p, "sensible",
        L"\Delta H~\mathrm{[W~m^{-2}]}"),
    ((latent_subplot_order[1], latent_subplot_order[2]), :wq, L_v * 1e-3, "latent",
        L"\Delta L_E~\mathrm{[W~m^{-2}]}"),
]

flux_lead_difference_correlations = block_analyze.lead_flux_difference_correlations(
    difference_time_series_specs, block_data, station_name, surface_type, heights)
CSV.write(joinpath(output_folder, "$(station_file_stem)_block_flux_lead_difference_correlations.csv"),
    flux_lead_difference_correlations)

fig_diff_ts, axs_diff_ts = PyPlot.subplots(1, 3, figsize=(15, 4.2), sharex=true)
fig_diff_ts.suptitle("$(station_label) - Heat flux and lead fraction differences")
for (ax, (flux_names, flux_column, conversion_factor, flux_kind, flux_ylabel)) in
        zip(vec(axs_diff_ts), difference_time_series_specs)
    block_analyze.plot_block_difference_timeseries_panel!(ax, block_data, flux_names,
        flux_column, conversion_factor, flux_kind, flux_ylabel, heights,
        surface_type, wT_limits_diff_timeseries, wq_limits_diff_timeseries)
    ax.xaxis_date()
end
fig_diff_ts.autofmt_xdate()
PyPlot.tight_layout(rect=[0, 0, 1, 0.91])
fig_diff_ts.savefig(joinpath(output_folder,
    "$(station_file_stem)_block_flux_difference_timeseries.pdf"),
    bbox_inches="tight")

#######################################################
#scatter plot
fig_diff_scatter, axs_diff_scatter = PyPlot.subplots(1, 3, figsize=(15, 4.2), sharex=true)
fig_diff_scatter.suptitle("$(station_label) - Heat flux difference vs lead fraction difference")
for (ax, (flux_names, flux_column, conversion_factor, flux_kind, flux_ylabel)) in
        zip(vec(axs_diff_scatter), difference_time_series_specs)
    flux_ylim = flux_column == :wT ? wT_limits_diff_timeseries : wq_limits_diff_timeseries
    color = flux_column == :wT ? "black" : "C0"
    block_analyze.plot_lead_flux_difference_panel!(ax, block_data, flux_names,
        flux_column, conversion_factor, flux_kind, flux_ylabel, flux_ylim,
        surface_type, heights; color=color)
end
PyPlot.tight_layout(rect=[0, 0, 1, 0.91])
fig_diff_scatter.savefig(joinpath(output_folder,
    "$(station_file_stem)_block_flux_lead_difference_scatter.pdf"),
    bbox_inches="tight")
