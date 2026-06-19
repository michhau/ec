######################################################
###              ANALYSE THE BLOCK DATA            ###
###  (FLUXES, FOOTPRINTS, SURFACE FRACTIONS)       ###
###            author: Michi Haugeneder            ###
######################################################
#=
Run load_data.jl, turb_fluxes.jl, ffp_per_flux_value.jl
in this order before this script to generate
netcdf_files containing the block data. If they already
exist, no need to run the previous scripts, simply load
them.
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
station_name = "3c"
save_figs = true
save_csv = true
#######################################################
#read station specifics
station_config = stationcfg.load_station_config(station_name)
station_label = stationcfg.station_label(station_config)
station_file_stem = stationcfg.station_file_stem(station_config)
station_id = String(stationcfg.require_key(station_config, "id"))
surface_type = String.(stationcfg.require_key(station_config, "surface_type"))
heights = Float64.(stationcfg.require_key(station_config, "heights"))
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
output_folder = stationcfg.plot_dir(station_config, "footprints", "blocks", station_name)
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
##
#######################################################
#plot time series of fluxes and feature (lead, pond, ridge) fraction
time_series_columns = [
    (subplot_order[1], subplot_order[3]),
    (subplot_order[2], subplot_order[4]),
]

fig_ts, axs_ts = PyPlot.subplots(2, 2, figsize=(12, 7.5), sharex=true)
fig_ts.suptitle("$(station_label) - Block fluxes and $(surface_feature) fraction")
for (col, flux_names) in enumerate(time_series_columns)
    axs_ts[1, col].set_title(block_analyze.column_title(flux_names, surface_type))

    block_analyze.plot_block_timeseries_panel!(axs_ts[1, col], block_data,
        flux_names, :wT, ρ_air * c_p, L"\overline{w'T_s'}~\mathrm{[W~m^{-2}]}",
        heights;
        feature_label=surface_feature,
        show_ylabel=col == 1, show_right_ylabel=col == 2, show_legend=true)
    axs_ts[1, col].set_ylim(wT_limits)

    block_analyze.plot_block_timeseries_panel!(axs_ts[2, col], block_data,
        flux_names, :wq, L_v * 1e-3, L"\overline{w'q'}~\mathrm{[W~m^{-2}]}",
        heights;
        feature_label=surface_feature,
        show_xlabel=true, show_ylabel=col == 1, show_right_ylabel=col == 2,
        show_legend=false)
    axs_ts[2, col].set_ylim(wq_limits)
    axs_ts[2, col].xaxis_date()
end
fig_ts.autofmt_xdate()
PyPlot.tight_layout(rect=[0, 0, 1, 0.95])
if save_figs
    fig_ts.savefig(joinpath(output_folder, "$(station_file_stem)_block_flux_$(feature_file_label)_timeseries.pdf"), bbox_inches="tight")
end
#######################################################
#plot correlation fluxes vs. feature fraction
fig_wT, axs_wT = PyPlot.subplots(2, 2, figsize=(9, 7), sharex=true, sharey=true)
fig_wT.suptitle("$(station_label) - Sensible heat flux vs $(surface_feature) fraction")
wT_panel_specs = [
    (axs_wT[1, 1], subplot_order[1], false, true),
    (axs_wT[1, 2], subplot_order[2], false, false),
    (axs_wT[2, 1], subplot_order[3], true, true),
    (axs_wT[2, 2], subplot_order[4], true, false),
]
for (ax, flux_name, show_xlabel, show_ylabel) in wT_panel_specs
    block_analyze.plot_feature_flux_panel!(ax, block_data, flux_name, :wT, ρ_air * c_p,
        L"\overline{w'T_s'}~\mathrm{[W~m^{-2}]}", wT_limits, surface_type, heights;
        color="black", show_xlabel=show_xlabel, show_ylabel=show_ylabel,
        feature_label=surface_feature)
end
PyPlot.tight_layout(rect=[0, 0, 1, 0.95])
if save_figs
    fig_wT.savefig(joinpath(output_folder, "$(station_file_stem)_block_wT_$(feature_file_label)_fraction.pdf"),bbox_inches="tight")
end
fig_wq, axs_wq = PyPlot.subplots(1, 2, figsize=(9, 3.8), sharex=true, sharey=true)

for (ix, (ax, flux_name)) in enumerate(zip(vec(axs_wq), latent_subplot_order))
    block_analyze.plot_feature_flux_panel!(ax, block_data, flux_name, :wq, L_v * 1e-3,
        L"\overline{w'q'}~\mathrm{[W~m^{-2}]}", wq_limits, surface_type, heights;
        color="C0", show_ylabel=ix == 1, feature_label=surface_feature)
end
fig_wq.suptitle("$(station_label) - Latent heat flux vs $(surface_feature) fraction")
PyPlot.tight_layout(rect=[0, 0, 1, 0.92])
if save_figs
    fig_wq.savefig(joinpath(output_folder, "$(station_file_stem)_block_wq_$(feature_file_label)_fraction.pdf"), bbox_inches="tight")
end


#######################################################
##
#colored correlation scatter plots
colored_scatter_variables = ["time", "u", "wT", "wq", "wind_dir", "T"]
colored_correlation_cmap = "viridis"
colored_correlation_specs = [
    (subplot_order[1], :wT, ρ_air * c_p, "sensible",
        L"\overline{w'T_s'}~\mathrm{[W~m^{-2}]}", wT_limits),
    (subplot_order[2], :wT, ρ_air * c_p, "sensible",
        L"\overline{w'T_s'}~\mathrm{[W~m^{-2}]}", wT_limits),
    (subplot_order[3], :wT, ρ_air * c_p, "sensible",
        L"\overline{w'T_s'}~\mathrm{[W~m^{-2}]}", wT_limits),
    (subplot_order[4], :wT, ρ_air * c_p, "sensible",
        L"\overline{w'T_s'}~\mathrm{[W~m^{-2}]}", wT_limits),
    (latent_subplot_order[1], :wq, L_v * 1e-3, "latent",
        L"\overline{w'q'}~\mathrm{[W~m^{-2}]}", wq_limits),
    (latent_subplot_order[2], :wq, L_v * 1e-3, "latent",
        L"\overline{w'q'}~\mathrm{[W~m^{-2}]}", wq_limits),
]
colored_correlation_flux_names = [
    flux_name for (flux_name, _, _, _, _, _) in colored_correlation_specs
]

for color_variable_raw in colored_scatter_variables
    colored_scatter_variable = block_analyze.normalize_colored_scatter_variable(color_variable_raw)
    colored_scatter_origin = colored_scatter_variable == "time" ?
        block_analyze.colored_flux_time_origin(block_data, colored_correlation_flux_names) : nothing
    colored_panel_data = [
        block_analyze.colored_correlation_data(
            block_data,
            flux_name,
            flux_column,
            conversion_factor,
            colored_scatter_variable,
            colored_scatter_origin,
        )
        for (flux_name, flux_column, conversion_factor, _, _, _) in colored_correlation_specs
    ]

    all_color_values = Float64[]
    for panel_data in colored_panel_data
        append!(all_color_values, panel_data.color_values)
    end

    if isempty(all_color_values)
        @warn(
            "No colored correlation scatter plot saved for color variable \"$(colored_scatter_variable)\". " *
            block_analyze.colored_scatter_prerequisite_note(colored_scatter_variable)
        )
    else
        color_min = minimum(all_color_values)
        color_max = maximum(all_color_values)
        if colored_scatter_variable == "wind_dir"
            color_min = 0.0
            color_max = 360.0
        end
        if color_min == color_max
            color_min -= 0.5
            color_max += 0.5
        end

        fig_corr_colored, axs_corr_colored =
            PyPlot.subplots(3, 2, figsize=(12, 10.2), sharex=true)
        fig_corr_colored.suptitle(
            "$(station_label) - Heat flux vs $(surface_feature) fraction, colored by $(colored_scatter_variable)"
        )
        corr_axes = [
            axs_corr_colored[1, 1], axs_corr_colored[1, 2],
            axs_corr_colored[2, 1], axs_corr_colored[2, 2],
            axs_corr_colored[3, 1], axs_corr_colored[3, 2],
        ]
        scatter_handle = nothing

        for panel_ix in eachindex(colored_correlation_specs)
            ax = corr_axes[panel_ix]
            flux_name, _, _, flux_kind, flux_ylabel, flux_ylim =
                colored_correlation_specs[panel_ix]
            panel_data = colored_panel_data[panel_ix]
            row_ix = cld(panel_ix, 2)
            col_ix = isodd(panel_ix) ? 1 : 2

            ax.set_title("$(flux_kind): $(block_analyze.panel_title(flux_name, surface_type, heights))")
            ax.set_xlabel(row_ix == 3 ? "$(surface_feature) fraction" : "")
            ax.set_ylabel(col_ix == 1 ? flux_ylabel : "")
            ax.set_xlim(0, 1)
            ax.set_ylim(flux_ylim)
            ax.set_axisbelow(true)
            ax.grid(alpha=0.9)

            if panel_data.available && !isempty(panel_data.color_values)
                scatter_handle = ax.scatter(
                    panel_data.feature_fraction,
                    panel_data.flux,
                    c=panel_data.color_values,
                    cmap=colored_correlation_cmap,
                    vmin=color_min,
                    vmax=color_max,
                    s=12,
                    alpha=0.75,
                    edgecolors="none",
                    zorder=2,
                )
            else
                message = panel_data.available ?
                    "No finite values for color variable \"$(colored_scatter_variable)\"." :
                    panel_data.message
                @warn(message)
                ax.text(
                    0.5,
                    0.5,
                    message,
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=8,
                    wrap=true,
                )
            end
        end

        PyPlot.tight_layout(rect=[0, 0, 0.90, 0.94])
        colorbar_ax = fig_corr_colored.add_axes([0.92, 0.16, 0.015, 0.70])
        colorbar = fig_corr_colored.colorbar(scatter_handle, cax=colorbar_ax)
        colorbar.set_label(block_analyze.colored_scatter_colorbar_label(
            colored_scatter_variable, colored_scatter_origin))

        if save_figs
            fig_corr_colored.savefig(joinpath(output_folder,
                "$(station_file_stem)_block_flux_$(feature_file_label)_correlation_scatter_colored_by_$(colored_scatter_variable).pdf"),
                bbox_inches="tight")
        end
    end
end
##
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

fig_diff_ts, axs_diff_ts = PyPlot.subplots(1, 3, figsize=(15, 4.2), sharex=true)
fig_diff_ts.suptitle("$(station_label) - Heat flux and $(surface_feature) fraction differences")
for (ax, (flux_names, flux_column, conversion_factor, flux_kind, flux_ylabel)) in
        zip(vec(axs_diff_ts), difference_time_series_specs)
    block_analyze.plot_block_difference_timeseries_panel!(ax, block_data, flux_names,
        flux_column, conversion_factor, flux_kind, flux_ylabel, heights,
        surface_type, wT_limits_diff_timeseries, wq_limits_diff_timeseries;
        feature_label=surface_feature)
    ax.xaxis_date()
end
fig_diff_ts.autofmt_xdate()
PyPlot.tight_layout(rect=[0, 0, 1, 0.91])
if save_figs
    fig_diff_ts.savefig(joinpath(output_folder,
        "$(station_file_stem)_block_flux_$(feature_file_label)_difference_timeseries.pdf"),
        bbox_inches="tight")
end

#######################################################
#scatter plot
fig_diff_scatter, axs_diff_scatter = PyPlot.subplots(1, 3, figsize=(15, 4.2), sharex=true)
fig_diff_scatter.suptitle("$(station_label) - Heat flux difference vs $(surface_feature) fraction difference")
for (ax, (flux_names, flux_column, conversion_factor, flux_kind, flux_ylabel)) in
        zip(vec(axs_diff_scatter), difference_time_series_specs)
    flux_ylim = flux_column == :wT ? wT_limits_diff_timeseries : wq_limits_diff_timeseries
    color = flux_column == :wT ? "black" : "C0"
    block_analyze.plot_feature_flux_difference_panel!(ax, block_data, flux_names,
        flux_column, conversion_factor, flux_kind, flux_ylabel, flux_ylim,
        surface_type, heights; color=color, feature_label=surface_feature)
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

if save_figs
    fig_diff_scatter.savefig(joinpath(output_folder,
    "$(station_file_stem)_block_flux_$(feature_file_label)_difference_scatter.pdf"), bbox_inches="tight")
end
##

#######################################################
##
#scatter plot with color-coded points
for colored_scatter_variable in ["time", "u", "wT", "wq", "wind_dir", "T"]
    #colored_scatter_variable = "wq" # choose "u", "T", "wind_dir", "wT", "wq", or "time"
    colored_scatter_cmap = "viridis"

    colored_scatter_variable = block_analyze.normalize_colored_scatter_variable(colored_scatter_variable)
    colored_scatter_origin = colored_scatter_variable == "time" ?
        block_analyze.colored_scatter_time_origin(block_data, difference_time_series_specs) : nothing
    colored_panel_data = [
        block_analyze.colored_difference_data(
            block_data,
            flux_names,
            flux_column,
            conversion_factor,
            colored_scatter_variable,
            colored_scatter_origin,
        )
        for (flux_names, flux_column, conversion_factor, _, _) in difference_time_series_specs
    ]

    all_color_values = Float64[]
    for panel_data in colored_panel_data
        append!(all_color_values, panel_data.color_values)
    end

    if isempty(all_color_values)
        @warn(
            "No colored difference scatter plot saved for color variable \"$(colored_scatter_variable)\". " *
            block_analyze.colored_scatter_prerequisite_note(colored_scatter_variable)
        )
    else
        color_min = minimum(all_color_values)
        color_max = maximum(all_color_values)
        if colored_scatter_variable == "wind_dir"
            color_min = 0.0
            color_max = 360.0
        end
        if color_min == color_max
            color_min -= 0.5
            color_max += 0.5
        end

        fig_diff_scatter_colored, axs_diff_scatter_colored =
            PyPlot.subplots(1, 3, figsize=(15.8, 4.2), sharex=true)
        fig_diff_scatter_colored.suptitle(
            "$(station_label) - Heat flux difference vs $(surface_feature) fraction difference, colored by $(colored_scatter_variable)"
        )
        scatter_handle = nothing

        for panel_ix in eachindex(difference_time_series_specs)
            ax = vec(axs_diff_scatter_colored)[panel_ix]
            flux_names, flux_column, _, flux_kind, flux_ylabel = difference_time_series_specs[panel_ix]
            panel_data = colored_panel_data[panel_ix]
            flux_ylim = flux_column == :wT ? wT_limits_diff_timeseries : wq_limits_diff_timeseries

            ax.axhline(0, color="grey", linewidth=0.8, alpha=0.55)
            ax.axvline(0, color="grey", linewidth=0.8, alpha=0.55)
            ax.set_title("$(flux_kind): $(block_analyze.sensor_difference_label(flux_names, heights, surface_type))")
            ax.set_xlabel("$(surface_feature) fraction difference")
            ax.set_ylabel(flux_ylabel)
            ax.set_ylim(flux_ylim)
            ax.set_axisbelow(true)
            ax.grid(alpha=0.9)

            if panel_data.available && !isempty(panel_data.color_values)
                scatter_handle = ax.scatter(
                    panel_data.feature_difference,
                    panel_data.flux_difference,
                    c=panel_data.color_values,
                    cmap=colored_scatter_cmap,
                    vmin=color_min,
                    vmax=color_max,
                    s=12,
                    alpha=0.75,
                    edgecolors="none",
                    zorder=2,
                )

                prediction = predictions[panel_ix]
                ax.fill_between(
                    prediction.x_grid,
                    prediction.predictions_confidence.lower,
                    prediction.predictions_confidence.upper;
                    color="0.45",
                    alpha=0.16,
                    linewidth=0,
                    label="95% CI",
                )
                ax.plot(
                    prediction.x_grid,
                    prediction.predictions_confidence.prediction;
                    color="0.2",
                    linewidth=1.0,
                )
                ax.legend(loc="best")
            else
                message = panel_data.available ?
                    "No finite values for color variable \"$(colored_scatter_variable)\"." :
                    panel_data.message
                @warn(message)
                ax.text(
                    0.5,
                    0.5,
                    message,
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=8,
                    wrap=true,
                )
            end
        end

        PyPlot.tight_layout(rect=[0, 0, 0.90, 0.91])
        colorbar_ax = fig_diff_scatter_colored.add_axes([0.92, 0.18, 0.015, 0.66])
        colorbar = fig_diff_scatter_colored.colorbar(scatter_handle, cax=colorbar_ax)
        colorbar.set_label(block_analyze.colored_scatter_colorbar_label(colored_scatter_variable, colored_scatter_origin))

        if save_figs
            fig_diff_scatter_colored.savefig(joinpath(output_folder,
                "$(station_file_stem)_block_flux_$(feature_file_label)_difference_scatter_colored_by_$(colored_scatter_variable).pdf"),
                bbox_inches="tight")
        end
    end
end
##
