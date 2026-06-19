######################################################
###       SYNTHESIS OF BLOCK-DIFFERENCE FITS       ###
######################################################
#=
Create a three-panel scatter plot that pools block-level flux-difference
relationships from multiple stations. Each station is shown with a different
color, while the linear fit in each panel treats every valid block observation
equally across all stations.
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings
import PyPlot, CSV

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

# Variables to convert fluxes to W/m^2.
rho_air = 1.2 # kg m^{-3}
c_p = 1004 # J kg^{-1} K^{-1}
L_v = 2500e3 # J kg^{-1}

#######################################################
station_names = ["1b_2", "2b"]
save_figs = true
save_csv = true
#######################################################

function load_synthesis_station(station_name::AbstractString)
    station_config = stationcfg.load_station_config(station_name)
    station_label = stationcfg.station_label(station_config)
    station_file_stem = stationcfg.station_file_stem(station_config)
    surface_type = String.(stationcfg.require_key(station_config, "surface_type"))
    heights = Float64.(stationcfg.require_key(station_config, "heights"))
    surface_feature = String(stationcfg.require_key(station_config, "block_footprints", "feature"))
    block_data = ffp_block.read_block_fluxes_netcdf(station_name)

    subplot_order = block_analyze.read_subplot_order(station_config, "subplot_order",
        [:fx4, :fx2, :fx3, :fx1])
    latent_subplot_order = block_analyze.read_subplot_order(station_config, "latent_subplot_order",
        [:fx3, :fx1])
    difference_specs = block_analyze.standard_difference_specs(
        subplot_order, latent_subplot_order, rho_air * c_p, L_v * 1e-3)

    wT_limits = Tuple(Float64.(stationcfg.optional_key(
        station_config, [-10.0, 10.0], "block_footprints", "wT_limits_diff_timeseries")))
    wq_limits = Tuple(Float64.(stationcfg.optional_key(
        station_config, [-10.0, 10.0], "block_footprints", "wq_limits_diff_timeseries")))

    return (
        station_name=String(station_name),
        station_config=station_config,
        station_label=station_label,
        station_file_stem=station_file_stem,
        surface_type=surface_type,
        heights=heights,
        surface_feature=surface_feature,
        block_data=block_data,
        difference_specs=difference_specs,
        wT_limits=wT_limits,
        wq_limits=wq_limits,
    )
end

function combined_ylim(station_contexts, flux_column::Symbol)
    limits = flux_column == :wT ?
        [context.wT_limits for context in station_contexts] :
        [context.wq_limits for context in station_contexts]
    return (minimum([limit[1] for limit in limits]), maximum([limit[2] for limit in limits]))
end

station_contexts = [load_synthesis_station(station_name) for station_name in station_names]
feature_labels = unique([context.surface_feature for context in station_contexts])
length(feature_labels) == 1 || error(
    "Synthesis stations must use the same block_footprints.feature; got $(join(feature_labels, ", "))."
)
surface_feature = feature_labels[1]
feature_file_label = block_analyze.feature_file_label(surface_feature)
comparison_label = join([context.station_label for context in station_contexts], " + ")
comparison_file_stem = join([context.station_file_stem for context in station_contexts], "_")

fig_diff_scatter, axs_diff_scatter = PyPlot.subplots(1, 3, figsize=(15, 4.2), sharex=true)
fig_diff_scatter.suptitle(
    "$(comparison_label) - Heat flux difference vs $(surface_feature) fraction difference"
)

station_colors = ["C0", "C3", "C2", "C4", "C5", "C6"]
linear_fit_params = []
station_count_rows = []

for panel_ix in eachindex(station_contexts[1].difference_specs)
    ax = vec(axs_diff_scatter)[panel_ix]
    reference_spec = station_contexts[1].difference_specs[panel_ix]
    _, flux_column, _, flux_kind, flux_ylabel = reference_spec
    all_x = Float64[]
    all_y = Float64[]

    for (station_ix, context) in enumerate(station_contexts)
        flux_names, station_flux_column, conversion_factor, station_flux_kind, _ =
            context.difference_specs[panel_ix]
        station_flux_column == flux_column || error(
            "Panel $panel_ix mixes flux columns $(flux_column) and $(station_flux_column)."
        )
        station_flux_kind == flux_kind || error(
            "Panel $panel_ix mixes flux kinds $(flux_kind) and $(station_flux_kind)."
        )

        x, y = block_analyze.feature_flux_difference_xy(
            context.block_data, flux_names, flux_column, conversion_factor)
        append!(all_x, x)
        append!(all_y, y)

        color = station_colors[mod1(station_ix, length(station_colors))]
        ax.scatter(x, y; s=12, alpha=0.65, color=color, edgecolors="none",
            zorder=2, label=context.station_label)

        push!(station_count_rows, (
            panel_ix=panel_ix,
            station=context.station_name,
            station_label=context.station_label,
            feature=surface_feature,
            flux_variable=String(flux_column),
            flux_label=String(flux_kind),
            difference_label=block_analyze.sensor_difference_label(
                flux_names, context.heights, context.surface_type),
            n_valid=length(x),
        ))
    end

    fit = block_analyze.fit_feature_flux_difference(all_x, all_y)
    block_analyze.plot_feature_flux_difference_fit!(ax, fit; color="black",
        label="pooled linear fit", ci_label="pooled 95% CI", ci_alpha=0.14, linewidth=1.2)

    fit_summary = block_analyze.feature_flux_difference_fit_summary(
        fit, surface_feature, flux_column, flux_kind, "pooled panel $(panel_ix)")
    push!(linear_fit_params, (;
        panel_ix=panel_ix,
        stations=join(station_names, "+"),
        fit_summary...,
    ))

    ax.axhline(0, color="grey", linewidth=0.8, alpha=0.55)
    ax.axvline(0, color="grey", linewidth=0.8, alpha=0.55)
    ax.set_title("$(flux_kind): panel $(panel_ix)")
    ax.set_xlabel("$(surface_feature) fraction difference")
    ax.set_ylabel(flux_ylabel)
    ax.set_ylim(combined_ylim(station_contexts, flux_column))
    ax.set_axisbelow(true)
    ax.grid(alpha=0.9)
    ax.text(0.03, 0.97, "n=$(fit.n_valid), R^2=$(round(fit.r2, digits=2))";
        transform=ax.transAxes, ha="left", va="top", fontsize=9)
    ax.legend(loc="best")
end

PyPlot.tight_layout(rect=[0, 0, 1, 0.91])

output_folder = stationcfg.plot_dir(
    station_contexts[1].station_config, "footprints", "blocks", "comparison")
mkpath(output_folder)

if save_csv
    CSV.write(joinpath(output_folder,
        "$(comparison_file_stem)_block_flux_$(feature_file_label)_difference_linear_fit.csv"),
        DataFrame(linear_fit_params))
    CSV.write(joinpath(output_folder,
        "$(comparison_file_stem)_block_flux_$(feature_file_label)_difference_counts.csv"),
        DataFrame(station_count_rows))
end

if save_figs
    fig_diff_scatter.savefig(joinpath(output_folder,
        "$(comparison_file_stem)_block_flux_$(feature_file_label)_difference_scatter.pdf"),
        bbox_inches="tight")
end
