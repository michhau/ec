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
if !@isdefined stationcfg
    include(joinpath(importdir, "src", "station_config.jl"))
    import .stationcfg
end
import .ffp_block

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
function valid_number(value)
    return !ismissing(value) && value isa Number && isfinite(Float64(value))
end

function flux_index(flux_name::Symbol)
    name = String(flux_name)
    startswith(name, "fx") || error("Expected flux names like :fx1, got :$(name).")
    return parse(Int, name[3:end])
end

function panel_title(flux_name::Symbol)
    ix = flux_index(flux_name)
    return "$(surface_type[ix]), $(heights[ix]) m"
end

function scatter_xy(data::AbstractDataFrame, flux_column::Symbol, conversion_factor::Real)
    lead_fraction = data[!, :footprint_lead_fraction]
    flux = data[!, flux_column]
    valid = [
        valid_number(lead_fraction[ix]) && valid_number(flux[ix])
        #&& 0.0 <= Float64(lead_fraction[ix]) <= 1.0
        for ix in eachindex(lead_fraction)
    ]
    return Float64.(lead_fraction[valid]), Float64.(flux[valid]) .* conversion_factor
end

function add_axes_diagonal!(ax)
    # Lead fraction and flux have different units; draw this as an axes-space guide.
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="grey",
        alpha=0.35, linewidth=1.2, zorder=1)
end

function plot_lead_flux_panel!(ax, flux_name::Symbol, flux_column::Symbol,
        conversion_factor::Real, ylabel, ylim; color="C0",
        show_xlabel=true, show_ylabel=true)
    data = block_data[flux_name]
    lead_fraction, flux = scatter_xy(data, flux_column, conversion_factor)

    #add_axes_diagonal!(ax)
    ax.scatter(lead_fraction, flux, s=12, alpha=0.65, color=color,
        edgecolors="none", zorder=2)
    ax.set_title(panel_title(flux_name))
    ax.set_xlabel(show_xlabel ? "Lead fraction" : "")
    ax.set_ylabel(show_ylabel ? ylabel : "")
    ax.set_xlim(0, 1)
    ax.set_ylim(ylim)
    ax.set_axisbelow(true)
    ax.grid(alpha=0.9)

    return ax
end

wT_limits = Tuple(Float64.(stationcfg.optional_key(station_config, [-20.0, 20.0], "plot", "wT_limits")))
wq_limits = Tuple(Float64.(stationcfg.optional_key(station_config, [-5.0, 5.0], "plot", "wq_limits")))
output_folder = stationcfg.plot_dir(station_config, "footprints", "blocks", station_name)
mkpath(output_folder)

fig_wT, axs_wT = PyPlot.subplots(2, 2, figsize=(9, 7), sharex=true, sharey=true)
fig_wT.suptitle("$(station_label) - Sensible heat flux vs lead fraction")
wT_panel_specs = [
    (axs_wT[1, 1], subplot_order[1], false, true),
    (axs_wT[1, 2], subplot_order[2], false, false),
    (axs_wT[2, 1], subplot_order[3], true, true),
    (axs_wT[2, 2], subplot_order[4], true, false),
]
for (ax, flux_name, show_xlabel, show_ylabel) in wT_panel_specs
    plot_lead_flux_panel!(ax, flux_name, :wT, ρ_air * c_p,
        L"\overline{w'T_s'}~\mathrm{[W~m^{-2}]}", wT_limits;
        color="black", show_xlabel=show_xlabel, show_ylabel=show_ylabel)
end
PyPlot.tight_layout(rect=[0, 0, 1, 0.95])
fig_wT.savefig(joinpath(output_folder, "$(station_file_stem)_block_wT_lead_fraction.pdf"),
    bbox_inches="tight")

fig_wq, axs_wq = PyPlot.subplots(1, 2, figsize=(9, 3.8), sharex=true, sharey=true)

for (ix, (ax, flux_name)) in enumerate(zip(vec(axs_wq), latent_subplot_order))
    plot_lead_flux_panel!(ax, flux_name, :wq, L_v * 1e-3,
        L"\overline{w'q'}~\mathrm{[W~m^{-2}]}", wq_limits;
        color="C0", show_ylabel=ix == 1)
end
fig_wq.suptitle("$(station_label) - Latent heat flux vs lead fraction")
PyPlot.tight_layout(rect=[0, 0, 1, 0.92])
fig_wq.savefig(joinpath(output_folder, "$(station_file_stem)_block_wq_lead_fraction.pdf"),
    bbox_inches="tight")

