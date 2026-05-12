######################################################
###      CALCULATE FOOTPRINTS PER BLOCK FLUX       ###
###            author: Michi Haugeneder            ###
######################################################
#=
Calculate flux footprints according to
N. Kljun et al. (2015) A simple two-dimensional para-
meterisation for Flux Footprint Prediction (FFP), gmd
need to run turb_fluxes before for eg. Obukhov-length

The flux DataFrames fx1-4 contain centered moving averages.
This script samples every ra1-4-th value, starting at the
first full moving-average window center, so consecutive
sampled values represent non-overlapping block averages.
=#
using DataFrames

importdir = joinpath(@__DIR__, "..", "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "kljun_ffp.jl"))
include(joinpath(importdir, "src", "ffp_block.jl"))
if !@isdefined stationcfg
    include(joinpath(importdir, "src", "station_config.jl"))
    import .stationcfg
end
import .turb
import .ffp_block

if !@isdefined station_config
    error("Run load_data.jl before ffp_per_flux_value.jl so station_config is available.")
end
station_label = stationcfg.station_label(station_config)

#variables
names = [:evaldf1, :evaldf2, :evaldf3, :evaldf4]
meas_heights = Float64.(stationcfg.optional_key(station_config, heights, "footprint", "measurement_heights"))
pbl_height = 1000.0
fluxes = [:fx1, :fx2, :fx3, :fx4]
reyavg_periods = [:ra1, :ra2, :ra3, :ra4]
wd = [:wd1, :wd2, :wd3, :wd4] #wind directions
block_flux_names = [:block_fluxes1, :block_fluxes2, :block_fluxes3, :block_fluxes4]
footprint_input_names = [:ffp_inputs1, :ffp_inputs2, :ffp_inputs3, :ffp_inputs4]
outnames = [:ffp_blocks1, :ffp_blocks2, :ffp_blocks3, :ffp_blocks4]
nrelemgrid = 1000
wind_direction_offsets = Float64.(stationcfg.optional_key(
    station_config,
    fill(0.0, length(names)),
    "footprint",
    "wind_direction_offsets",
))

block_fluxes_all = Dict{Symbol, DataFrame}()
ffp_inputs_all = Dict{Symbol, DataFrame}()
ffp_blocks_all = Dict{Symbol, Vector{Dict{String, Any}}}()

#optional input
rs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] #levels for later plotting/analysis
rslayer = false #measurement within roughness sublayer (theory not working properly)
crop = false #crop output to maximum defined rs (max 0.9)

println()
println("Calculating block fluxes and per-block footprints for station ", station_label)

footprints = Vector{Dict{String, Any}}(undef, 0)

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

    footprints = ffp_block.calculate_footprints(inputs, meas_heights[ix], rs, rslayer, nrelemgrid, crop)

    block_fluxes_all[fluxes[ix]] = block_flux
    ffp_inputs_all[fluxes[ix]] = inputs
    ffp_blocks_all[fluxes[ix]] = footprints

    @eval $(block_flux_names[ix]) = $block_flux
    @eval $(footprint_input_names[ix]) = $inputs
    @eval $(outnames[ix]) = $footprints
end

###############################################
# Plotting can be added after the desired per-block display/storage format is defined.
plot_footprints = false
if plot_footprints
    @info("Plotting of per-block footprints is deferred until the desired display/storage format is defined.")
end
println("------------D-O-N-E---------------")
println()
