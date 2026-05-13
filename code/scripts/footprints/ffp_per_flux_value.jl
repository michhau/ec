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
station_file_stem = stationcfg.station_file_stem(station_config)

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
contour_indices = Int.(stationcfg.optional_key(
    station_config,
    fill(-1, length(names)),
    "footprint",
    "contour_indices",
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
# Plot the per-block footprints as an animation.
plot_footprints = false
if plot_footprints
    fileorthomosaic = String(stationcfg.require_key(station_config, "footprint", "orthomosaic"))
    fluxloc = stationcfg.toml_matrix(stationcfg.require_key(station_config, "footprint", "fluxloc"); T=Float64)
    bgextend_m = Float64.(stationcfg.require_key(station_config, "footprint", "bgextend_m"))
    bgextend_pxl = Float64.(stationcfg.require_key(station_config, "footprint", "bgextend_pxl"))
    figorigin = Float64.(stationcfg.require_key(station_config, "footprint", "figorigin"))
    footprint_labels = @isdefined(instr_labels) ? instr_labels : stationcfg.station_labels(station_config)

    footprint_sets = [ffp_blocks_all[flux] for flux in fluxes]
    input_sets = [ffp_inputs_all[flux] for flux in fluxes]
    animation_output_dir = stationcfg.plot_dir(station_config, "footprints", "blocks", station_file_stem)
    animation_output_file = joinpath(animation_output_dir, "$(station_file_stem)_ffp_blocks.mp4")

    animation_file = ffp_block.save_block_footprint_animation(
        footprint_sets,
        input_sets,
        fileorthomosaic,
        fluxloc,
        bgextend_m,
        bgextend_pxl,
        figorigin,
        footprint_labels,
        contour_indices,
        animation_output_file;
        station_label=station_label,
    )
    println("Saved block footprint animation to ", animation_file)
end
##################################################
# Combine with a drone overview image to classify

mpimg = pyimport("matplotlib.image")
classesimg = mpimg.imread(String(station_config["block_footprints"]["classes_file"]))

using Images
classes_img = load(String(station_config["block_footprints"]["classes_file"]))

#create classes matrix
"""
    create_classes_array(img::Matrix{RGB{N0f8}}, )

Create a class array with surface classes from an input image. Criterions need to be hard coded!
"""
function create_classes_array(img::Matrix{RGB{N0f8}})::Array{Union{Missing, Bool}}
    classes_array = Array{Union{Missing, Bool}}(undef, size(img))#fill(false, size(img))
    #missing: do not use
    #true: ice
    #false: lead
    red_vals = round.(Int, float.(red.(img))*255)
    green_vals = round.(Int, float.(green.(img))*255)
    blue_vals = round.(Int, float.(blue.(img))*255)
    for i in eachindex(classes_array)
        if red_vals[i] > 222 && green_vals[i] > 222 && blue_vals[i] > 222
            classes_array[i] = true
        elseif red_vals[i] == green_vals[i] == blue_vals[i] == 0
            classes_array[i] = false
        else
            classes_array[i] = missing
        end
    end
    return classes_array
end

#read same as for footprints. Change, if pictures are not the same!
classes_array = create_classes_array(classes_img)
classes_fluxloc_pxl = stationcfg.toml_matrix(stationcfg.require_key(station_config, "footprint", "fluxloc"); T=Float64)
classes_extend_m = Float64.(stationcfg.require_key(station_config, "footprint", "bgextend_m"))
classes_extend_pxl = Float64.(stationcfg.require_key(station_config, "footprint", "bgextend_pxl"))
classes_gorigin = Float64.(stationcfg.require_key(station_config, "footprint", "figorigin"))

footprint_sets = [ffp_blocks_all[flux] for flux in fluxes]

for ix in eachindex(footprint_sets)
    pixel_center = [classes_fluxloc_pxl[ix, 2], classes_fluxloc_pxl[ix, 1]]
    
end

a = copy(classes_array)
replace!(a, missing => true)
#=
#todo similar to the footprint plotting
- norm footprint (multiply by (x_ci[2]-x_ci[1])^2 to obtain absolute value, otherwise 1/area)
- extract the footprint weighted ice fraction
- extract the footprint weighted lead fraction
=#
PyPlot.figure()
PyPlot.imshow(classesimg)
PyPlot.gcf()