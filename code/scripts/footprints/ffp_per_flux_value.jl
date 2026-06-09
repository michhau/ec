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
using Dates, DataFrames, PyCall, PyPlot

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
station_config = stationcfg.load_station_config(station_name)
station_label = stationcfg.station_label(station_config)
station_file_stem = stationcfg.station_file_stem(station_config)
station_id = String(stationcfg.require_key(station_config, "id"))

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

    block_flux[!, :u] = [
        ffp_block.nanmean(
            sqrt.(ecdata.u[start_indices[j]:end_indices[j]] .^ 2 .+
                  ecdata.v[start_indices[j]:end_indices[j]] .^ 2)
        )
        for j in eachindex(start_indices)
    ]
    block_flux[!, :wind_dir] = [
        turb.mean_winddir(wd_tmp[ecdata.time[start_indices[j]] .<= wd_tmp.time .<
            ecdata.time[end_indices[j]], Symbol("α")])
        for j in eachindex(start_indices)
    ]

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
plot_footprints = true
if plot_footprints
    fileorthomosaic = String(stationcfg.require_key(station_config, "footprint", "orthomosaic"))
    fluxloc = stationcfg.toml_matrix(stationcfg.require_key(station_config, "footprint", "fluxloc"); T=Float64)
    bgextend_m_config = stationcfg.optional_key(station_config, nothing, "footprint", "bgextend_m")
    bgextend_pxl_config = stationcfg.optional_key(station_config, nothing, "footprint", "bgextend_pxl")
    bgextend_m = isnothing(bgextend_m_config) ? nothing : Float64.(bgextend_m_config)
    bgextend_pxl = isnothing(bgextend_pxl_config) ? nothing : Float64.(bgextend_pxl_config)
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
using Images
classes_file = String(station_config["block_footprints"]["classes_file"])
classes_reference_file = String(stationcfg.require_key(station_config, "footprint", "orthomosaic"))
classes_scale_file = isfile(ffp_block.footprint_world_file(classes_file)) ? classes_file : classes_reference_file
classes_img = load(classes_file)
classes_img_size = size(classes_img)
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

#=
PyPlot.pygui(true)
check_plot_array = replace(classes_array, missing => NaN)
PyPlot.figure()
PyPlot.imshow(check_plot_array)
=#

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

block_fluxes_output_file = ffp_block.save_block_fluxes_netcdf(block_fluxes_all, station_id)
println("Saved block fluxes to ", block_fluxes_output_file)

function class_background_geometry(fluxloc::AbstractMatrix,
        bgextend_m::Union{AbstractVector, Nothing},
        bgextend_pxl::Union{AbstractVector, Nothing}, figorigin::AbstractVector;
        image_file::Union{AbstractString, Nothing}=nothing, image_size=nothing)
    return ffp_block.footprint_background_geometry(
        fluxloc,
        bgextend_m,
        bgextend_pxl,
        figorigin;
        image_file=image_file,
        image_size=image_size,
    )
end

function format_fraction_percent(value)
    return ismissing(value) ? "missing" : "$(round(value * 100, digits=1))%"
end

function class_fraction_frame_text(fraction_sets::AbstractVector{<:AbstractDataFrame},
        labels::AbstractVector, feature_name::String, frame_index::Integer)
    lines = String[]
    for ix in eachindex(fraction_sets)
        fractions = fraction_sets[ix]
        if frame_index <= nrow(fractions)
            ice = format_fraction_percent(fractions.ice_fraction[frame_index])
            feature = format_fraction_percent(fractions.feature_fraction[frame_index])
            classified = format_fraction_percent(fractions.classified_weight_fraction[frame_index])
            outside = format_fraction_percent(fractions.outside_weight_fraction[frame_index])
            missing_class = format_fraction_percent(fractions.missing_class_weight_fraction[frame_index])
            push!(lines, "$(labels[ix]): ice=$(ice), $(feature_name)=$(feature), classified=$(classified), outside=$(outside), missing=$(missing_class)")
        else
            push!(lines, "$(labels[ix]): no footprint")
        end
    end
    return join(lines, "\n")
end

function save_class_fraction_animation(footprint_sets::AbstractVector,
        fraction_sets::AbstractVector{<:AbstractDataFrame}, class_image,
        class_img_size::Tuple{Int64, Int64}, feature::String,
        fluxloc::AbstractMatrix, bgextend_m::Union{AbstractVector, Nothing},
        bgextend_pxl::Union{AbstractVector, Nothing}, figorigin::AbstractVector,
        labels::AbstractVector, contour_indices::AbstractVector{<:Integer},
        output_file::AbstractString; station_label::AbstractString="",
        class_scale_file::Union{AbstractString, Nothing}=nothing,
        interval::Integer=250, fps::Integer=4, dpi::Integer=150, figsize=(10, 8))

    length(footprint_sets) == length(fraction_sets) || error("footprint_sets and fraction_sets must have the same length.")
    length(footprint_sets) == length(labels) || error("footprint_sets and labels must have the same length.")
    length(footprint_sets) == length(contour_indices) || error("footprint_sets and contour_indices must have the same length.")
    nframes = maximum(length.(footprint_sets))
    nframes > 0 || error("No footprint blocks available for animation.")

    animation = pyimport("matplotlib.animation")
    fluxloc_final, bgextend_final = class_background_geometry(
        fluxloc,
        bgextend_m,
        bgextend_pxl,
        figorigin;
        image_file=class_scale_file,
        image_size=class_img_size,
    )

    fig = PyPlot.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.imshow(class_image, extent=bgextend_final)
    ax.set_xlabel("meter")
    ax.set_ylabel("meter")

    ctab10 = PyPlot.cm.tab10
    lines = Any[]
    for ix in eachindex(footprint_sets)
        ax.plot(fluxloc_final[ix, 2], fluxloc_final[ix, 1], ".", color=ctab10(ix - 1), ms=15)
        line = ax.plot(Float64[], Float64[], color=ctab10(ix - 1), label=labels[ix])[1]
        push!(lines, line)
    end
    ax.legend(loc="lower right")

    fraction_text = ax.text(
        0.02,
        0.98,
        "",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox=Dict("facecolor" => "white", "alpha" => 0.75, "edgecolor" => "none"),
    )

    function update(frame)
        j = Int(frame) + 1
        for ix in eachindex(footprint_sets)
            if j <= length(footprint_sets[ix])
                x, y = ffp_block.footprint_contour_xy(footprint_sets[ix][j], contour_indices[ix])
                lines[ix].set_data(x .+ fluxloc_final[ix, 2], y .+ fluxloc_final[ix, 1])
            else
                lines[ix].set_data(Float64[], Float64[])
            end
        end

        fraction_text.set_text(class_fraction_frame_text(fraction_sets, labels, feature, j))
        timestamp = ffp_block.first_frame_time(fraction_sets, j)
        title_time = isnothing(timestamp) ? "" : " - $(Dates.format(timestamp, "yyyy-mm-dd HH:MM:SS"))"
        ax.set_title("Station $(station_label) class footprints - block $(j)/$(nframes)$(title_time)")
        return vcat(lines, [fraction_text])
    end

    update(0)
    PyPlot.tight_layout()
    mkpath(dirname(output_file))

    anim = animation.FuncAnimation(fig, update, frames=nframes, interval=interval, blit=false)
    saved_file = String(output_file)
    try
        anim.save(saved_file, writer="ffmpeg", fps=fps, dpi=dpi)
    catch err
        if endswith(lowercase(saved_file), ".mp4")
            fallback_file = replace(saved_file, r"\.mp4$"i => ".gif")
            @warn("Could not save mp4 animation with ffmpeg. Falling back to gif.", exception=(err, catch_backtrace()), output=fallback_file)
            anim.save(fallback_file, writer="pillow", fps=fps, dpi=dpi)
            saved_file = fallback_file
        else
            rethrow()
        end
    finally
        PyPlot.close(fig)
    end

    return saved_file
end

plot_class_fraction_animation = true
if plot_class_fraction_animation
    footprint_sets = [ffp_blocks_all[flux] for flux in fluxes]
    fraction_sets = [footprint_class_fractions_all[flux] for flux in fluxes]
    footprint_labels = @isdefined(instr_labels) ? instr_labels : stationcfg.station_labels(station_config)
    class_animation_output_dir = stationcfg.plot_dir(station_config, "footprints", "blocks", station_file_stem)
    class_animation_output_file = joinpath(class_animation_output_dir, "$(station_file_stem)_ffp_blocks_classes.mp4")

    class_animation_file = save_class_fraction_animation(
        footprint_sets,
        fraction_sets,
        classes_img_pyplot,
        classes_img_size,
        surface_feature,
        classes_fluxloc_pxl,
        classes_extend_m,
        classes_extend_pxl,
        classes_figorigin_pxl,
        footprint_labels,
        contour_indices,
        class_animation_output_file;
        station_label=station_label,
        class_scale_file=classes_scale_file,
    )
    println("Saved class footprint animation to ", class_animation_file)
end

ffp_block.save_block_fluxes_netcdf(block_fluxes_all, station_label)
