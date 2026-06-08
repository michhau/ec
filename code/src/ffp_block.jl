######################################################
###       HELPERS FOR PER-BLOCK FOOTPRINTS         ###
###            author: Michi Haugeneder            ###
######################################################
module ffp_block

using Dates, DataFrames, Statistics, ProgressMeter, PyCall, NCDatasets
using Images
import PyPlot

import ..kljun
import ..turb
include(joinpath(@__DIR__, "footprint_plotting.jl"))

export block_footprint_inputs, calculate_footprints, extract_block_fluxes,
    footprint_dict, full_block_indices, nan_footprint, nanmean, nanstd,
    rows_per_period, footprint_background_geometry, footprint_meterperpxl,
    footprint_world_file, footprint_contour_xy,
    save_block_footprint_animation, block_fluxes_netcdf_path,
    save_block_fluxes_netcdf, read_block_fluxes_netcdf#, valid_footprint_inputs

const DEFAULT_BLOCK_FLUXES_DIR = "/home/haugened/Documents/data/CONTRASTS/block_fluxes"
const NETCDF_NAME_SEPARATOR = "\t"

function nanmean(values)
    valid = [value for value in values if !ismissing(value) && isa(value, Number) && !isnan(value)]
    return isempty(valid) ? NaN : mean(valid)
end

function nanstd(values)
    valid = [value for value in values if !ismissing(value) && isa(value, Number) && !isnan(value)]
    return length(valid) <= 1 ? NaN : std(valid)
end

function rows_per_period(time::AbstractVector{DateTime}, peri::Period)
    length(time) >= 2 || error("At least two timestamps are required to determine the sampling interval.")
    timestep = time[2] - time[1]
    return maximum([1, round(Int, Millisecond(peri) / Millisecond(timestep))])
end

function full_block_indices(time::AbstractVector{DateTime}, peri::Period)
    block_rows = rows_per_period(time, peri)
    back_delta = div(block_rows, 2)
    forward_delta = isodd(block_rows) ? div(block_rows, 2) : div(block_rows, 2) - 1

    first_center = 1 + back_delta
    last_center = length(time) - forward_delta
    if last_center < first_center
        return Int[], Int[], Int[], block_rows
    end

    center_indices = collect(first_center:block_rows:last_center)
    start_indices = center_indices .- back_delta
    end_indices = center_indices .+ forward_delta
    return center_indices, start_indices, end_indices, block_rows
end

function extract_block_fluxes(fluxdata::DataFrame, peri::Period)
    center_indices, start_indices, end_indices, block_rows = full_block_indices(fluxdata.time, peri)
    out = DataFrame(fluxdata[center_indices, :])
    insertcols!(out, 2,
        :block_start => fluxdata.time[start_indices],
        :block_end => fluxdata.time[end_indices],
        :block_rows => fill(block_rows, length(center_indices)),
    )
    return out, center_indices, start_indices, end_indices
end

function nan_footprint(rs, reason::AbstractString, inputs::Dict{String, Any})
    nr = isnothing(rs) ? 1 : length(rs)
    return Dict{String, Any}(
        "x_ci_max" => NaN,
        "x_ci" => [NaN],
        "f_ci" => [NaN],
        "x_2d" => fill(NaN, 1, 1),
        "y_2d" => fill(NaN, 1, 1),
        "f_2d" => fill(NaN, 1, 1),
        "rs" => rs,
        "fr" => fill(NaN, nr),
        "xr" => fill(NaN, 1, nr),
        "yr" => fill(NaN, 1, nr),
        "flag_err" => true,
        "error" => reason,
        "inputs" => inputs,
    )
end

function footprint_dict(result, inputs::Dict{String, Any})
    x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d, rs_out, fr, xr, yr, flag_err = result
    return Dict{String, Any}(
        "x_ci_max" => x_ci_max,
        "x_ci" => x_ci,
        "f_ci" => f_ci,
        "x_2d" => x_2d,
        "y_2d" => y_2d,
        "f_2d" => f_2d,
        "rs" => rs_out,
        "fr" => fr,
        "xr" => xr,
        "yr" => yr,
        "flag_err" => flag_err,
        "error" => flag_err ? "kljun.ffp returned flag_err=true" : "",
        "inputs" => inputs,
    )
end

function block_footprint_inputs(ecdata::DataFrame, block_flux::DataFrame, wd_tmp::DataFrame,
        start_indices::Vector{Int}, end_indices::Vector{Int}, zm::Real, pbl_height::Real,
        wind_direction_offset::Real)
    nblocks = nrow(block_flux)
    inputs = DataFrame(
        time=block_flux.time,
        block_start=block_flux.block_start,
        block_end=block_flux.block_end,
        umean=fill(NaN, nblocks),
        h=fill(Float64(pbl_height), nblocks),
        ol=block_flux.L_highfreq,
        sigmav=fill(NaN, nblocks),
        ustar=block_flux.u_star,
        wind_dir=fill(NaN, nblocks),
    )

    alpha_col = Symbol("α")
    for j in 1:nblocks
        six = start_indices[j]
        eix = end_indices[j]

        speed = sqrt.(ecdata.u[six:eix] .^ 2 .+ ecdata.v[six:eix] .^ 2 .+ ecdata.w[six:eix] .^ 2)
        inputs.umean[j] = nanmean(speed)
        inputs.sigmav[j] = nanstd(ecdata.v[six:eix])

        wind_dir_raw = wd_tmp[ecdata.time[six] .<= wd_tmp.time .< ecdata.time[eix], alpha_col]
        wind_dir = turb.mean_winddir(wind_dir_raw)
        inputs.wind_dir[j] = isfinite(wind_dir) ? (wind_dir + wind_direction_offset) % 360 : NaN
    end

    return inputs
end

function calculate_footprints(inputs::DataFrame, zm::Real, rs, rslayer::Bool, nrelemgrid::Int, crop::Bool)
    footprints = Vector{Dict{String, Any}}(undef, nrow(inputs))

    @showprogress "Calculating block footprints..." for j in 1:nrow(inputs)
        input_dict = Dict{String, Any}(
            "time" => inputs.time[j],
            "block_start" => inputs.block_start[j],
            "block_end" => inputs.block_end[j],
            "zm" => zm,
            "umean" => inputs.umean[j],
            "h" => inputs.h[j],
            "ol" => inputs.ol[j],
            "sigmav" => inputs.sigmav[j],
            "ustar" => inputs.ustar[j],
            "wind_dir" => inputs.wind_dir[j],
        )

        required_inputs = (
        input_dict["zm"],
        input_dict["umean"],
        input_dict["h"],
        input_dict["ol"],
        input_dict["sigmav"],
        input_dict["ustar"],
        input_dict["wind_dir"])

        if !all(value -> value isa Real && isfinite(Float64(value)), required_inputs)
            footprints[j] = nan_footprint(rs, "non-finite footprint input", input_dict)
            continue
        end

        try
            result = kljun.ffp(
                zm,
                nothing,
                inputs.umean[j],
                inputs.h[j],
                inputs.ol[j],
                inputs.sigmav[j],
                inputs.ustar[j],
                inputs.wind_dir[j],
                rs,
                rslayer,
                nrelemgrid,
                crop,
            )
            footprints[j] = footprint_dict(result, input_dict)
        catch err
            footprints[j] = nan_footprint(rs, sprint(showerror, err), input_dict)
        end
    end

    return footprints
end

function contour_level_index(contours::AbstractMatrix, contour_offset::Integer)
    nlevels = size(contours, 2)
    contour_index = nlevels + contour_offset
    1 <= contour_index <= nlevels || error(
        "Footprint contour offset $contour_offset is outside available contour range."
    )
    return contour_index
end

function footprint_contour_xy(footprint::AbstractDict, contour_offset::Integer)
    get(footprint, "flag_err", false) && return Float64[], Float64[]
    haskey(footprint, "xr") || return Float64[], Float64[]
    haskey(footprint, "yr") || return Float64[], Float64[]

    xr = footprint["xr"]
    yr = footprint["yr"]
    if xr isa AbstractMatrix && yr isa AbstractMatrix
        contour_index = contour_level_index(xr, contour_offset)
        x = xr[:, contour_index]
        y = yr[:, contour_index]
    elseif xr isa AbstractVector && yr isa AbstractVector
        x = xr
        y = yr
    else
        return Float64[], Float64[]
    end

    finite = isfinite.(x) .& isfinite.(y)
    return Float64.(x[finite]), Float64.(y[finite])
end

function first_frame_time(input_sets::AbstractVector, frame_index::Integer)
    for inputs in input_sets
        if inputs isa DataFrame && frame_index <= nrow(inputs)
            return inputs.time[frame_index]
        end
    end
    return nothing
end

"""
    block_fluxes_netcdf_path(station_name; output_dir=DEFAULT_BLOCK_FLUXES_DIR)

Build the default NetCDF output path for block fluxes.
"""
function block_fluxes_netcdf_path(station_name::AbstractString;
        output_dir::AbstractString=DEFAULT_BLOCK_FLUXES_DIR)
    return joinpath(String(output_dir), "$(station_name).nc")
end

function resolve_block_fluxes_netcdf_path(source::AbstractString;
        input_dir::AbstractString=DEFAULT_BLOCK_FLUXES_DIR)
    source_string = String(source)
    if endswith(lowercase(source_string), ".nc") || isabspath(source_string) || occursin("/", source_string)
        return source_string
    end
    return block_fluxes_netcdf_path(source_string; output_dir=input_dir)
end

function split_netcdf_names(value)
    value_string = String(value)
    isempty(value_string) && return String[]
    return String.(split(value_string, NETCDF_NAME_SEPARATOR))
end

function normalized_netcdf_column(values, column_name::AbstractString)
    data = collect(values)
    eltype(data) !== Any && return data

    valid = collect(skipmissing(data))
    isempty(valid) && return Union{Missing, Float64}[missing for _ in data]

    if all(value -> value isa DateTime, valid)
        return Union{Missing, DateTime}[ismissing(value) ? missing : value for value in data]
    elseif all(value -> value isa Bool, valid)
        return Union{Missing, Bool}[ismissing(value) ? missing : value for value in data]
    elseif all(value -> value isa Integer, valid)
        return Union{Missing, Int64}[ismissing(value) ? missing : Int64(value) for value in data]
    elseif all(value -> value isa Number, valid)
        return Union{Missing, Float64}[ismissing(value) ? missing : Float64(value) for value in data]
    elseif all(value -> value isa AbstractString, valid)
        return Union{Missing, String}[ismissing(value) ? missing : String(value) for value in data]
    end

    error("Cannot save column '$column_name' to NetCDF. Unsupported element type $(eltype(data)).")
end

"""
    save_block_fluxes_netcdf(block_fluxes_all, station_name;
        output_dir=DEFAULT_BLOCK_FLUXES_DIR, deflatelvl=5, overwrite=true)

Save a `Dict{Symbol, DataFrame}` of block fluxes to a station NetCDF file.
Each flux key is stored as a NetCDF group and each DataFrame column as a
variable in that group. Returns the written file path.
"""
function save_block_fluxes_netcdf(block_fluxes_all::AbstractDict, station_name::AbstractString;
        output_dir::AbstractString=DEFAULT_BLOCK_FLUXES_DIR, deflatelvl::Integer=5,
        overwrite::Bool=true)
    target = block_fluxes_netcdf_path(station_name; output_dir=output_dir)
    mkpath(dirname(target))

    if isfile(target)
        overwrite || error("NetCDF file already exists: $target")
        rm(target)
    end

    @info("Saving block fluxes to NetCDF-file", target)
    flux_keys = sort(collect(keys(block_fluxes_all)); by=String)
    ds = NCDataset(target, "c")
    try
        ds.attrib["station_name"] = String(station_name)
        ds.attrib["flux_names"] = join(String.(flux_keys), NETCDF_NAME_SEPARATOR)
        ds.attrib["created"] = string(now())

        for flux_key in flux_keys
            flux_name = String(flux_key)
            data = block_fluxes_all[flux_key]
            data isa AbstractDataFrame || error("Block flux '$flux_name' must be a DataFrame.")

            group = defGroup(ds, flux_name)
            defDim(group, "block", nrow(data))

            column_names = names(data)
            variable_names = ["col_$(lpad(i, 3, '0'))" for i in eachindex(column_names)]
            group.attrib["column_names"] = join(column_names, NETCDF_NAME_SEPARATOR)
            group.attrib["variable_names"] = join(variable_names, NETCDF_NAME_SEPARATOR)

            for (column_name, variable_name) in zip(column_names, variable_names)
                column_data = normalized_netcdf_column(data[!, column_name], column_name)
                variable = defVar(
                    group,
                    variable_name,
                    column_data,
                    ("block",);
                    shuffle=true,
                    deflatelevel=deflatelvl,
                )
                variable.attrib["column_name"] = column_name
            end
        end
    finally
        close(ds)
    end

    return target
end

"""
    read_block_fluxes_netcdf(source; input_dir=DEFAULT_BLOCK_FLUXES_DIR)

Read block fluxes saved with `save_block_fluxes_netcdf`. `source` can be a
station name such as `"2a"` or a full NetCDF file path. Returns
`Dict{Symbol, DataFrame}`.
"""
function read_block_fluxes_netcdf(source::AbstractString;
        input_dir::AbstractString=DEFAULT_BLOCK_FLUXES_DIR)
    path = resolve_block_fluxes_netcdf_path(source; input_dir=input_dir)
    @info("Reading block fluxes from NetCDF-file", path)

    ds = NCDataset(path, "r")
    block_fluxes_all = Dict{Symbol, DataFrame}()
    try
        flux_names = split_netcdf_names(ds.attrib["flux_names"])
        for flux_name in flux_names
            group = ds.group[flux_name]
            column_names = split_netcdf_names(group.attrib["column_names"])
            variable_names = split_netcdf_names(group.attrib["variable_names"])
            length(column_names) == length(variable_names) || error(
                "Column metadata mismatch in NetCDF group '$flux_name'."
            )

            data = DataFrame()
            for (column_name, variable_name) in zip(column_names, variable_names)
                data[!, Symbol(column_name)] = group[variable_name][:]
            end
            block_fluxes_all[Symbol(flux_name)] = data
        end
    finally
        close(ds)
    end

    return block_fluxes_all
end

function save_block_footprint_animation(footprint_sets::AbstractVector, input_sets::AbstractVector,
        orthomosaic_file::AbstractString, fluxloc::AbstractMatrix,
        bgextend_m::Union{AbstractVector, Nothing}, bgextend_pxl::Union{AbstractVector, Nothing},
        figorigin::AbstractVector, labels::AbstractVector,
        contour_indices::AbstractVector{<:Integer}, output_file::AbstractString;
        station_label::AbstractString="", interval::Integer=250, fps::Integer=4,
        dpi::Integer=150, figsize=(10, 8))

    length(footprint_sets) == length(input_sets) || error("footprint_sets and input_sets must have the same length.")
    length(footprint_sets) == length(labels) || error("footprint_sets and labels must have the same length.")
    length(footprint_sets) == length(contour_indices) || error("footprint_sets and contour_indices must have the same length.")
    nframes = maximum(length.(footprint_sets))
    nframes > 0 || error("No footprint blocks available for animation.")

    animation = pyimport("matplotlib.animation")
    mpimg = pyimport("matplotlib.image")

    orthomosaic = mpimg.imread(String(orthomosaic_file))
    orthomosaic_jl = load(String(orthomosaic_file))
    fluxloc_final, bgextend_final = footprint_background_geometry(
        fluxloc,
        bgextend_m,
        bgextend_pxl,
        figorigin;
        image_file=orthomosaic_file,
        image_size=size(orthomosaic_jl),
    )

    fig = PyPlot.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.imshow(orthomosaic, extent=bgextend_final)
    ax.set_xlabel("meter")
    ax.set_ylabel("meter")

    ctab10 = PyPlot.cm.tab10
    lines = PyObject[]
    for ix in eachindex(footprint_sets)
        ax.plot(fluxloc_final[ix, 2], fluxloc_final[ix, 1], ".", color=ctab10(ix - 1), ms=15)
        line = ax.plot(Float64[], Float64[], color=ctab10(ix - 1), label=labels[ix])[1]
        push!(lines, line)
    end
    ax.legend()

    function update(frame)
        j = Int(frame) + 1
        for ix in eachindex(footprint_sets)
            if j <= length(footprint_sets[ix])
                x, y = footprint_contour_xy(footprint_sets[ix][j], contour_indices[ix])
                lines[ix].set_data(x .+ fluxloc_final[ix, 2], y .+ fluxloc_final[ix, 1])
            else
                lines[ix].set_data(Float64[], Float64[])
            end
        end

        timestamp = first_frame_time(input_sets, j)
        title_time = isnothing(timestamp) ? "" : " - $(Dates.format(timestamp, "yyyy-mm-dd HH:MM:SS"))"
        ax.set_title("Station $(station_label) Flux footprints - block $(j)/$(nframes)$(title_time)")
        return lines
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

end # module ffp_block
