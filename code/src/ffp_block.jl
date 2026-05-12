######################################################
###       HELPERS FOR PER-BLOCK FOOTPRINTS         ###
###            author: Michi Haugeneder            ###
######################################################
module ffp_block

using Dates, DataFrames, Statistics, ProgressMeter

import ..kljun
import ..turb

export block_footprint_inputs, calculate_footprints, extract_block_fluxes,
    footprint_dict, full_block_indices, nan_footprint, nanmean, nanstd,
    rows_per_period, valid_footprint_inputs

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

function valid_footprint_inputs(zm, umean, h, ol, sigmav, ustar, wind_dir)
    values = (zm, umean, h, ol, sigmav, ustar, wind_dir)
    any(ismissing, values) && return false
    all(x -> isa(x, Number) && isfinite(x), values) || return false
    zm > 0 || return false
    umean > 0 || return false
    h > 10 || return false
    zm < h || return false
    ol != 0 || return false
    zm / ol > -15.5 || return false
    sigmav > 0 || return false
    ustar > 0.1 || return false
    0 <= wind_dir <= 360 || return false
    return true
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
        valid=fill(false, nblocks),
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

        inputs.valid[j] = valid_footprint_inputs(
            zm,
            inputs.umean[j],
            inputs.h[j],
            inputs.ol[j],
            inputs.sigmav[j],
            inputs.ustar[j],
            inputs.wind_dir[j],
        )
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

        if !inputs.valid[j]
            footprints[j] = nan_footprint(rs, "invalid footprint inputs", input_dict)
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

end # module ffp_per_flux_value
