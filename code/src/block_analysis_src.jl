######################################################
###          HELPERS FOR BLOCK ANALYSIS            ###
###            author: Michi Haugeneder            ###
######################################################

module block_analyze

using DataFrames, Statistics, Dates, LaTeXStrings
using GLM, StatsModels

export add_axes_diagonal!, column_title, finite_series, flux_index,
    feature_file_label, feature_flux_correlations, feature_flux_difference_correlations,
    feature_flux_difference_fit_summary, feature_flux_difference_xy, feature_flux_xy,
    fit_feature_flux_difference, panel_title,
    plot_block_difference_timeseries_panel!, plot_block_timeseries_panel!,
    plot_feature_flux_difference_fit!, plot_feature_flux_difference_panel!,
    plot_feature_flux_panel!, safe_cor, standard_difference_specs,
    read_subplot_order, scatter_xy, sensor_difference_label, tied_ranks, valid_feature_value,
    valid_number, normalize_colored_scatter_variable, colored_scatter_prerequisite_note,
    colored_scatter_time_origin, colored_flux_time_origin, colored_difference_data,
    colored_correlation_data, colored_scatter_colorbar_label

const FEATURE_FRACTION_COLUMN = :footprint_feature_fraction
const LEGACY_FRACTION_COLUMN = :footprint_lead_fraction
const CLASSIFIED_WEIGHT_FRACTION_COLUMN = :footprint_classified_weight_fraction

function valid_number(value)
    return !ismissing(value) && value isa Number && isfinite(Float64(value))
end

function flux_index(flux_name::Symbol)
    name = String(flux_name)
    startswith(name, "fx") || error("Expected flux names like :fx1, got :$(name).")
    return parse(Int, name[3:end])
end

function feature_file_label(feature_label::AbstractString)
    label = replace(lowercase(strip(feature_label)), r"[^a-z0-9]+" => "_")
    label = replace(label, r"^_+|_+$" => "")
    return isempty(label) ? "feature" : label
end

function _flux_symbol(value, config_key::AbstractString)
    if value isa Symbol
        return value
    elseif value isa AbstractString
        name = startswith(value, ":") ? value[2:end] : value
        isempty(name) && error("Config key block_analysis.$config_key contains an empty flux name.")
        return Symbol(name)
    end

    error("Config key block_analysis.$config_key must contain flux names as strings.")
end

function read_subplot_order(config::AbstractDict, key, default_order;
        expected_length::Integer=length(default_order))
    key_string = String(key)
    block_config = get(config, "block_analysis", nothing)

    raw_order = if isnothing(block_config)
        default_order
    elseif block_config isa AbstractDict
        get(block_config, key_string, default_order)
    else
        error("Config key block_analysis must be a table.")
    end

    raw_order isa AbstractVector || error(
        "Config key block_analysis.$key_string must be a list of flux names."
    )

    order = [_flux_symbol(value, key_string) for value in raw_order]
    length(order) == expected_length || error(
        "Config key block_analysis.$key_string must contain $expected_length flux names; got $(length(order))."
    )
    foreach(flux_index, order)
    return order
end

function panel_title(flux_name::Symbol, surface_type, heights)
    ix = flux_index(flux_name)
    return "$(surface_type[ix]), $(heights[ix]) m"
end

function feature_fraction_column(data::AbstractDataFrame)
    if FEATURE_FRACTION_COLUMN in propertynames(data)
        return FEATURE_FRACTION_COLUMN
    elseif LEGACY_FRACTION_COLUMN in propertynames(data)
        return LEGACY_FRACTION_COLUMN
    end

    error("Block data must contain :$(FEATURE_FRACTION_COLUMN).")
end

function feature_fraction_series(data::AbstractDataFrame)
    feature_column = feature_fraction_column(data)
    feature_fraction = data[!, feature_column]

    if feature_column == FEATURE_FRACTION_COLUMN &&
            CLASSIFIED_WEIGHT_FRACTION_COLUMN in propertynames(data)
        # Stored feature fractions are conditional on classified pixels; scale them back
        # to the whole footprint so outside/missing surface stays unknown downstream.
        classified_fraction = data[!, CLASSIFIED_WEIGHT_FRACTION_COLUMN]
        return [
            valid_number(feature_fraction[ix]) && valid_number(classified_fraction[ix]) ?
                Float64(feature_fraction[ix]) * Float64(classified_fraction[ix]) : missing
            for ix in eachindex(feature_fraction)
        ]
    end

    return feature_fraction
end

function scatter_xy(data::AbstractDataFrame, flux_column::Symbol, conversion_factor::Real)
    feature_fraction = feature_fraction_series(data)
    flux = data[!, flux_column]
    valid = [
        valid_number(feature_fraction[ix]) && valid_number(flux[ix])
        #&& 0.0 <= Float64(feature_fraction[ix]) <= 1.0
        for ix in eachindex(feature_fraction)
    ]
    return Float64.(feature_fraction[valid]), Float64.(flux[valid]) .* conversion_factor
end

function add_axes_diagonal!(ax)
    # Feature fraction and flux have different units; draw this as an axes-space guide.
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="grey",
        alpha=0.35, linewidth=1.2, zorder=1)
end

function plot_feature_flux_panel!(ax, block_data::AbstractDict, flux_name::Symbol,
        flux_column::Symbol, conversion_factor::Real, ylabel, ylim, surface_type,
        heights; feature_label::AbstractString="feature", color="C0",
        show_xlabel=true, show_ylabel=true)
    data = block_data[flux_name]
    feature_fraction, flux = scatter_xy(data, flux_column, conversion_factor)

    #add_axes_diagonal!(ax)
    ax.scatter(feature_fraction, flux, s=12, alpha=0.65, color=color,
        edgecolors="none", zorder=2)
    ax.set_title(panel_title(flux_name, surface_type, heights))
    ax.set_xlabel(show_xlabel ? "$(feature_label) fraction" : "")
    ax.set_ylabel(show_ylabel ? ylabel : "")
    ax.set_xlim(0, 1)
    ax.set_ylim(ylim)
    ax.set_axisbelow(true)
    ax.grid(alpha=0.9)

    return ax
end

function finite_series(values)
    return [valid_number(value) ? Float64(value) : NaN for value in values]
end

function valid_feature_value(value)
    return valid_number(value) && 0.0 <= Float64(value) <= 1.0
end

function feature_flux_xy(data::AbstractDataFrame, flux_column::Symbol, conversion_factor::Real)
    feature_fraction = feature_fraction_series(data)
    flux = data[!, flux_column]
    valid = [
        valid_feature_value(feature_fraction[ix]) && valid_number(flux[ix])
        for ix in eachindex(feature_fraction)
    ]
    return Float64.(feature_fraction[valid]), Float64.(flux[valid]) .* conversion_factor
end

function tied_ranks(values::AbstractVector{<:Real})
    order = sortperm(values)
    ranks = Vector{Float64}(undef, length(values))
    ix = 1

    while ix <= length(values)
        jx = ix
        while jx < length(values) && values[order[jx + 1]] == values[order[ix]]
            jx += 1
        end

        rank = (ix + jx) / 2
        for kx in ix:jx
            ranks[order[kx]] = rank
        end
        ix = jx + 1
    end

    return ranks
end

function safe_cor(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) >= 2 || return NaN
    std(x) > 0 && std(y) > 0 || return NaN
    return cor(x, y)
end

function feature_flux_correlations(correlation_specs, block_data::AbstractDict,
        station_name::AbstractString, surface_type, heights; feature_label::AbstractString="feature")
    rows = []
    for (flux_name, flux_column, flux_label, conversion_factor) in correlation_specs
        data = block_data[flux_name]
        feature_fraction, flux = feature_flux_xy(data, flux_column, conversion_factor)
        sensor_ix = flux_index(flux_name)

        push!(rows, (
            station=station_name,
            feature=feature_label,
            flux_name=String(flux_name),
            flux_variable=String(flux_column),
            flux_label=flux_label,
            surface_type=surface_type[sensor_ix],
            measurement_height_m=heights[sensor_ix],
            n_valid=length(feature_fraction),
            pearson_r=safe_cor(feature_fraction, flux),
            spearman_r=safe_cor(tied_ranks(feature_fraction), tied_ranks(flux)),
        ))
    end

    return DataFrame(rows)
end

function column_title(flux_names, surface_type)
    titles = unique(surface_type[collect(flux_index.(flux_names))])
    return join(titles, " / ")
end

function sensor_difference_label(flux_names, heights, surface_type)
    length(flux_names) == 2 || error("Expected exactly two flux names.")
    first_name, second_name = flux_names
    first_label = "$(String(surface_type[flux_index(first_name)])) $(heights[flux_index(first_name)]) m"
    second_label = "$(String(surface_type[flux_index(second_name)])) $(heights[flux_index(second_name)]) m"
    return "$(first_label) - $(second_label)"
end

function standard_difference_specs(subplot_order, latent_subplot_order,
        sensible_conversion_factor::Real, latent_conversion_factor::Real)
    length(subplot_order) >= 4 || error("subplot_order must contain at least four flux names.")
    length(latent_subplot_order) >= 2 || error("latent_subplot_order must contain at least two flux names.")

    return [
        ((subplot_order[1], subplot_order[2]), :wT, sensible_conversion_factor, "sensible",
            L"\Delta H~\mathrm{[W~m^{-2}]}"),
        ((subplot_order[3], subplot_order[4]), :wT, sensible_conversion_factor, "sensible",
            L"\Delta H~\mathrm{[W~m^{-2}]}"),
        ((latent_subplot_order[1], latent_subplot_order[2]), :wq, latent_conversion_factor, "latent",
            L"\Delta L_E~\mathrm{[W~m^{-2}]}"),
    ]
end

function difference_timeseries(block_data::AbstractDict, flux_names,
        flux_column::Symbol, conversion_factor::Real)
    length(flux_names) == 2 || error("Expected exactly two flux names.")
    first_name, second_name = flux_names
    first_data = block_data[first_name]
    second_data = block_data[second_name]

    first_df = DataFrame(
        time=first_data[!, :time],
        flux_first=first_data[!, flux_column],
        feature_first=feature_fraction_series(first_data),
    )
    second_df = DataFrame(
        time=second_data[!, :time],
        flux_second=second_data[!, flux_column],
        feature_second=feature_fraction_series(second_data),
    )
    joined = innerjoin(first_df, second_df, on=:time)

    flux_difference = (finite_series(joined[!, :flux_first]) .-
                       finite_series(joined[!, :flux_second])) .* conversion_factor
    feature_difference = finite_series(joined[!, :feature_first]) .-
                         finite_series(joined[!, :feature_second])

    return joined.time, flux_difference, feature_difference
end

function feature_flux_difference_xy(block_data::AbstractDict, flux_names,
        flux_column::Symbol, conversion_factor::Real)
    _, flux_difference, feature_difference = difference_timeseries(
        block_data, flux_names, flux_column, conversion_factor)
    valid = isfinite.(feature_difference) .& isfinite.(flux_difference)
    return feature_difference[valid], flux_difference[valid]
end

function fit_feature_flux_difference(feature_difference::AbstractVector,
        flux_difference::AbstractVector; xgrid=nothing, grid_length::Integer=100)
    length(feature_difference) == length(flux_difference) || error(
        "feature_difference and flux_difference must have the same length."
    )

    valid = isfinite.(feature_difference) .& isfinite.(flux_difference)
    x = Float64.(feature_difference[valid])
    y = Float64.(flux_difference[valid])

    length(x) >= 2 || error("At least two valid points are required for a linear fit.")
    length(unique(x)) >= 2 || error("At least two unique feature-difference values are required for a linear fit.")

    fit_df = DataFrame(feature_difference=x, flux_difference=y)
    model = lm(@formula(flux_difference ~ feature_difference), fit_df)
    fit_vars = coef(model)
    fit_vars_confint = confint(model, level=0.95)
    fit_stderror = stderror(model)

    x_grid = isnothing(xgrid) ?
        collect(range(minimum(x), maximum(x), length=grid_length)) :
        Float64.(collect(xgrid))
    preds_confidence = predict(model, DataFrame(feature_difference=x_grid), interval=:confidence)

    return (
        model=model,
        feature_difference=x,
        flux_difference=y,
        n_valid=length(x),
        coef=fit_vars,
        confint=fit_vars_confint,
        stderror=fit_stderror,
        r2=r2(model),
        loglikelihood=loglikelihood(model),
        x_grid=x_grid,
        predictions_confidence=preds_confidence,
    )
end

function feature_flux_difference_fit_summary(fit, feature_label::AbstractString,
        flux_column::Symbol, flux_kind::AbstractString, difference_label::AbstractString)
    return (
        feature=feature_label,
        flux_variable=String(flux_column),
        flux_label=String(flux_kind),
        difference_label=difference_label,
        n_valid=fit.n_valid,
        intercept=fit.coef[1],
        intercept_ci_low=fit.confint[1, 1],
        intercept_ci_high=fit.confint[1, 2],
        intercept_stderr=fit.stderror[1],
        slope=fit.coef[2],
        slope_low=fit.confint[2, 1],
        slope_high=fit.confint[2, 2],
        slope_stderr=fit.stderror[2],
        r2=fit.r2,
        loglikelihood=fit.loglikelihood,
    )
end

function plot_feature_flux_difference_fit!(ax, fit; color="black",
        label::AbstractString="linear fit", ci_label::AbstractString="95% CI",
        ci_alpha::Real=0.18, linewidth::Real=1.0)
    ax.fill_between(fit.x_grid,
        fit.predictions_confidence.lower,
        fit.predictions_confidence.upper;
        color=color, alpha=ci_alpha, linewidth=0, label=ci_label)
    ax.plot(fit.x_grid, fit.predictions_confidence.prediction;
        color=color, linewidth=linewidth, label=label)

    return ax
end

function feature_flux_difference_correlations(difference_specs, block_data::AbstractDict,
        station_name::AbstractString, surface_type, heights; feature_label::AbstractString="feature")
    rows = []
    for (panel_ix, (flux_names, flux_column, conversion_factor, flux_kind, _)) in
            enumerate(difference_specs)
        _, flux_difference, feature_difference = difference_timeseries(
            block_data, flux_names, flux_column, conversion_factor)
        valid = isfinite.(feature_difference) .& isfinite.(flux_difference)
        feature_valid = feature_difference[valid]
        flux_valid = flux_difference[valid]
        first_name, second_name = flux_names

        push!(rows, (
            feature=feature_label,
            flux_variable=String(flux_column),
            flux_label=String(flux_kind),
            difference_label=sensor_difference_label(flux_names, heights, surface_type),
            first_surface_type=surface_type[flux_index(first_name)],
            second_surface_type=surface_type[flux_index(second_name)],
            first_measurement_height_m=heights[flux_index(first_name)],
            second_measurement_height_m=heights[flux_index(second_name)],
            n_valid=length(feature_valid),
            pearson_r=safe_cor(feature_valid, flux_valid),
            spearman_r=safe_cor(tied_ranks(feature_valid), tied_ranks(flux_valid)),
        ))
    end

    return DataFrame(rows)
end

function symmetric_limits(values; fallback=(-1.0, 1.0), padding=0.05)
    finite_values = values[isfinite.(values)]
    isempty(finite_values) && return fallback

    max_abs = maximum(abs.(finite_values))
    max_abs > 0 || return fallback
    limit = max_abs * (1 + padding)
    return (-limit, limit)
end

function plot_block_difference_timeseries_panel!(ax, block_data::AbstractDict,
        flux_names, flux_column::Symbol, conversion_factor::Real, flux_kind::AbstractString,
        flux_ylabel, heights, surface_type, wT_limits, wq_limits;
        feature_label::AbstractString="feature", show_xlabel=true, show_ylabel=true,
        show_right_ylabel=true)
    ax_right = ax.twinx()
    time, flux_difference, feature_difference = difference_timeseries(
        block_data, flux_names, flux_column, conversion_factor)
    difference_label = sensor_difference_label(flux_names, heights, surface_type)

    if any(isfinite, flux_difference)
        flux_line = ax.plot(time, flux_difference, color="black", linewidth=1.1, alpha=0.8,
            label="Flux difference")
    end
    if any(isfinite, feature_difference)
        fraction_line = ax_right.plot(time, feature_difference, color="C3", linestyle="--", linewidth=1.0,
            alpha=0.8, label="$(feature_label) fraction difference")
    end

    ax.axhline(0, color="grey", linewidth=0.8, alpha=0.55)
    ax.set_title("$(flux_kind): $(difference_label)")
    ax.set_xlabel(show_xlabel ? "Time" : "")
    ax.set_ylabel(show_ylabel ? "$(flux_ylabel)" : "")
    ax.set_ylim(symmetric_limits(flux_difference))
    ax.grid(alpha=0.9)
    ax_right.set_ylabel(show_right_ylabel ? "$(feature_label) fraction difference" : "")
    flux_y_label = flux_column==:wT ? wT_limits : wq_limits
    ax.set_ylim(flux_y_label)
    #ax_right.set_ylim(-1.01, 1.01)
    ax_right.legend()

    return ax, ax_right
end

function plot_feature_flux_difference_panel!(ax, block_data::AbstractDict,
        flux_names, flux_column::Symbol, conversion_factor::Real, flux_kind::AbstractString,
        flux_ylabel, ylim, surface_type, heights; feature_label::AbstractString="feature",
        color="black", show_xlabel=true, show_ylabel=true)
    _, flux_difference, feature_difference = difference_timeseries(
        block_data, flux_names, flux_column, conversion_factor)
    valid = isfinite.(feature_difference) .& isfinite.(flux_difference)
    difference_label = sensor_difference_label(flux_names, heights, surface_type)

    ax.scatter(feature_difference[valid], flux_difference[valid], s=12, alpha=0.65,
        color=color, edgecolors="none", zorder=2)
    ax.axhline(0, color="grey", linewidth=0.8, alpha=0.55)
    ax.axvline(0, color="grey", linewidth=0.8, alpha=0.55)
    ax.set_title("$(flux_kind): $(difference_label)")
    ax.set_xlabel(show_xlabel ? "$(feature_label) fraction difference" : "")
    ax.set_ylabel(show_ylabel ? flux_ylabel : "")
    #ax.set_xlim()
    ax.set_ylim(ylim)
    ax.set_axisbelow(true)
    ax.grid(alpha=0.9)

    return ax
end

function plot_block_timeseries_panel!(ax, block_data::AbstractDict, flux_names,
        flux_column::Symbol, conversion_factor::Real, ylabel, heights;
        feature_label::AbstractString="feature", show_xlabel=false, show_ylabel=true,
        show_right_ylabel=true, show_legend=true)
    ax_right = ax.twinx()
    colors = ["C0", "C1"]

    for (ix, flux_name) in enumerate(flux_names)
        data = block_data[flux_name]
        color = colors[ix]

        flux = finite_series(data[!, flux_column]) .* conversion_factor
        if any(isfinite, flux)
            label = "$(heights[flux_index(flux_name)]) m"
            ax.plot(data.time, flux, color=color, linewidth=1.1, label=label)
        end

        feature_fraction = finite_series(feature_fraction_series(data))
        if any(isfinite, feature_fraction) && any(isfinite, flux)
            ax_right.plot(data.time, feature_fraction, color=color, linestyle="--",
                linewidth=1.0, alpha=0.8)
        end
    end

    ax.set_ylabel(show_ylabel ? ylabel : "")
    ax.set_xlabel(show_xlabel ? "Time" : "")
    ax.grid(alpha=0.9)
    ax_right.set_ylabel(show_right_ylabel ? "$(feature_label) fraction" : "")
    ax_right.set_ylim(0, 1.01)

    if show_legend && !isempty(ax.get_legend_handles_labels()[1])
        ax.legend(loc="best")
    end

    return ax, ax_right
end

function normalize_colored_scatter_variable(value)
    raw = lowercase(strip(String(value)))
    if raw == "u"
        return "u"
    elseif raw in ("t", "temperature", "sonic_temperature", "sonic_temp")
        return "T"
    elseif raw == "wt"
        return "wT"
    elseif raw == "wq"
        return "wq"
    elseif raw == "time"
        return "time"
    elseif raw in ("wind_dir", "winddir", "wind_direction", "wd")
        return "wind_dir"
    end

    error("colored_scatter_variable must be one of \"u\", \"T\", \"wind_dir\", \"wT\", \"wq\", or \"time\".")
end

function colored_scatter_prerequisite_note(color_variable::AbstractString)
    if color_variable == "u"
        return "Source variable: load_data.jl. Rerun load_data.jl, turb_fluxes.jl, and ffp_per_flux_value.jl so :u is written to the block NetCDF."
    elseif color_variable == "T"
        return "Source variable: load_data.jl. Rerun load_data.jl, turb_fluxes.jl, and ffp_per_flux_value.jl so :T is written to the block NetCDF."
    elseif color_variable == "wind_dir"
        return "Source variable: load_data.jl. Rerun load_data.jl, turb_fluxes.jl, and ffp_per_flux_value.jl so :wind_dir is written to the block NetCDF."
    elseif color_variable == "wT" || color_variable == "wq"
        return "Source variable: turb_fluxes.jl. Rerun load_data.jl, turb_fluxes.jl, and ffp_per_flux_value.jl so fluxes are written to the block NetCDF."
    elseif color_variable == "time"
        return "Time should be written by ffp_per_flux_value.jl; rerun it if the block NetCDF is stale."
    end

    return ""
end

function colored_scatter_color_column(color_variable::AbstractString)
    color_variable == "u" && return :u
    color_variable == "T" && return :T
    color_variable == "wT" && return :wT
    color_variable == "wq" && return :wq
    color_variable == "wind_dir" && return :wind_dir
    error("No data column for color variable \"$(color_variable)\".")
end

function colored_scatter_conversion_factor(color_variable::AbstractString)
    ρ_air = 1.2 #kg m^{-3}
    c_p = 1004 #J kg^{-1} K^{-1}
    L_v = 2500e3 #J kg^{-1} (approx @0°C)
    color_variable == "wT" && return ρ_air * c_p
    color_variable == "wq" && return L_v * 1e-3
    return 1.0
end

function colored_scatter_colorbar_label(color_variable::AbstractString, time_origin)
    if color_variable == "u"
        return L"\overline{u}~\mathrm{[m~s^{-1}]}"
    elseif color_variable == "T"
        return L"\overline{T_{sonic}}~\mathrm{[^\circ C]}"
    elseif color_variable == "wind_dir"
        return L"\overline{\alpha}~\mathrm{[^\circ]}"
    elseif color_variable == "wT"
        return L"\overline{w'T'}~\mathrm{[W~m^{-2}]}"
    elseif color_variable == "wq"
        return L"\overline{w'q'}~\mathrm{[W~m^{-2}]}"
    elseif color_variable == "time" && !isnothing(time_origin)
        return "time since $(Dates.format(time_origin, "yyyy-mm-dd HH:MM")) [h]"
    end

    return color_variable
end

function colored_scatter_time_origin(block_data::AbstractDict, difference_specs)
    times = DateTime[]
    for (flux_names, _, _, _, _) in difference_specs
        for flux_name in flux_names
            data = block_data[flux_name]
            :time in propertynames(data) || continue
            append!(times, [time for time in data.time if time isa DateTime])
        end
    end

    return isempty(times) ? nothing : minimum(times)
end

function colored_flux_time_origin(block_data::AbstractDict, flux_names)
    times = DateTime[]
    for flux_name in flux_names
        data = block_data[flux_name]
        :time in propertynames(data) || continue
        append!(times, [time for time in data.time if time isa DateTime])
    end

    return isempty(times) ? nothing : minimum(times)
end

function colored_scatter_elapsed_hours(value, time_origin)
    value isa DateTime && time_origin isa DateTime || return NaN
    return Dates.value(value - time_origin) / (1000 * 60 * 60)
end

function mean_wind_direction_pair(first_value, second_value)
    valid_number(first_value) && valid_number(second_value) || return NaN

    angles_rad = deg2rad.([Float64(first_value), Float64(second_value)])
    east_component = mean(sin.(angles_rad))
    north_component = mean(cos.(angles_rad))
    return mod(rad2deg(atan(east_component, north_component)), 360)
end

function colored_scatter_pair_value(first_value, second_value, color_variable::AbstractString)
    block_analyze.valid_number(first_value) && block_analyze.valid_number(second_value) || return NaN

    first_number = Float64(first_value)
    second_number = Float64(second_value)
    if color_variable == "wind_dir"
        return mean_wind_direction_pair(first_number, second_number)
    end
    if color_variable == "wT" || color_variable == "wq"
        return 0.5 * (first_number + second_number) *
            colored_scatter_conversion_factor(color_variable)
    end

    return 0.5 * (first_number + second_number)
end

function empty_colored_difference_data(message::AbstractString)
    return (
        available=false,
        message=String(message),
        time=DateTime[],
        feature_difference=Float64[],
        flux_difference=Float64[],
        color_values=Float64[],
    )
end

function empty_colored_correlation_data(message::AbstractString)
    return (
        available=false,
        message=String(message),
        time=DateTime[],
        feature_fraction=Float64[],
        flux=Float64[],
        color_values=Float64[],
    )
end

function colored_correlation_data(block_data::AbstractDict, flux_name::Symbol,
        flux_column::Symbol, conversion_factor::Real, color_variable::AbstractString,
        time_origin)
    data = block_data[flux_name]

    if color_variable == "time" && isnothing(time_origin)
        return empty_colored_correlation_data(
            "Color variable \"time\" is unavailable. $(colored_scatter_prerequisite_note(color_variable))"
        )
    end

    feature_fraction = finite_series(feature_fraction_series(data))
    flux = finite_series(data[!, flux_column]) .* conversion_factor

    color_values = if color_variable == "time"
        [colored_scatter_elapsed_hours(value, time_origin) for value in data.time]
    else
        color_column = colored_scatter_color_column(color_variable)
        if !(color_column in propertynames(data))
            return empty_colored_correlation_data(
                "Color variable \"$(color_variable)\" is unavailable ($(String(flux_name)).$(String(color_column)) missing). " *
                colored_scatter_prerequisite_note(color_variable)
            )
        end

        [
            valid_number(value) ?
                Float64(value) * colored_scatter_conversion_factor(color_variable) : NaN
            for value in data[!, color_column]
        ]
    end

    valid = isfinite.(feature_fraction) .& isfinite.(flux) .& isfinite.(color_values)
    return (
        available=true,
        message="",
        time=data.time[valid],
        feature_fraction=feature_fraction[valid],
        flux=flux[valid],
        color_values=color_values[valid],
    )
end

function colored_difference_data(block_data::AbstractDict, flux_names,
        flux_column::Symbol, conversion_factor::Real, color_variable::AbstractString,
        time_origin)
    length(flux_names) == 2 || error("Expected exactly two flux names.")
    first_name, second_name = flux_names
    first_data = block_data[first_name]
    second_data = block_data[second_name]

    if color_variable == "time" && isnothing(time_origin)
        return empty_colored_difference_data(
            "Color variable \"time\" is unavailable. $(colored_scatter_prerequisite_note(color_variable))"
        )
    end

    first_df = DataFrame(
        time=first_data[!, :time],
        flux_first=first_data[!, flux_column],
        feature_first=block_analyze.feature_fraction_series(first_data),
    )
    second_df = DataFrame(
        time=second_data[!, :time],
        flux_second=second_data[!, flux_column],
        feature_second=block_analyze.feature_fraction_series(second_data),
    )

    if color_variable != "time"
        color_column = colored_scatter_color_column(color_variable)
        missing_columns = String[]
        color_column in propertynames(first_data) || push!(missing_columns,
            "$(String(first_name)).$(String(color_column))")
        color_column in propertynames(second_data) || push!(missing_columns,
            "$(String(second_name)).$(String(color_column))")

        if !isempty(missing_columns)
            return empty_colored_difference_data(
                "Color variable \"$(color_variable)\" is unavailable ($(join(missing_columns, ", ")) missing). " *
                colored_scatter_prerequisite_note(color_variable)
            )
        end

        first_df[!, :color_first] = first_data[!, color_column]
        second_df[!, :color_second] = second_data[!, color_column]
    end

    joined = innerjoin(first_df, second_df, on=:time)
    flux_difference = (block_analyze.finite_series(joined[!, :flux_first]) .-
                       block_analyze.finite_series(joined[!, :flux_second])) .* conversion_factor
    feature_difference = block_analyze.finite_series(joined[!, :feature_first]) .-
                         block_analyze.finite_series(joined[!, :feature_second])

    color_values = if color_variable == "time"
        [colored_scatter_elapsed_hours(value, time_origin) for value in joined.time]
    else
        [colored_scatter_pair_value(joined.color_first[ix], joined.color_second[ix], color_variable)
            for ix in 1:nrow(joined)]
    end

    valid = isfinite.(feature_difference) .& isfinite.(flux_difference) .& isfinite.(color_values)
    return (
        available=true,
        message="",
        time=joined.time[valid],
        feature_difference=feature_difference[valid],
        flux_difference=flux_difference[valid],
        color_values=color_values[valid],
    )
end

end # module block_analyze
