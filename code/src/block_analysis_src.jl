######################################################
###          HELPERS FOR BLOCK ANALYSIS            ###
###            author: Michi Haugeneder            ###
######################################################

module block_analyze

using DataFrames, Statistics

export add_axes_diagonal!, column_title, finite_series, flux_index,
    lead_flux_correlations, lead_flux_xy, panel_title, plot_block_timeseries_panel!,
    plot_lead_flux_panel!, safe_cor, scatter_xy, tied_ranks, valid_lead_value,
    valid_number

function valid_number(value)
    return !ismissing(value) && value isa Number && isfinite(Float64(value))
end

function flux_index(flux_name::Symbol)
    name = String(flux_name)
    startswith(name, "fx") || error("Expected flux names like :fx1, got :$(name).")
    return parse(Int, name[3:end])
end

function panel_title(flux_name::Symbol, surface_type, heights)
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

function plot_lead_flux_panel!(ax, block_data::AbstractDict, flux_name::Symbol,
        flux_column::Symbol, conversion_factor::Real, ylabel, ylim, surface_type,
        heights; color="C0", show_xlabel=true, show_ylabel=true)
    data = block_data[flux_name]
    lead_fraction, flux = scatter_xy(data, flux_column, conversion_factor)

    #add_axes_diagonal!(ax)
    ax.scatter(lead_fraction, flux, s=12, alpha=0.65, color=color,
        edgecolors="none", zorder=2)
    ax.set_title(panel_title(flux_name, surface_type, heights))
    ax.set_xlabel(show_xlabel ? "Lead fraction" : "")
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

function valid_lead_value(value)
    return valid_number(value) && 0.0 <= Float64(value) <= 1.0
end

function lead_flux_xy(data::AbstractDataFrame, flux_column::Symbol, conversion_factor::Real)
    lead_fraction = data[!, :footprint_lead_fraction]
    flux = data[!, flux_column]
    valid = [
        valid_lead_value(lead_fraction[ix]) && valid_number(flux[ix])
        for ix in eachindex(lead_fraction)
    ]
    return Float64.(lead_fraction[valid]), Float64.(flux[valid]) .* conversion_factor
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

function lead_flux_correlations(correlation_specs, block_data::AbstractDict,
        station_name::AbstractString, surface_type, heights)
    rows = []
    for (flux_name, flux_column, flux_label, conversion_factor) in correlation_specs
        data = block_data[flux_name]
        lead_fraction, flux = lead_flux_xy(data, flux_column, conversion_factor)
        sensor_ix = flux_index(flux_name)

        push!(rows, (
            station=station_name,
            flux_name=String(flux_name),
            flux_variable=String(flux_column),
            flux_label=flux_label,
            surface_type=surface_type[sensor_ix],
            measurement_height_m=heights[sensor_ix],
            n_valid=length(lead_fraction),
            pearson_r=safe_cor(lead_fraction, flux),
            spearman_r=safe_cor(tied_ranks(lead_fraction), tied_ranks(flux)),
        ))
    end

    return DataFrame(rows)
end

function column_title(flux_names, surface_type)
    titles = unique(surface_type[collect(flux_index.(flux_names))])
    return join(titles, " / ")
end

function plot_block_timeseries_panel!(ax, block_data::AbstractDict, flux_names,
        flux_column::Symbol, conversion_factor::Real, ylabel, heights;
        show_xlabel=false, show_ylabel=true, show_right_ylabel=true, show_legend=true)
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

        lead_fraction = finite_series(data[!, :footprint_lead_fraction])
        if any(isfinite, lead_fraction) && any(isfinite, flux)
            ax_right.plot(data.time, lead_fraction, color=color, linestyle="--",
                linewidth=1.0, alpha=0.8)
        end
    end

    ax.set_ylabel(show_ylabel ? ylabel : "")
    ax.set_xlabel(show_xlabel ? "Time" : "")
    ax.grid(alpha=0.9)
    ax_right.set_ylabel(show_right_ylabel ? "Lead fraction" : "")
    ax_right.set_ylim(0, 1.01)

    if show_legend && !isempty(ax.get_legend_handles_labels()[1])
        ax.legend(loc="best")
    end

    return ax, ax_right
end

end # module block_analyze
