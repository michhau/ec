######################################################
###        HISTOGRAMS OF TURBULENT FLUXES          ###
###            author: Michi Haugeneder            ###
######################################################
#=
REPL workflow for flux histograms.

Run the sections below interactively, for example from VS Code with the `##`
cell markers. The expensive path is split into loading, flux calculation,
averaging, saving, grouping, and plotting. To avoid reading the high-frequency
NetCDF files again, run the settings section and then the cache-reading section.
=#

using Dates, DataFrames, Statistics, LaTeXStrings, ProgressMeter
import CSV, PyPlot

importdir = joinpath(@__DIR__, "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "station_config.jl"))
import .turb
import .gen
import .stationcfg

const RHO_AIR = 1.2       # kg m^-3
const C_P = 1004.0        # J kg^-1 K^-1
const L_V = 2500e3        # J kg^-1

const FLUX_SURFACE_TYPES = ("ice", "lead", "pond", "ridge")
const DEFAULT_SENSIBLE_BINS = collect(-50.0:1.0:50.0)
const DEFAULT_LATENT_BINS = collect(-25.0:0.5:25.0)

"""
    read_dship_meteo_csv(filename::String, datatypes::Vector{DataType}, dateformat_inp=DateFormat("yyyymmddTHHMMSS"))

Read a CSV file extracted from DSHIP.

# Arguments
- `filename::String`: Path to the CSV file

# Returns
- DataFrame with the data
"""
function read_dship_meteo_csv(filename::String, datatypes::Vector{DataType}, dateformat_inp=DateFormat("yyyymmddTHHMMSS"))
    a = CSV.read(filename, DataFrame, header=1, skipto=4,  dateformat=dateformat_inp, types=datatypes, maxwarnings=10)

    #extract and rename colnames (without the generic "Polarstern" leading text)
    colnames_raw = names(a)
    colnames_final = Vector{String}(undef, length(colnames_raw))
    colnames_final[1] = "time"
    reg_test = r"\.([^.]*)$"
    for i in 2:length(colnames_raw)
        colnames_final[i] = first(Tuple(match(reg_test, colnames_raw[i])))
    end
    rename!(a, colnames_final)
        return a
end

period_seconds(period::Period) = period / Second(1)

function selected_station_names(stations=nothing)
    if isnothing(stations)
        return stationcfg.available_station_names()
    elseif stations isa AbstractString
        return [strip(name) for name in split(stations, ",") if !isempty(strip(name))]
    end

    return String.(stations)
end

load_station_configs(station_names) = [stationcfg.load_station_config(name) for name in station_names]

function default_cache_file(station_configs, avg_period::Period)
    first_config = first(station_configs)
    data_root = normpath(String(stationcfg.require_key(first_config, "data_root")))
    cache_root = basename(data_root) == "cut" ? dirname(data_root) : data_root
    return joinpath(cache_root, "for_hist", "fluxes_histogram_$(round(Int, period_seconds(avg_period)))s.csv")
end

default_plot_dir(station_configs) = stationcfg.plot_dir(first(station_configs), "wT_wq")

function normalize_flux_surface_type(value::AbstractString)
    value in FLUX_SURFACE_TYPES && return value
    error("Unknown flux histogram surface type '$value'. Expected one of $(join(FLUX_SURFACE_TYPES, ", ")).")
end

function flux_surface_types(station_config)
    configured = stationcfg.optional_key(station_config, nothing, "flux_histogram", "surface_type")
    values = isnothing(configured) ? stationcfg.require_key(station_config, "surface_type") : configured
    types = normalize_flux_surface_type.(String.(values))
    length(types) == 4 || error("flux_histogram.surface_type must define exactly four sensor categories.")
    return types
end

function sensor_height_bands(station_config)
    heights = Float64.(stationcfg.require_key(station_config, "heights"))
    length(heights) == 4 || error("heights must define exactly four sensor heights.")

    bands = Vector{String}(undef, 4)
    for (lower_candidate, upper_candidate) in ((1, 2), (3, 4))
        if heights[lower_candidate] <= heights[upper_candidate]
            bands[lower_candidate] = "1m"
            bands[upper_candidate] = "2m"
        else
            bands[lower_candidate] = "2m"
            bands[upper_candidate] = "1m"
        end
    end
    return bands
end

function station_nan_files(station_config)
    nan_root = String(stationcfg.optional_key(
        station_config,
        "/home/haugened/Documents/data/CONTRASTS/nan_periods",
        "manual_nanmask",
        "root",
    ))
    nan_files = String.(stationcfg.optional_key(
        station_config,
        ["t1irg_nan.csv", "t1csat_nan.csv", "t2irg_nan.csv", "t2csat_nan.csv"],
        "manual_nanmask",
        "files",
    ))
    return joinpath.(nan_root, nan_files)
end

function apply_manual_nanmask!(evaldf::DataFrame, station_config, sensor_index::Integer)
    stationcfg.optional_key(station_config, true, "manual_nanmask", "enabled") || return evaldf

    nanfiles = station_nan_files(station_config)
    sensor_index <= length(nanfiles) || return evaldf

    nanmask_columns = String.(stationcfg.optional_key(
        station_config,
        ["u", "v", "w", "T"],
        "manual_nanmask",
        "columns",
    ))
    manual_period = turb.read_nantimes_csv(nanfiles[sensor_index])
    isbad_quality = turb.manual_nanmask(manual_period, evaldf.time)
    evaldf[isbad_quality, nanmask_columns] .= NaN
    return evaldf
end

function load_sensor_data(station_config, data_file::AbstractString, sensor_index::Integer)
    datapath = String(stationcfg.require_key(station_config, "data_root"))
    evalstart = stationcfg.station_datetime(station_config, "evalstart")
    evalend = stationcfg.station_datetime(station_config, "evalend")

    evaldf = turb.readturbasnetcdf(joinpath(datapath, data_file), evalstart, evalend)
    turb.printmissstats(evaldf)
    evaldf = turb.interpolatemissing(evaldf)
    turb.drdf!(evaldf, periodwise=false)
    apply_manual_nanmask!(evaldf, station_config, sensor_index)
    turb.missing2nan!(evaldf)
    return evaldf
end

function load_station_data(station_config)
    data_files = String.(stationcfg.require_key(station_config, "data_files"))
    length(data_files) == 4 || error("Station $(stationcfg.require_key(station_config, "id")) must define four data_files.")

    println()
    println("Loading station ", stationcfg.station_label(station_config), " from ", station_config["config_file"])
    return [load_sensor_data(station_config, data_file, sensor_index)
            for (sensor_index, data_file) in enumerate(data_files)]
end

calculate_fluxes(evaldfs, avg_period::Period) = [turb.turbflux(evaldf, avg_period, 1.0) for evaldf in evaldfs]
average_fluxes(fx_raws, avg_period::Period) = [turb.avgflux(fx_raw, avg_period, true, 0.1) for fx_raw in fx_raws]

function low_frequency_flux(fx::DataFrame, avg_period::Period)
    low = DataFrame()
    time_out, _ = gen.block_average(fx.time, fx.wT, avg_period)
    low[!, :time] = time_out

    for column_name in names(fx)[2:end]
        _, values = gen.block_average(fx.time, fx[!, column_name], avg_period)
        low[!, Symbol(column_name)] = values
    end

    return low
end

low_frequency_fluxes(fx_avgs, avg_period::Period) = [low_frequency_flux(fx_avg, avg_period) for fx_avg in fx_avgs]

function station_flux_records(station_config, low_fluxes)
    length(low_fluxes) == 4 || error("Expected four low-frequency flux DataFrames.")

    station_id = String(stationcfg.require_key(station_config, "id"))
    station_label = stationcfg.station_label(station_config)
    surface_types = String.(stationcfg.require_key(station_config, "surface_type"))
    histogram_surface_types = flux_surface_types(station_config)
    heights = Float64.(stationcfg.require_key(station_config, "heights"))
    height_bands = sensor_height_bands(station_config)
    instrument_types = String.(stationcfg.optional_key(
        station_config,
        ["IRG", "CSAT", "IRG", "CSAT"],
        "instrument_type",
    ))
    instrument_labels = stationcfg.station_labels(station_config)

    station_records = DataFrame[]
    for (sensor_index, fx_low) in enumerate(low_fluxes)
        n = nrow(fx_low)
        records = DataFrame(
            station = fill(station_id, n),
            station_label = fill(station_label, n),
            sensor = fill(sensor_index, n),
            instrument_type = fill(instrument_types[sensor_index], n),
            instrument_label = fill(instrument_labels[sensor_index], n),
            surface_type = fill(surface_types[sensor_index], n),
            flux_surface_type = fill(histogram_surface_types[sensor_index], n),
            height = fill(heights[sensor_index], n),
            height_band = fill(height_bands[sensor_index], n),
            time = fx_low.time,
        )

        for column_name in names(fx_low)[2:end]
            records[!, Symbol(column_name)] = fx_low[!, column_name]
        end

        records[!, :H] = records.wT .* (RHO_AIR * C_P)
        records[!, :LE] = records.wq .* (L_V * 1e-3)
        push!(station_records, records)
    end

    return vcat(station_records...; cols=:union)
end

function parse_datetime_value(value)
    value isa DateTime && return value
    value isa Date && return DateTime(value)
    return DateTime(replace(String(value), " " => "T"))
end

function read_flux_cache(cache_file::AbstractString)
    data = CSV.read(cache_file, DataFrame)
    if :time in propertynames(data) && !(eltype(data.time) <: DateTime)
        data[!, :time] = parse_datetime_value.(data.time)
    end
    return data
end

function write_flux_cache(cache_file::AbstractString, data::DataFrame)
    mkpath(dirname(cache_file))
    CSV.write(cache_file, data)
    return cache_file
end

function finite_values(values)
    out = Float64[]
    for value in values
        if !ismissing(value) && isfinite(Float64(value))
            push!(out, Float64(value))
        end
    end
    return out
end

data_hours(values, avg_period::Period) = length(finite_values(values)) * period_seconds(avg_period) / 3600

function add_hours_text!(ax, text; x=0.97, y=0.95)
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox=Dict("facecolor" => "white", "alpha" => 0.75, "edgecolor" => "none", "pad" => 2),
    )
end

function safe_hist!(ax, values; bins, color, alpha=0.75, label=nothing)
    isempty(values) && return nothing

    if isnothing(label)
        ax.hist(values, bins=bins, density=true, color=color, alpha=alpha)
    else
        ax.hist(values, bins=bins, density=true, color=color, alpha=alpha, label=label)
    end
end

function plot_total_histogram(data::DataFrame, output_file::AbstractString, avg_period::Period;
                              title_suffix="")
    fig, axes = PyPlot.subplots(1, 2, figsize=(6.8, 3.2))

    ax_h = axes[1]
    ax_h.axvline(0, color="grey", alpha=0.4)
    safe_hist!(ax_h, finite_values(data.H); bins=DEFAULT_SENSIBLE_BINS, color="C0")
    ax_h.set_xlabel(L"\overline{w'T'}~\mathrm{[W~m^{-2}]}")
    ax_h.set_ylabel("density")
    ax_h.tick_params(axis="y", labelleft=false)
    ax_h.set_title("Sensible")
    add_hours_text!(ax_h, "$(round(data_hours(data.H, avg_period), digits=1)) h")
    ax_h.grid(alpha=0.3)

    ax_le = axes[2]
    ax_le.axvline(0, color="grey", alpha=0.4)
    safe_hist!(ax_le, finite_values(data.LE); bins=DEFAULT_LATENT_BINS, color="C1")
    ax_le.set_xlabel(L"\overline{w'q'}~\mathrm{[W~m^{-2}]}")
    ax_le.tick_params(axis="y", labelleft=false)
    ax_le.set_title("Latent")
    add_hours_text!(ax_le, "$(round(data_hours(data.LE, avg_period), digits=1)) h")
    ax_le.grid(alpha=0.3)

    !isempty(title_suffix) && fig.suptitle(title_suffix)
    PyPlot.tight_layout()
    mkpath(dirname(output_file))
    PyPlot.savefig(output_file, bbox_inches="tight")
    PyPlot.close(fig)
    return output_file
end

subset_surface(data::DataFrame, surface_type::AbstractString) =
    surface_type == "total" ? data : data[data.flux_surface_type .== surface_type, :]

function plot_surface_histogram(data::DataFrame, output_file::AbstractString, avg_period::Period;
                                title_suffix="")
    fig, axes = PyPlot.subplots(5, 2, figsize=(8.4, 12.0), sharex="col")

    row_specs = [
        ("total", "total"),
        ("lead", "lead"),
        ("pond", "pond"),
        ("ice", "ice"),
        ("ridge", "ridge"),
    ]

    for (row_index, (label, surface_type)) in enumerate(row_specs)
        row_data = subset_surface(data, surface_type)
        ax_h = axes[row_index, 1]
        ax_le = axes[row_index, 2]

        data_1m = row_data[row_data.height_band .== "1m", :]
        data_2m = row_data[row_data.height_band .== "2m", :]

        ax_h.axvline(0, color="grey", alpha=0.4)
        safe_hist!(ax_h, finite_values(data_1m.H); bins=DEFAULT_SENSIBLE_BINS, color="C0", alpha=0.60, label="~1m")
        safe_hist!(ax_h, finite_values(data_2m.H); bins=DEFAULT_SENSIBLE_BINS, color="C1", alpha=0.60, label="~2m")
        add_hours_text!(
            ax_h,
            "1m: $(round(data_hours(data_1m.H, avg_period), digits=1)) h\n2m: $(round(data_hours(data_2m.H, avg_period), digits=1)) h",
            x=0.24,
        )
        ax_h.set_ylabel(label, fontsize=11, fontweight="bold", rotation=0, labelpad=36, va="center")
        ax_h.tick_params(axis="y", labelleft=false)
        ax_h.grid(alpha=0.3)
        row_index == 1 && ax_h.legend(fontsize=8)

        ax_le.axvline(0, color="grey", alpha=0.4)
        safe_hist!(ax_le, finite_values(row_data.LE); bins=DEFAULT_LATENT_BINS, color="C2")
        add_hours_text!(ax_le, "$(round(data_hours(row_data.LE, avg_period), digits=1)) h", x=0.20)
        ax_le.tick_params(axis="y", labelleft=false)
        ax_le.grid(alpha=0.3)

        if row_index == length(row_specs)
            ax_h.set_xlabel(L"\overline{w'T'}~\mathrm{[W~m^{-2}]}")
            ax_le.set_xlabel(L"\overline{w'q'}~\mathrm{[W~m^{-2}]}")
        end
    end

    axes[1, 1].set_title("Sensible")
    axes[1, 2].set_title("Latent")
    !isempty(title_suffix) && fig.suptitle(title_suffix)

    PyPlot.tight_layout()
    mkpath(dirname(output_file))
    PyPlot.savefig(output_file, bbox_inches="tight")
    PyPlot.close(fig)
    return output_file
end

function output_suffix(group_value)
    group_string = lowercase(replace(String(group_value), r"[^A-Za-z0-9]+" => "_"))
    group_string = strip(group_string, ['_'])
    isempty(group_string) || group_string == "all" ? "" : "_$group_string"
end

function plot_histograms(data::DataFrame, plot_dir::AbstractString, avg_period::Period)
    outputs = String[]
    groups = sort(unique(String.(data.histogram_group)))

    for group_value in groups
        group_data = data[data.histogram_group .== group_value, :]
        suffix = output_suffix(group_value)
        title_suffix = group_value == "all" ? "" : group_value

        push!(outputs, plot_total_histogram(
            group_data,
            joinpath(plot_dir, "hist_all_fluxes$(suffix).pdf"),
            avg_period;
            title_suffix=title_suffix,
        ))
        push!(outputs, plot_surface_histogram(
            group_data,
            joinpath(plot_dir, "hist_surface_types$(suffix).pdf"),
            avg_period;
            title_suffix=title_suffix,
        ))
    end

    return outputs
end


function finite_numeric_values(values)
    valid = Float64[]
    for value in values
        if !ismissing(value) && value isa Number
            value_float = Float64(value)
            isfinite(value_float) && push!(valid, value_float)
        end
    end
    return valid
end

function is_dship_wind_direction_column(column_name::AbstractString)
    normalized = lowercase(column_name)
    return normalized in ("winddirection", "wind_direction", "wind_dir", "true_wind_direction") ||
           endswith(normalized, "_wind_direction")
end

function dship_interval_mean(values, column_name::AbstractString)
    valid = finite_numeric_values(values)
    isempty(valid) && return NaN
    return is_dship_wind_direction_column(column_name) ? turb.mean_winddir(valid) : mean(valid)
end

function aggregate_dship_meteo_to_flux_times(
        dship_meteo::DataFrame, flux_times::AbstractVector, avg_interval::Period)
    unique_flux_times = unique(flux_times)
    aggregate_by_time = DataFrame(time=unique_flux_times)

    for column_name in names(dship_meteo)
        column_name == "time" && continue

        aggregate_by_time[!, Symbol(column_name)] = [
            begin
                interval_start = flux_time - avg_interval / 2
                interval_end = flux_time + avg_interval / 2
                in_interval = (interval_start .<= dship_meteo.time) .& (dship_meteo.time .< interval_end)
                dship_interval_mean(dship_meteo[in_interval, column_name], column_name)
            end
            for flux_time in unique_flux_times
        ]
    end

    row_by_time = Dict(flux_time => i for (i, flux_time) in pairs(unique_flux_times))
    aggregated = DataFrame(time=flux_times)
    for column_name in names(aggregate_by_time)
        column_name == "time" && continue
        values_by_time = aggregate_by_time[!, column_name]
        aggregated[!, Symbol(column_name)] = [values_by_time[row_by_time[flux_time]] for flux_time in flux_times]
    end

    return aggregated
end

function add_dship_meteo_columns!(flux_data::DataFrame, dship_meteo::DataFrame, avg_interval::Period)
    dship_aggregated = aggregate_dship_meteo_to_flux_times(dship_meteo, flux_data.time, avg_interval)
    for column_name in names(dship_aggregated)
        column_name == "time" && continue
        flux_data[!, Symbol(column_name)] = dship_aggregated[!, column_name]
    end
    return flux_data
end

## Settings
# Run this section first, then either rebuild from high-frequency data or read the cache.

station_names = selected_station_names() # or ["1a", "1b_1"]
avg_period = Second(400)
group_name = "all" # "all" or "daynight"; extend custom_grouping_labels for more

station_configs = load_station_configs(station_names)
cache_file = default_cache_file(station_configs, avg_period)
plot_dir = default_plot_dir(station_configs)

## Loading
# This loads one station at a time, similar to load_data.jl. Change station_ix in the REPL.
#=
station_ix = 1
station_config = station_configs[station_ix]
evaldfs = load_station_data(station_config)

## Flux calculation

fx_raws = calculate_fluxes(evaldfs, avg_period)

## Averaging to low-frequency data

fx_avgs = average_fluxes(fx_raws, avg_period)
low_fluxes = low_frequency_fluxes(fx_avgs, avg_period)
station_flux_data = station_flux_records(station_config, low_fluxes)
flux_data = copy(station_flux_data)
=#
## Accumulate all stations without keeping all high-frequency data in memory
# Run this section instead of the one-station sections above when rebuilding the full cache.
#=
flux_data = DataFrame()
@showprogress for station_config_loop in station_configs
    evaldfs = load_station_data(station_config_loop)
    fx_raws = calculate_fluxes(evaldfs, avg_period)
    fx_avgs = average_fluxes(fx_raws, avg_period)
    low_fluxes = low_frequency_fluxes(fx_avgs, avg_period)
    append!(flux_data, station_flux_records(station_config_loop, low_fluxes); cols=:union)

    evaldfs = nothing
    fx_raws = nothing
    fx_avgs = nothing
    low_fluxes = nothing
    GC.gc()
end

## Saving or loading low-frequency data

write_flux_cache(cache_file, flux_data)
=#
# To reuse the expensive result later, run the settings section and then this line:
flux_data = read_flux_cache(cache_file)

## Add meteo data

datatypes = [DateTime, Float32, Float32, Int32, Int32, Float32, Float32, Float32, Float32, Float32, Float32, Int32, Float32]
filename = "/home/haugened/Documents/data/CONTRASTS/weather_raw_data/meteo_dship_extracted_260619/meteo_dship_260619.dat"

dship_meteo = read_dship_meteo_csv(filename, datatypes) #no worries, there will be warnings for parsing "#" (missing data)

dship_meteo = dship_meteo[flux_data.time[1]-avg_period/2 .<= dship_meteo.time .<= flux_data.time[end]+avg_period/2, :]

#aggregate dship_meteo to flux_data time
dship_meteo_agg = aggregate_dship_meteo_to_flux_times(dship_meteo, flux_data.time, avg_period)

add_dship_meteo_columns!(flux_data, dship_meteo_agg, avg_period)

## Plotting

plot_outputs = plot_histograms(flux_data, plot_dir, avg_period)
