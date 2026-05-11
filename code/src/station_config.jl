module stationcfg

using Dates
using TOML

export available_station_names, load_station_config, optional_key, plot_dir,
       plot_root, require_key, station_datetime, station_file_stem,
       station_label, station_labels, toml_matrix

const DEFAULT_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "config", "stations"))

_key(key) = String(key)

function available_station_names(; config_dir::AbstractString=DEFAULT_CONFIG_DIR)
    isdir(config_dir) || return String[]
    names = String[]
    for filename in readdir(config_dir)
        if endswith(filename, ".toml")
            push!(names, first(splitext(filename)))
        end
    end
    return sort(names)
end

function load_station_config(station_name::AbstractString;
                             config_dir::AbstractString=DEFAULT_CONFIG_DIR)
    station_id = String(station_name)
    config_file = joinpath(config_dir, "$station_id.toml")
    if !isfile(config_file)
        available = join(available_station_names(config_dir=config_dir), ", ")
        error("Unknown station '$station_id'. Available station configs: $available")
    end

    config = TOML.parsefile(config_file)
    configured_id = get(config, "id", station_id)
    configured_id == station_id || error(
        "Config id '$configured_id' does not match file station '$station_id'."
    )
    config["config_file"] = config_file
    return config
end

function require_key(config::AbstractDict, keys...)
    isempty(keys) && return config

    value = config
    traversed = String[]
    for key in keys
        key_string = _key(key)
        push!(traversed, key_string)
        if !(value isa AbstractDict) || !haskey(value, key_string)
            error("Missing required station config key: $(join(traversed, "."))")
        end
        value = value[key_string]
    end
    return value
end

function optional_key(config::AbstractDict, default, keys...)
    value = config
    for key in keys
        key_string = _key(key)
        if !(value isa AbstractDict) || !haskey(value, key_string)
            return default
        end
        value = value[key_string]
    end
    return value
end

function station_datetime(config::AbstractDict, key::AbstractString)
    value = require_key(config, key)
    if value isa DateTime
        return value
    elseif value isa Date
        return DateTime(value)
    elseif value isa AbstractString
        return DateTime(replace(String(value), " " => "T"))
    end
    error("Station config key '$key' must be a DateTime or ISO datetime string.")
end

station_label(config::AbstractDict) = String(optional_key(config, require_key(config, "id"), "label"))
station_file_stem(config::AbstractDict) = String(optional_key(config, require_key(config, "id"), "file_stem"))

function station_labels(config::AbstractDict)
    configured_labels = optional_key(config, nothing, "instrument_labels")
    isnothing(configured_labels) || return String.(configured_labels)

    surface_type = String.(require_key(config, "surface_type"))
    heights = require_key(config, "heights")
    length(surface_type) == length(heights) || error(
        "surface_type and heights must have the same length."
    )
    return ["$(surface_type[i]) $(heights[i])m" for i in eachindex(surface_type)]
end

plot_root(config::AbstractDict) = String(
    optional_key(config, "/home/haugened/Documents/data/CONTRASTS/plots", "paths", "plot_root")
)

plot_dir(config::AbstractDict, parts::AbstractString...) = joinpath(plot_root(config), parts...)

function toml_matrix(value; T=Float64)
    value isa AbstractMatrix && return T.(value)
    rows = [T.(row) for row in value]
    isempty(rows) && return Matrix{T}(undef, 0, 0)
    return reduce(vcat, permutedims.(rows))
end

end
