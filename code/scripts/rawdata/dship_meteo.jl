######################################################
###       INVESTIGATE THE DSHIP METEO DATA         ###
###            author: Michi Haugeneder            ###
######################################################
#=
Use to test the fog/clear sky/overcast&cloudy classification schemes.
Used in fluxes_histogram.jl
=#

using Dates, DataFrames, Statistics, LaTeXStrings, ProgressMeter
import CSV, PyCall, PyPlot
pydates = PyCall.pyimport("matplotlib.dates")
pypatches = PyCall.pyimport("matplotlib.patches")

importdir = joinpath(@__DIR__, "..", "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "station_config.jl"))
import .turb
import .gen
import .stationcfg
PyPlot.pygui(true)

"""
    read_dship_meteo_csv(filename::String, datatypes::Vector{DataType}, dateformat_inp=DateFormat("yyyymmddTHHMMSS"))

Read a CSV file extracted from DSHIP.

# Arguments
- `filename::String`: Path to the CSV file

# Returns
- DataFrame with the data
"""
function read_dship_meteo_csv(filename::String, datatypes::Vector{DataType}, dateformat_inp=DateFormat("yyyymmddTHHMMSS"))
    a = CSV.read(filename, DataFrame, header=1, skipto=4, dateformat=dateformat_inp,
                 types=datatypes, missingstring="#", silencewarnings=true)

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

"""
    classify_dship_weather(row)::Symbol

Classify DSHIP meteo records into broad visual-weather classes.

The DSHIP cloud-base field uses `99999` for no detected ceiling and `-9999`
for invalid data. Fog follows the standard visibility threshold (< 1000 m).
"""
function classify_dship_weather(row)::Symbol
    vis = row.visibility
    ceiling = row.ceiling_m

    if ismissing(vis) || ismissing(ceiling) || ceiling < 0
        return :unknown
    elseif vis <= 2000
        return :foggy
    elseif vis >= 10_000 && ceiling >= 90_000
        return :clear_sky
    elseif vis >= 5000 && ceiling >= 100
        return :overcast_cloudy
    else
        return :transient
    end
end

function add_dship_weather_classification!(dship_meteo::DataFrame)
    dship_meteo.weather_class = [classify_dship_weather(row) for row in eachrow(dship_meteo)]
    return dship_meteo
end

## Add meteo data

datatypes = [DateTime, Float32, Float32, Int32, Int32, Float32, Float32, Float32, Float32, Float32, Float32, Int32, Float32]
filename = "/home/haugened/Documents/data/CONTRASTS/weather_raw_data/meteo_dship_extracted_260619/meteo_dship_260619.dat"

dship_meteo = read_dship_meteo_csv(filename, datatypes) #no worries, there will be warnings for parsing "#" (missing data)

starttime = DateTime(2025, 07, 09, 15, 10, 00)
endtime = DateTime(2025, 08, 27, 15, 11, 40)
dship_meteo = dship_meteo[starttime .<= dship_meteo.time .<= endtime, :]
add_dship_weather_classification!(dship_meteo)

##
#plot DSHIP meteo time series with weather-class background
class_order = [:clear_sky, :overcast_cloudy, :foggy, :transient, :unknown]
class_labels = Dict(
    :clear_sky => "clear sky",
    :overcast_cloudy => "overcast/cloudy",
    :foggy => "foggy",
    :transient => "transient",
    :unknown => "unknown",
)
class_colors = Dict(
    :clear_sky => "tab:orange",
    :overcast_cloudy => "tab:green",
    :foggy => "tab:grey",
    :transient => "tab:red",
    :unknown => "white",
)
class_alpha = 0.18
class_percentages = Dict(
    class_name => 100 * count(==(class_name), dship_meteo.weather_class) / nrow(dship_meteo)
    for class_name in class_order
)

visibility_plot = [
    ismissing(v) || (v isa AbstractFloat && isnan(v)) ? NaN : Float64(v)
    for v in dship_meteo.visibility
]
ceiling_plot = [
    ismissing(v) || v <= -9000 ? NaN : min(Float64(v), 12000.0)
    for v in dship_meteo.ceiling_m
]
global_radiation_plot = [
    ismissing(v) || (v isa AbstractFloat && isnan(v)) ? NaN : Float64(v)
    for v in dship_meteo.global_radiation
]
rel_humidity_plot = [
    ismissing(v) || (v isa AbstractFloat && isnan(v)) ? NaN : Float64(v)
    for v in dship_meteo.rel_humidity
]

fig_weather, axs_weather = PyPlot.subplots(4, 1, figsize=(12, 9), sharex=true)
axs_weather = vec(axs_weather)

time_num = pydates.date2num.(dship_meteo.time)
dt_num = length(time_num) > 1 ? median(diff(time_num)) : 1 / (24 * 60)
span_start = first(time_num) - dt_num / 2
span_class = first(dship_meteo.weather_class)

for i in 2:length(time_num)
    global span_start, span_class

    if dship_meteo.weather_class[i] != span_class
        span_end = (time_num[i - 1] + time_num[i]) / 2
        for ax_tmp in axs_weather
            ax_tmp.axvspan(span_start, span_end,
                           color=get(class_colors, span_class, class_colors[:unknown]),
                           alpha=class_alpha, lw=0, zorder=0)
        end
        span_start = span_end
        span_class = dship_meteo.weather_class[i]
    end
end

span_end = last(time_num) + dt_num / 2
for ax_tmp in axs_weather
    ax_tmp.axvspan(span_start, span_end,
                   color=get(class_colors, span_class, class_colors[:unknown]),
                   alpha=class_alpha, lw=0, zorder=0)
end

axs_weather[1].plot(dship_meteo.time, visibility_plot, color="black", lw=0.8, zorder=2)
axs_weather[1].set_ylabel(L"visibility~\mathrm{[m]}")
axs_weather[1].set_yscale("log")
axs_weather[1].set_ylim(10, 80_000)
axs_weather[1].axhline(1000, color="black", ls="--", lw=0.8, alpha=0.45)
axs_weather[1].axhline(10_000, color="black", ls=":", lw=0.8, alpha=0.45)

axs_weather[2].plot(dship_meteo.time, ceiling_plot, color="black", lw=0.8, zorder=2)
axs_weather[2].set_ylabel(L"ceiling~\mathrm{[m]}")
axs_weather[2].set_ylim(0, 5200)
#axs_weather[2].axhline(5000, color="black", ls=":", lw=0.8, alpha=0.45)

axs_weather[3].plot(dship_meteo.time, global_radiation_plot, color="black", lw=0.8, zorder=2)
axs_weather[3].set_ylabel(L"global~radiation~\mathrm{[W~m^{-2}]}")
axs_weather[3].set_ylim(0,500)

axs_weather[4].plot(dship_meteo.time, rel_humidity_plot, color="black", lw=0.8, zorder=2)
axs_weather[4].set_ylabel(L"relative~humidity~\mathrm{[\%]}")
axs_weather[4].set_ylim(70, 105)
axs_weather[4].set_xlabel("date")

date_format = pydates.DateFormatter("%d.%m.")
majorloc = pydates.DayLocator(interval=5)
minorloc = pydates.DayLocator(interval=1)
for ax_tmp in axs_weather
    ax_tmp.xaxis_date()
    ax_tmp.xaxis.set_major_locator(majorloc)
    ax_tmp.xaxis.set_minor_locator(minorloc)
    ax_tmp.xaxis.set_major_formatter(date_format)
    ax_tmp.set_xlim(first(dship_meteo.time), last(dship_meteo.time))
    ax_tmp.grid(true, alpha=0.25, zorder=1)
end

legend_handles = [
    pypatches.Patch(facecolor=class_colors[class_name], alpha=class_alpha,
                    label=class_labels[class_name])
    for class_name in class_order
]
axs_weather[1].legend(handles=legend_handles, loc="center left", bbox_to_anchor=(1.01, 0.5))
percentage_text = join([
    string(class_labels[class_name], ": ",
           round(class_percentages[class_name]; digits=1), "%")
    for class_name in class_order
], "\n")
fig_weather.text(0.83, 0.58, percentage_text, ha="left", va="top", fontsize=10,
                 bbox=Dict("boxstyle" => "round", "facecolor" => "white",
                           "edgecolor" => "0.7", "alpha" => 0.9))
fig_weather.suptitle("DSHIP meteo time series and weather classification", fontweight="bold")
fig_weather.tight_layout(rect=(0, 0, 0.81, 0.96))
fig_weather.autofmt_xdate(bottom=0.06, rotation=30, ha="right")
fig_weather.subplots_adjust(right=0.81, bottom=0.07, top=0.93, hspace=0.08)
##

#aggregate dship_meteo to flux_data time
#dship_meteo_agg = aggregate_dship_meteo_to_flux_times(dship_meteo, flux_data.time, avg_period)

