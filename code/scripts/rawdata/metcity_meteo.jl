####################################################################
###          READING AND EVALUATING DATA FROM METCITY            ###
####################################################################

using Dates, DataFrames, Statistics, LaTeXStrings
import CSV

const METCITY_PATH = "/home/haugened/Documents/data/CONTRASTS/MetCity"
const TURBULENCE_FLUX_FILE = "/home/haugened/Documents/data/CONTRASTS/EC_offline_preproc/cut/for_hist/fluxes_histogram_400s.csv"
const METCITY_DATE_FORMAT = dateformat"yyyy-mm-dd HH:MM:SS"

const FLAG_BASE_NAMES = Dict(
    "RH2" => "Tower_RH2",
    "RH6" => "Tower_RH6",
    "RH10" => "Tower_RH10",
    "LuftT2" => "Tower_AirT2",
    "LuftT6" => "Tower_AirT6",
    "LuftT10" => "Tower_AirT10",
    "TKE" => "TKE_mean",
)

const RADIATION_COLUMNS_TO_NEGATE = [
    :CNR4_sw_up,
    :CNR4_lw_up_cor,
    :CNR4_lw_up_uncor,
    :CNR4_sw_net_ice,
    :CNR4_lw_net_ice,
    :CNR4_net_ice,
]

function standardize_flag_names!(flags)
    rename!(flags, :datetime => :time)

    for name in names(flags)[2:end]
        base_name = replace(name, r"_QAQC_flag$|_Flag$|_flag$" => "")
        data_name = get(FLAG_BASE_NAMES, base_name, base_name)
        rename!(flags, name => "$(data_name)_Flag")
    end

    fill_missing_flags!(flags)
    disallowmissing!(flags)
    return flags
end

function fill_missing_flags!(data)
    for column in names(data, r"_Flag$")
        data[!, column] = coalesce.(data[!, column], -1)
    end

    return data
end

function apply_radiation_sign_convention!(data)
    for column in RADIATION_COLUMNS_TO_NEGATE
        data[!, column] .*= -1
    end

    return data
end

function read_metcity_station(directory, station)
    stem = "CONTRASTS_IceStation$(station)_Level0"

    data = CSV.read(
        joinpath(directory, "$(stem)_data.csv"),
        DataFrame;
        header=1,
        skipto=3,
        dateformat=METCITY_DATE_FORMAT,
        types=Dict(1 => DateTime),
    )
    rename!(data, names(data)[1:4] .=> [:time, :Tower_RH2, :Tower_RH6, :Tower_RH10])
    insertcols!(data, 2, :ice_station => fill(station, nrow(data)))
    apply_radiation_sign_convention!(data)

    flags = CSV.read(
        joinpath(directory, "$(stem)_flag.csv"),
        DataFrame;
        dateformat=METCITY_DATE_FORMAT,
        types=Dict(:datetime => DateTime),
    )
    standardize_flag_names!(flags)

    meteo = leftjoin(data, flags; on=:time)
    fill_missing_flags!(meteo)
    return meteo
end

function read_metcity_meteo(directory=METCITY_PATH)
    meteo = vcat(
        [read_metcity_station(directory, station) for station in 1:3]...;
        cols=:union,
    )
    sort!(meteo, :time)
    return meteo
end

function qc_timeseries(meteo, column, acceptable_flags)
    flag_column = Symbol(column, "_Flag")
    accepted = in.(meteo[!, flag_column], Ref(acceptable_flags))
    return select(meteo[accepted, :], :time, column)
end

function qc_dataframe(meteo, columns, acceptable_flags)
    accepted = reduce(.&, [
        in.(meteo[!, Symbol(column, "_Flag")], Ref(acceptable_flags))
        for column in columns
    ])
    return select(meteo[accepted, :], :time, :ice_station, columns)
end

function valid_turbulence_windows(filename=TURBULENCE_FLUX_FILE)
    fluxes = CSV.read(filename, DataFrame; types=Dict(:time => DateTime))
    fluxes.ice_station = parse.(Int, first.(fluxes.station))
    fluxes.time = round.(fluxes.time, Minute(10))

    valid_flux = isfinite.(fluxes.wT) .| isfinite.(fluxes.wq)
    return unique(select(fluxes[valid_flux, :], :ice_station, :time))
end

function filter_to_valid_turbulence(meteo, filename=TURBULENCE_FLUX_FILE)
    valid_windows = valid_turbulence_windows(filename)
    return semijoin(meteo, valid_windows; on=[:ice_station, :time])
end

metcity_meteo_all = read_metcity_meteo()
metcity_meteo = filter_to_valid_turbulence(metcity_meteo_all)

metcity_station1 = copy(metcity_meteo[metcity_meteo.ice_station .== 1, :])
metcity_station2 = copy(metcity_meteo[metcity_meteo.ice_station .== 2, :])
metcity_station3 = copy(metcity_meteo[metcity_meteo.ice_station .== 3, :])

#################################################################
#plotting
using PyPlot
PyPlot.pygui(true)

###########################
#albedo
##
albedo1_good = qc_timeseries(metcity_station1, "CNR4_albedo", [0,1])
albedo2_good = qc_timeseries(metcity_station2, "CNR4_albedo", [0,1])
albedo3_good = qc_timeseries(metcity_station3, "CNR4_albedo", [0,1])

fig=PyPlot.figure()
ax = fig.add_subplot(111)
ax.plot(metcity_station1.time, metcity_station1.CNR4_albedo, label="regime 1")
ax.plot(metcity_station2.time, metcity_station2.CNR4_albedo, label="regime 2")
ax.plot(metcity_station3.time, metcity_station3.CNR4_albedo, label="regime 3")
ax.legend()
ax.grid(true)
ax.set_xlabel("time")
ax.set_ylabel("CNR4 albedo")
ax.set_ylim(0.3,1.1)
#fig.savefig("/home/haugened/Documents/data/CONTRASTS/plots/MetCity/albedo.pdf", bbox_inches="tight")
##

###########################
#radiation components
##
radiation_columns = [
    ("CNR4_sw_dn", "shortwave upwelling"),
    ("CNR4_sw_up", "shortwave downwelling"),
    ("CNR4_lw_dn_cor", "longwave upwelling"),
    ("CNR4_lw_up_cor", "longwave downwelling"),
]
radiation_components_df = qc_dataframe(metcity_meteo, first.(radiation_columns), [0,1])
metcity_stations = [metcity_station1, metcity_station2, metcity_station3]

fig_radiation, axs_radiation = PyPlot.subplots(4, 1, figsize=(12, 10), sharex=true)
axs_radiation = vec(axs_radiation)

for (ax_radiation, (column, label)) in zip(axs_radiation, radiation_columns)
    for (regime, station) in enumerate(metcity_stations)
        radiation_good = qc_timeseries(station, column, [0,1])
        ax_radiation.plot(
            radiation_good.time,
            radiation_good[!, column],
            label="regime $(regime)",
        )
    end

    ax_radiation.set_ylabel("$(label)\n[W m⁻²]")
    ax_radiation.grid(true)
end

axs_radiation[1].legend()
axs_radiation[end].set_xlabel("time")
fig_radiation.suptitle("CNR4 radiation components")
fig_radiation.tight_layout()
#fig_radiation.savefig("/home/haugened/Documents/data/CONTRASTS/plots/MetCity/radiation_components.pdf", bbox_inches="tight")
##

###########################
#net radiation
##
net_radiation_columns = [
    ("CNR4_sw_net_ice", "net shortwave radiation"),
    ("CNR4_lw_net_ice", "net longwave radiation"),
    ("CNR4_net_ice", "total radiation budget"),
]
net_radiation_df = qc_dataframe(metcity_meteo, first.(net_radiation_columns), [0,1])

fig_net_radiation, axs_net_radiation = PyPlot.subplots(3, 1, figsize=(12, 8), sharex=true)
axs_net_radiation = vec(axs_net_radiation)

for (ax_net, (column, label)) in zip(axs_net_radiation, net_radiation_columns)
    for (regime, station) in enumerate(metcity_stations)
        radiation_good = qc_timeseries(station, column, [0,1])
        ax_net.plot(
            radiation_good.time,
            radiation_good[!, column],
            label="regime $(regime)",
        )
    end

    ax_net.axhline(0, color="black", linewidth=0.8)
    ax_net.set_ylabel("$(label)\n[W m⁻²]")
    ax_net.grid(true)
end

axs_net_radiation[1].legend()
axs_net_radiation[end].set_xlabel("time")
fig_net_radiation.suptitle("CNR4 net radiation")
fig_net_radiation.tight_layout()
#fig_net_radiation.savefig("/home/haugened/Documents/data/CONTRASTS/plots/MetCity/net_radiation.pdf", bbox_inches="tight")
##

###########################
#radiation and flux distributions
##
fluxes_for_comparison = CSV.read(
    TURBULENCE_FLUX_FILE,
    DataFrame;
    types=Dict(:time => DateTime),
)
fluxes_for_comparison.ice_station = parse.(Int, first.(fluxes_for_comparison.station))
fluxes_for_comparison.time = ceil.(fluxes_for_comparison.time, Minute(10))

flux_availability = combine(
    groupby(fluxes_for_comparison, [:ice_station, :time]),
    :H => (values -> any(isfinite, values)) => :has_H,
    :LE => (values -> any(isfinite, values)) => :has_LE,
)
common_flux_windows = flux_availability[
    flux_availability.has_H .& flux_availability.has_LE,
    [:ice_station, :time],
]

radiation_flux_windows = innerjoin(
    net_radiation_df,
    common_flux_windows;
    on=[:ice_station, :time],
)
radiation_flux_data = semijoin(
    fluxes_for_comparison,
    radiation_flux_windows;
    on=[:ice_station, :time],
)

radiation_flux_distributions = [
    ("net radiation", filter(isfinite, radiation_flux_windows.CNR4_net_ice)),
    ("sensible heat flux", filter(isfinite, radiation_flux_data.H)),
    ("latent heat flux", filter(isfinite, radiation_flux_data.LE)),
]

radiation_flux_statistics = DataFrame(
    variable=first.(radiation_flux_distributions),
    count=length.(last.(radiation_flux_distributions)),
    mean=mean.(last.(radiation_flux_distributions)),
    median=median.(last.(radiation_flux_distributions)),
    standard_deviation=std.(last.(radiation_flux_distributions)),
    minimum=minimum.(last.(radiation_flux_distributions)),
    maximum=maximum.(last.(radiation_flux_distributions)),
)

fig_radiation_flux, axs_radiation_flux = PyPlot.subplots(1, 3, figsize=(12, 4))
axs_radiation_flux = vec(axs_radiation_flux)
xlims = [(-150,50),(-50,50), (-30,30)]
axlabels = [L"\mathrm{net~radiation~[W~m^{-2}]}", L"Q_H~\mathrm{[W~m^{-2}]}", L"Q_E~\mathrm{[W~m^{-2}]}"]

for (ax_hist, (label, values), xlim, axlabel) in zip(axs_radiation_flux, radiation_flux_distributions, xlims, axlabels)
    total_width = last(xlim) - first(xlim)
    width_per_bin = total_width/100
    ax_hist.hist(values, bins=collect(first(xlim):width_per_bin:last(xlim)), density=true, alpha=0.8)
    ax_hist.axvline(mean(values), color="black", linestyle="--", label="mean")
    ax_hist.axvline(median(values), color="black", linestyle=":", label="median")
    #ax_hist.set_title(label)
    ax_hist.set_xlim(xlim)
    ax_hist.set_xlabel(axlabel)
    ax_hist.grid(true)
end

axs_radiation_flux[1].set_ylabel("probability density")
axs_radiation_flux[1].legend()
fig_radiation_flux.tight_layout()
#fig_radiation_flux.savefig("/home/haugened/Documents/data/CONTRASTS/plots/MetCity/radiation_flux_histograms.pdf", bbox_inches="tight")

###########################
#net radiation and turbulent heat flux
##
selected_irgason = (
    (fluxes_for_comparison.instrument_type .== "IRG") .&
    (
        ((fluxes_for_comparison.ice_station .< 3) .&
         (fluxes_for_comparison.surface_type .== "ice")) .|
        ((fluxes_for_comparison.ice_station .== 3) .&
         (fluxes_for_comparison.surface_type .!= "floe"))
    ) .&
    isfinite.(fluxes_for_comparison.H) .&
    isfinite.(fluxes_for_comparison.LE)
)

selected_irgason_fluxes = copy(fluxes_for_comparison[selected_irgason, :])
selected_irgason_fluxes.turbulent_flux_sum = (
    selected_irgason_fluxes.H .+ selected_irgason_fluxes.LE
)

turbulent_flux_sum_10min = combine(
    groupby(selected_irgason_fluxes, [:station, :ice_station, :time]),
    :turbulent_flux_sum => mean => :turbulent_flux_sum,
)

radiation_turbulent_comparison_df = innerjoin(
    net_radiation_df,
    turbulent_flux_sum_10min;
    on=[:ice_station, :time],
)

radiation_turbulent_distribution = [
    ("net radiation", radiation_turbulent_comparison_df.CNR4_net_ice),
    ("sensible + latent heat flux", radiation_turbulent_comparison_df.turbulent_flux_sum),
]

radiation_turbulent_statistics = DataFrame(
    variable=first.(radiation_turbulent_distribution),
    count=length.(last.(radiation_turbulent_distribution)),
    mean=mean.(last.(radiation_turbulent_distribution)),
    median=median.(last.(radiation_turbulent_distribution)),
    standard_deviation=std.(last.(radiation_turbulent_distribution)),
    minimum=minimum.(last.(radiation_turbulent_distribution)),
    maximum=maximum.(last.(radiation_turbulent_distribution)),
)

fig_net_turbulent, axs_net_turbulent = PyPlot.subplots(1, 2, figsize=(8, 4))
axs_net_turbulent = vec(axs_net_turbulent)
net_turbulent_xlims = [(-150, 50), (-50, 50)]
net_turbulent_labels = [
    L"\mathrm{net~radiation~[W~m^{-2}]}",
    L"Q_H + Q_E~\mathrm{[W~m^{-2}]}",
]

for (ax_hist, (_, values), xlim, axlabel) in zip(
    axs_net_turbulent,
    radiation_turbulent_distribution,
    net_turbulent_xlims,
    net_turbulent_labels,
)
    bin_width = (last(xlim) - first(xlim)) / 100
    bins = collect(first(xlim):bin_width:last(xlim))
    ax_hist.hist(values, bins=bins, density=true, alpha=0.8)
    ax_hist.axvline(mean(values), color="black", linestyle="--", label="mean")
    ax_hist.axvline(median(values), color="black", linestyle=":", label="median")
    ax_hist.set_xlim(xlim)
    ax_hist.set_xlabel(axlabel)
    y_max = last(ax_hist.get_ylim())
    ax_hist.set_yticks(collect(range(0, y_max; length=6)))
    ax_hist.tick_params(axis="y", labelleft=false)
    ax_hist.grid(true)
end

axs_net_turbulent[1].set_ylabel("probability density")
axs_net_turbulent[1].legend()
fig_net_turbulent.tight_layout()
#fig_net_turbulent.savefig("/home/haugened/Documents/data/CONTRASTS/plots/MetCity/net_radiation_turbulent_flux.pdf", bbox_inches="tight")
##
