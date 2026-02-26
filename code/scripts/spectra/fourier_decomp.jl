######################################################
###               PERFORMING FOURIER               ###
###                 DECOMPOSITION                  ###
###            author: Michi Haugeneder            ###
######################################################
#=
Performing Fourier spectral analysis for all 4 sonic anemometers.
According to R. Stull, "An Introduction to Boundary Layer Meteorology", 1988, p.307ff
Load data with load_data.jl first.
=#
using LsqFit, FFTW, Statistics, Dates
using PyPlot, PyCall, LaTeXStrings

importdir = joinpath(@__DIR__, "..", "..")

include(joinpath(importdir, "src", "fourier.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .ft
import .gen

######################################################
###              CHANGE VARIABLES HERE             ###
######################################################
log_avg_interval = 0.05  # log-averaging interval for plotting
mov_avg_window = 8       # moving average window for smoothing

######################################################

println()
println("-----------S-T-A-R-T-------------")

######################################################
###           COMPUTE AND PLOT SPECTRA             ###
######################################################

"""
    compute_spectrum(vec, timestep; log_interval=0.05, smooth_window=8)

Compute pre-multiplied spectral energy density f*S(f) for a vector.
Handles NaN values by filtering them out before processing.
Returns (freq, fSf) after log-averaging and smoothing.
"""
function compute_spectrum(vec::Vector, timestep::Period;
                          log_interval=0.05, smooth_window=8)
    # Filter NaN values
    valid = .!isnan.(vec)
    clean_vec = vec[valid]

    if length(clean_vec) < 30*20*60 #30 min
        return Float64[], Float64[]
    end

    # Duration in seconds based on number of valid samples
    dur = length(clean_vec) * Dates.value(timestep) / 1e3

    # Preprocessing: detrend and bell taper
    detrended = ft.detrend(clean_vec)
    tapered = ft.belltaper(detrended)

    # Fourier transform
    (FTraw, FTvec) = ft.fourierplus(tapered)
    (Soff, freq) = ft.spectralenergydensity(FTvec, dur)

    # Pre-multiply: f * S(f)
    fSf = freq .* Soff

    # Log-average for cleaner plotting
    (freq_avg, fSf_avg) = ft.logavg(freq, fSf, log_interval)

    # Smooth with moving average
    fSf_smooth = gen.movingaverage(fSf_avg, smooth_window)

    return freq_avg, fSf_smooth
end

##

println("Computing Fourier spectra and plotting...")
PyPlot.pygui(true)

# Panel layout: Row 1 = CSAT, Row 2 = IRG; Col 1 = T1, Col 2 = T2
# [evaldf2 evaldf4; evaldf1 evaldf3]
panel_dfs = [evaldf2 evaldf4; evaldf1 evaldf3]
panel_fxs = [fx2_raw fx4_raw; fx1_raw fx3_raw]  # flux DataFrames matching panel layout
panel_indices = [2 4; 1 3]  # indices into instr_labels, instr_type, etc.

# Timestep of flux data (from turb_fluxes.jl averaging)
fx_timestep = fx1.time[2] - fx1.time[1]

fig, axes = PyPlot.subplots(2, 2, figsize=(12, 8))
cmap = PyPlot.get_cmap("tab10")

for row in 1:2, col in 1:2
    ax = axes[row, col]
    df = panel_dfs[row, col]
    fx = panel_fxs[row, col]
    idx = panel_indices[row, col]

    title = "$(instr_labels[idx]) ($(instr_type[idx]))"

    # Variables to plot: u, w, T for all; h2o additionally for IRG
    vars = [("u", L"u"), ("w", L"w"), ("T", L"T_s")]
    if instr_type[idx] == "IRG"
        push!(vars, ("h2o", L"\rho_{H_2O}"))
    end

    for (i, (varname, varlabel)) in enumerate(vars)
        freq, fSf = compute_spectrum(df[!, varname], timestep;
                                     log_interval=log_avg_interval,
                                     smooth_window=mov_avg_window)
        if length(freq) > 0
            ax.plot(freq, fSf, color=cmap(i-1), linewidth=1.0, label=varlabel)
        end
    end

    #=
    # Flux spectrum wT for all sensors (from turb_fluxes.jl)
    freq_wT, fSf_wT = compute_spectrum(fx.wT, fx_timestep;
                                        log_interval=log_avg_interval,
                                        smooth_window=mov_avg_window)
    if length(freq_wT) > 0
        ax.plot(freq_wT, fSf_wT, color=cmap(2), linewidth=1.0, label=L"\overline{w'T'}")
    end

    # Flux spectrum uw for all sensors (from turb_fluxes.jl)
    freq_uw, fSf_uw = compute_spectrum(fx.uw, fx_timestep;
                                        log_interval=log_avg_interval,
                                        smooth_window=mov_avg_window)
    if length(freq_uw) > 0
        ax.plot(freq_uw, fSf_uw, color=cmap(3), linewidth=1.0, label=L"\overline{u'w'}")
    end

    # Flux spectrum wq for IRG sensors only (from turb_fluxes.jl)
    if instr_type[idx] == "IRG"
        freq_wq, fSf_wq = compute_spectrum(fx.wq, fx_timestep;
                                            log_interval=log_avg_interval,
                                            smooth_window=mov_avg_window)
        if length(freq_wq) > 0
            ax.plot(freq_wq, fSf_wq, color=cmap(4), linewidth=1.0, label=L"\overline{w'q'}")
        end
    end
    =#

    # Reference -2/3 slope (Kolmogorov)
    x_ref = collect(0.8:0.05:8)
    y_ref = x_ref .^ (-2 / 3)./10
    ax.plot(x_ref, y_ref, color="black", linestyle="--", linewidth=0.8,
            label=L"f^{-2/3}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(L"f~\mathrm{[Hz]}")
    ax.set_ylabel(L"fS(f)")
    ax.set_xlim(1e-4, 10)
    ax.set_ylim(1e-4, 2e-1)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(true, alpha=0.3)
    ax.legend(fontsize=8)
end

PyPlot.tight_layout()
PyPlot.gcf()
##
#PyPlot.savefig(joinpath("/home/haugened/Documents/data/CONTRASTS/plots/spectra/", "1a.pdf"), bbox_inches="tight")
