######################################################
###           CALCULATE TURBULENT FLUXES           ###
###            author: Michi Haugeneder            ###
######################################################
#=
Calculate and plot turbulent fluxes and advection.
Load data with load_data.jl
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")

importdir = joinpath(@__DIR__, "..")
#datapath = "/home/haugened/Documents/slf/CONTRASTS25/data/2a/"
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .gen
PyPlot.pygui(true)

timestep = Millisecond(50)
ρ_air = 1.2 #kg m^{-3}
c_p = 1004 #J kg^{-1} K^{-1}
L_v = 2500e3 #J kg^{-1} (approx @0°C) 

######################################################
###              TURBULENT FLUXES                  ###
######################################################
turb.missing2nan!(evaldf1)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf3)
turb.missing2nan!(evaldf4)

#Reynolds averaging times (determine usind MRD)
ra1 = Second(400) #Millisecond(2^11*timestep) #Second(330)
ra2 = Second(400)
ra3 = Second(400) #Millisecond(2^11*timestep) #Second(320)
ra4 = Second(400)

fx1_raw = turb.turbflux(evaldf1, ra1, 1.0)
fx2_raw = turb.turbflux(evaldf2, ra2, 1.0)
fx3_raw = turb.turbflux(evaldf3, ra3, 1.0)
fx4_raw = turb.turbflux(evaldf4, ra4, 1.0)

#averaging
fx1 = turb.avgflux(fx1_raw, ra1, true, 0.1)
fx2 = turb.avgflux(fx2_raw, ra2, true, 0.1)
fx3 = turb.avgflux(fx3_raw, ra3, true, 0.1)
fx4 = turb.avgflux(fx4_raw, ra4, true, 0.1)

######################################################
###                    PLOTS                       ###
######################################################
##
#=
#flux time series
fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.set_title("2a - Turbulent heat fluxes at the lead")
pt1 = ax.plot(fx1.time[1:20*600:end], fx1.wT[1:20*600:end] .* (ρ_air * c_p), label="CSAT ice")
pt2 = ax.plot(fx2.time[1:20*60:end], fx2.wq[1:20*60:end].* (L_v * 1e-3), label="IRG ice")
pt3 = ax.plot(fx3.time[1:20*60:end], fx3.wT[1:20*60:end] .* (ρ_air * c_p), label="CSAT lead")
pt4 = ax.plot(fx4.time[1:20*60:end], fx4.wq[1:20*60:end].* (L_v * 1e-3), label="IRG lead")
#pt5 = ax.plot(fx5.time[1:20*60:end], fx5.wT[1:20*60:end], label="KAIJO")
#pt6 = ax.plot(fx6.time[1:20*60:end], fx6.wT[1:20*60:end], label="TJK")
#pt2wq = ax.plot(fx2.time[1:20*60:end], fx2.wq[1:20*60:end], label="w'q'(T2IRG)")
#pt1wq = ax.plot(fx1.time[1:20*600:end], fx1.wq[1:20*600:end], label="w'q'(T1IRG)")
#ax.set_ylabel(L"u_\ast~\mathrm{[m~s^{-1}]}")
#ax.set_ylabel(L"\overline{u'w'}~\mathrm{[m^2~s^{-2}]}")
#ax.set_ylabel(L"\overline{w'q'}~\mathrm{[m~g~s^{-1}~m^{-3}]}")
#ax.set_xlabel("31.05.2021")
ax.set_ylabel(L"\overline{w'q'} ~\mathrm{or}~\overline{w'T_s'}  ~\mathrm{[W~m^{-2}]}")
#ax.set_ylim([-0.25,0.25])
ax.xaxis_date()
#majorlocator = pydates.HourLocator(interval=1)
#minorlocator = pydates.MinuteLocator(interval=15)
#ax.xaxis.set_major_locator(majorlocator)
#ax.xaxis.set_minor_locator(minorlocator)
#date_format = pydates.DateFormatter("%H:%M")
#ax.xaxis.set_major_formatter(date_format)
#fig.autofmt_xdate()
#ax.set_xlim(DateTime(2021, 05, 31, 10, 30, 00), DateTime(2021, 05, 31, 20, 00, 00))
ax.grid()
ax.legend()
=#
##
######################################################
# Flux time series with dual subplots for latent and sensible heat fluxes (CONTRASTS)
##
fig, (ax1, ax2) = PyPlot.subplots(2, 1, figsize=(10, 8), sharex=true)

# Define colors for consistent legend
colors = ["C0", "C1", "C2", "C3"]  # Default matplotlib color cycle

#define plotting step size
step = 20*60 #every 1min

#y-axis limits
wT_limits = (-30, 30)
wq_limits = (-15,30) #wT_limits

# Upper subplot - Buoyancy fluxes (sensible heat)
ax1.set_title("1a Turbulent Heat Fluxes")
wt1 = ax1.plot(fx1.time[1:step:end], fx1.wT[1:step:end] .* (ρ_air * c_p), color=colors[1])
wt2 = ax1.plot(fx2.time[1:step:end], fx2.wT[1:step:end] .* (ρ_air * c_p), color=colors[2])
wt3 = ax1.plot(fx3.time[1:step:end], fx3.wT[1:step:end] .* (ρ_air * c_p), color=colors[3])
wt4 = ax1.plot(fx4.time[1:step:end], fx4.wT[1:step:end] .* (ρ_air * c_p), color=colors[4])
ax1.set_ylabel(L"\overline{w'T_s'} ~\mathrm{[W~m^{-2}]}")
ax1.grid()
ax1.set_ylim(wT_limits)

# Lower subplot - Latent heat fluxes
wq1 = ax2.plot(fx1.time[1:step:end], fx1.wq[1:step:end] .* (L_v * 1e-3), color=colors[1])
wq2 = ax2.plot(fx3.time[1:step:end], fx3.wq[1:step:end] .* (L_v * 1e-3), color=colors[3])
ax2.set_ylabel(L"\overline{w'q'} ~\mathrm{[W~m^{-2}]}")
ax2.set_xlabel("Time")
ax2.xaxis_date()
ax2.grid()
ax2.set_ylim(wq_limits)

# Create a single legend for the entire figure
handles = [wt1[1], wt2[1], wt3[1], wt4[1]]  # Get line objects
labels = ["pond 1.1m", "pond 2.1m", "ice 1.3m", "ice 2.1m"]
ax1.legend(handles, labels)#, loc="upper right", bbox_to_anchor=(1.0, 1))

# Optional: Uncomment these lines if you want to set specific time limits
# ax2.set_xlim(DateTime(2021, 05, 31, 10, 30, 00), DateTime(2021, 05, 31, 20, 00))

# Format dates on x-axis if needed
# majorlocator = pydates.HourLocator(interval=1)
# minorlocator = pydates.MinuteLocator(interval=15)
# ax2.xaxis.set_major_locator(majorlocator)
# ax2.xaxis.set_minor_locator(minorlocator)
# date_format = pydates.DateFormatter("%H:%M")
# ax2.xaxis.set_major_formatter(date_format)
# fig.autofmt_xdate()

PyPlot.tight_layout()
##
output_folder = "/home/haugened/Documents/data/CONTRASTS/plots/wT_wq/"
PyPlot.savefig(joinpath(output_folder, "1a.pdf"), bbox_inches="tight")
######################################################
#=
##
#understanding fluxes - plot components
plotdf = evaldf6
fluxdf = fx6
fluxrawdf = fx6_raw
step = 1
ra1idx = round(Int, Millisecond(ra1) / Millisecond(50))

cmap = PyPlot.get_cmap("tab10")
fig = PyPlot.figure()
gs = gridspec.GridSpec(6, 1)
ax5 = fig.add_subplot(gs[6, 1])
ax5.set_ylabel(L"T~\mathrm{[^\circ C]}")
#ax5.set_xlim([DateTime(2021,05,31,13,46,00), DateTime(2021,05,31,13,51,20)])
ax5.plot(plotdf.time[1:step:end], plotdf.T[1:step:end], color=cmap(0))
ax5.plot(plotdf.time[1:step:end], gen.movingaverage(plotdf.T, ra1idx)[1:step:end], color=cmap(1), alpha=0.8)
#minorlocator = pydates.MinuteLocator(interval=1)
#ax5.xaxis.set_minor_locator(minorlocator)
ax5.grid()
ax1 = fig.add_subplot(gs[1, 1], sharex=ax5)
ax1.set_title("T2LCSAT")
ax1.set_ylabel(L"\overline{w'T'}_{1min}~\mathrm{[K~m~s^{-1}]}")
ax1.plot(fluxrawdf.time[1:step:end], gen.movingaverage(fluxrawdf.wT, 20 * 60)[1:step:end], color=cmap(4), label=L"1~\mathrm{min~movavg}")
#ax1.plot(fluxrawdf.time[1:step:end], gen.movingaverage(fluxrawdf.wT[1:step:end],200)[1:step:end], color=cmap(2), label=L"10~\mathrm{s~movavg}", alpha=0.6)
#ax1.plot(fluxdf.time[1:step:end], fluxdf.wT[1:step:end], color=cmap(3), label=L"10~\mathrm{s~movavg}", alpha = 0.5)
ax1.legend()
ax1.grid()
ax1.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(gs[2, 1], sharex=ax5)
ax2.set_ylabel(L"u~\mathrm{[m~s^{-1}]}")
ax2.plot(plotdf.time[1:step:end], plotdf.u[1:step:end], color=cmap(0))
ax2.plot(plotdf.time[1:step:end], gen.movingaverage(plotdf.u, ra1idx)[1:step:end], color=cmap(1), alpha=0.8)
ax2.grid()
ax2.tick_params(axis="x", labelbottom=false)
ax3 = fig.add_subplot(gs[3, 1], sharex=ax5)
ax3.set_ylabel(L"v~\mathrm{[m~s^{-1}]}")
ax3.plot(plotdf.time[1:step:end], plotdf.v[1:step:end], color=cmap(0))
ax3.plot(plotdf.time[1:step:end], gen.movingaverage(plotdf.v, ra1idx)[1:step:end], color=cmap(1), alpha=0.8)
ax3.grid()
ax3.tick_params(axis="x", labelbottom=false)
ax3b = fig.add_subplot(gs[4, 1], sharex=ax5)
ax3b.set_ylabel(L"\mathrm{dev.~360^\circ~[^\circ]}")
ax3b.plot(plotdf.time[1:step:end], rad2deg.(atan.(-plotdf.v, plotdf.u))[1:step:end], ".", lw=1, color=cmap(0))
ax3b.grid()
ax3b.tick_params(axis="x", labelbottom=false)
ax4 = fig.add_subplot(gs[5, 1], sharex=ax5)
ax4.set_ylabel(L"w~\mathrm{[m~s^{-1}]}")
ax4.plot(plotdf.time[1:step:end], plotdf.w[1:step:end], color=cmap(0))
ax4.plot(plotdf.time[1:step:end], gen.movingaverage(plotdf.w, ra1idx)[1:step:end], color=cmap(1), alpha=0.8)
ax4.grid()
ax4.tick_params(axis="x", labelbottom=false)
PyPlot.show()
##
=#
######################################################
#4-window fluxes
#=
fig = PyPlot.figure()
gs = gridspec.GridSpec(4, 1)
axl = fig.add_subplot(gs[4, 1])
axl.grid()
axl.plot(fx1.time[1:1200:end], fx1.wT[1:1200:end])
axl.set_ylabel(L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}")
ax1 = fig.add_subplot(gs[1, 1], sharex=axl)
ax1.grid()
ax1.plot(fx2.time[1:1200:end], fx2.wT[1:1200:end])
ax1.set_ylabel(L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}")
ax1.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(gs[2, 1], sharex=axl)
ax2.grid()
ax2.plot(fx3.time[1:1200:end], fx3.wT[1:1200:end])
ax2.set_ylabel(L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}")
ax2.tick_params(axis="x", labelbottom=false)
ax3 = fig.add_subplot(gs[3, 1], sharex=axl)
ax3.grid()
ax3.plot(fx4.time[1:1200:end], fx4.wT[1:1200:end])
ax3.set_ylabel(L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}")
ax3.tick_params(axis="x", labelbottom=false)
=#
######################################################
#=
#fluxes at different scales

scale1 = Second(600)
scale2 = Second(10)

fx1_s1 = turb.avgflux(turb.turbflux(evaldf1, scale1), scale1)
fx2_s1 = turb.avgflux(turb.turbflux(evaldf2, scale1), scale1)
fx3_s1 = turb.avgflux(turb.turbflux(evaldf3, scale1), scale1)
fx4_s1 = turb.avgflux(turb.turbflux(evaldf4, scale1), scale1)

fx1_s2 = turb.avgflux(turb.turbflux(evaldf1, scale2), scale2)
fx2_s2 = turb.avgflux(turb.turbflux(evaldf2, scale2), scale2)
fx3_s2 = turb.avgflux(turb.turbflux(evaldf3, scale2), scale2)
fx4_s2 = turb.avgflux(turb.turbflux(evaldf4, scale2), scale2)
##
fig = PyPlot.figure()
gs = gridspec.GridSpec(2, 1)
ax1 = fig.add_subplot(gs[1, 1])
ax1.set_ylabel(L"\overline{w'T'}_{10min-10s}~\mathrm{[K~m~s^{-1}]}")
ax1.plot(fx1_s1.time[1:10:end], fx1_s1.wT[1:10:end] .- gen.movingaverage(fx1_s2.wT, round(Int, 20 * scale1 / scale2))[1:10:end], label="T1IRG")
ax1.plot(fx2_s1.time[1:10:end], fx2_s1.wT[1:10:end] .- gen.movingaverage(fx2_s2.wT, round(Int, 20 * scale1 / scale2))[1:10:end], label="T2IRG")
ax1.plot(fx3_s1.time[1:10:end], fx3_s1.wT[1:10:end] .- gen.movingaverage(fx3_s2.wT, round(Int, 20 * scale1 / scale2))[1:10:end], label="T2UCSAT")
ax1.plot(fx4_s1.time[1:10:end], fx4_s1.wT[1:10:end] .- gen.movingaverage(fx4_s2.wT, round(Int, 20 * scale1 / scale2))[1:10:end], label="T2LCSAT")
ax1.grid()
ax1.legend()
ax2 = fig.add_subplot(gs[2, 1], sharex=ax1)
ax2.set_ylabel(L"\overline{w'T'}_{10s}~\mathrm{[K~m~s^{-1}]}")
ax2.plot(fx1_s2.time[1:10:end], fx1_s2.wT[1:10:end], label="T1IRG")
ax2.plot(fx2_s2.time[1:10:end], fx2_s2.wT[1:10:end], label="T2IRG")
ax2.plot(fx3_s2.time[1:10:end], fx3_s2.wT[1:10:end], label="T2LCSAT")
ax2.plot(fx4_s2.time[1:10:end], fx4_s2.wT[1:10:end], label="T2UCSAT")
ax2.grid()
#ax2.legend()
##
=#
######################################################
#=
#comparison u'w' and v'w'
plotstep = 120*20
fig = PyPlot.figure()
gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
fig.suptitle("Comparison u'w' and v'w' - TJK")
ax1 = fig.add_subplot(gs[1, 1])
ax1.set_ylabel(L"\overline{u'w'}~[\mathrm{m^2~s^{-2}}]")
ax1.plot(fx1.time[1:plotstep:end], fx1.uw[1:plotstep:end])
ax1.tick_params(axis="x", labelbottom=false)
ax1.grid()
ax1.legend()
ax2 = fig.add_subplot(gs[2, 1], sharex=ax1)
ax2.plot(fx3.time[1:plotstep:end], fx3.vw[1:plotstep:end])
ax2.set_ylabel(L"\overline{v'w'}~[\mathrm{m^2~s^{-2}}]")
#date_format = pydates.DateFormatter("%H:%M")
#ax2.xaxis.set_major_formatter(date_format)
ax2.grid()
#ax3 = fig.add_subplot(gs[3,1], sharex=ax1)
#ax3.set_ylabel(L"\overline{v'w'}/\overline{u'w'}")
#ax3.plot(fx6.time[1:plotstep:end], fx6.vw[1:plotstep:end]./fx6.uw[1:plotstep:end])
#ax3.grid()
#ax3.set_ylim(-1,1)
##
=#

######################################################
#=
#scatter plot with histograms for x and y
evaldf1condition = fx1
evaldf2condition = fx2

condition = evaldf1.u .< 0 .&& evaldf2.u .< 0 #.!(Time(08,00,00) .< Dates.Time.(evaldf1condition.time) .< Time(18,00,00))
titlepart = "upvalley"

evaldf1condition = evaldf1condition[condition, :]
evaldf2condition = evaldf2condition[condition, :]

#x = evaldf1condition.T - evaldf2condition.T
#y = (evaldf1condition.u + evaldf2condition.u)./2
x = evaldf1condition.wT
y = evaldf2condition.wT

#xlabel = L"\overline{\Delta T} = \overline{T_{IRG1}} - \overline{T_{IRG2}}~\mathrm{[K]}"
#ylabel = L"u~\mathrm{[m~s^{-1}]}"
xlabel = L"\overline{w'T'}_{T1IRG}~\mathrm{[K~m~s^{-1}]}"
ylabel = L"\overline{w'T'}_{T2IRG}~\mathrm{[K~m~s^{-1}]}"

xmin = -0.1
xstep = 0.0001
xmax = 0.4
ymin = -0.1
ystep = 0.0001
ymax = 0.4

PyPlot.pygui(true)
gig = PyPlot.figure(figsize=(10, 8))
gs = gridspec.GridSpec(2, 3, height_ratios=[5, 1], width_ratios=[1, 5, 0])
axg = gig.add_subplot(gs[1, 2])
axg.set_title(string(titlepart, " ", evaldf1.time[1], " - ", evaldf1.time[end]))
axg.grid()
sp2 = axg.scatter(x, y, s=0.2, alpha=0.5)#, c=pydates.date2num.(evaldf1valupwind.time[:]), cmap="turbo")
axg.set_xlim([xmin, xmax])
axg.set_ylim([ymin, ymax])
axg.tick_params(axis="x", labelbottom=false)
axg.tick_params(axis="y", labelleft=false)
#=caxg = gig.add_subplot(gs[1,3])
cb2 = PyPlot.colorbar(sp2, cax=caxg)
cb2.set_label("date/time")
loc = pydates.AutoDateLocator()
cb2.ax.yaxis.set_major_locator(loc)
date_format = pydates.DateFormatter("%d.%m %H")
cb2.ax.yaxis.set_major_formatter(date_format)
axgl = gig.add_subplot(gs[1, 1], sharey=axg)=#
axgl.hist(y, bins=collect(ymin:ystep:ymax), density=true, orientation="horizontal")
axgl.set_ylabel(ylabel)
axgl.tick_params(axis="x", labelbottom=false)
axgb = gig.add_subplot(gs[2, 2], sharex=axg)
axgb.set_xlabel(xlabel)
axgb.hist(x, bins=collect(xmin:xstep:xmax), density=true)
axgb.tick_params(axis="y", labelleft=false)
=#
######################################################
#Scatter plot fluxes CONTRASTS
mpl_scatter_density = pyimport("mpl_scatter_density")

scatter1b = fx1.wT[.!isnan.(fx1.wT) .&& .!isnan.(fx3.wT)]
scatter1a = fx3.wT[.!isnan.(fx1.wT) .&& .!isnan.(fx3.wT)]

scatter1a .*= ρ_air .* c_p
scatter1b .*= ρ_air .* c_p

scatter2b = fx1.wq[.!isnan.(fx1.wq) .&& .!isnan.(fx3.wq)]
scatter2a = fx3.wq[.!isnan.(fx1.wq) .&& .!isnan.(fx3.wq)]

scatter2a .*= L_v .* 1e-3
scatter2b .*= L_v .* 1e-3

scatter3b = fx2.wT[.!isnan.(fx2.wT) .&& .!isnan.(fx4.wT)]
scatter3a = fx4.wT[.!isnan.(fx2.wT) .&& .!isnan.(fx4.wT)]

scatter3a .*= ρ_air .* c_p
scatter3b .*= ρ_air .* c_p

ratio1 = scatter1b ./ scatter1a
ratio2 = scatter2b ./ scatter2a
ratio3 = scatter3b ./ scatter3a

##
fig = PyPlot.figure(figsize=(12, 9.5))
fig.suptitle("1b position 1 pond - 25.07.2025 to 27.07.2025")
gs = gridspec.GridSpec(3, 3, height_ratios=(2, 1, 1))

# --- Row 1: Scatter density (unchanged) ---
ax1 = fig.add_subplot(gs[1, 1], projection="scatter_density")
ax1.set_title("1m Sensible Heat Fluxes")
ax1.axhline(0, color="grey", alpha=0.6)
ax1.axvline(0, color="grey", alpha=0.6)
ax1.plot([-25,25], [-25,25], color="grey", alpha=0.3)
ax1.scatter_density(scatter1a, scatter1b, color="black", vmin=00, vmax=2000)
ax1.set_xlim(0, 10)
ax1.set_ylim(-5,10)
#ax1.set_yticks(collect(-10:5:10))
ax1.set_xlabel(L"\overline{w'T'}_{ice}~\mathrm{[W~m^{-2}]}")
ax1.set_ylabel(L"\overline{w'T'}_{pond}~\mathrm{[W~m^{-2}]}")
ax1.grid()
ax1.set_aspect("equal")

ax2 = fig.add_subplot(gs[1, 2], projection="scatter_density")
ax2.set_title("1m Latent Heat Fluxes")
ax2.axhline(0, color="grey")
ax2.axvline(0, color="grey")
ax2.plot([-25,25], [-25,25], color="grey", alpha=0.3)
ax2.scatter_density(scatter2a, scatter2b, color="blue", vmin=0, vmax=2000)
ax2.set_xlim(0,7)
ax2.set_ylim(0,7)
#ax2.set_yticks(collect(-10:5:10))
ax2.set_xlabel(L"\overline{w'q'}_{ice}~\mathrm{[W~m^{-2}]}")
ax2.set_ylabel(L"\overline{w'q'}_{pond}~\mathrm{[W~m^{-2}]}")
ax2.grid()
ax2.set_aspect("equal")

ax3 = fig.add_subplot(gs[1, 3], projection="scatter_density")
ax3.plot([-25,30], [-25,30], color="grey", alpha=0.3)
ax3.scatter_density(scatter3a, scatter3b, color="black", vmin=0, vmax=2000)
ax3.set_title("2m Sensible Heat Fluxes")
ax3.axhline(0, color="grey")
ax3.axvline(0, color="grey")
ax3.set_xlim(-10,10)
ax3.set_ylim(-10,10)
#ax3.set_yticks(collect(-5:5:10))
ax3.set_xlabel(L"\overline{w'T'}_{ice}~\mathrm{[W~m^{-2}]}")
ax3.set_ylabel(L"\overline{w'T'}_{pond}~\mathrm{[W~m^{-2}]}")
ax3.grid()
ax3.set_aspect("equal")

# --- Row 2: Overlaid histograms of individual fluxes ---
ax4 = fig.add_subplot(gs[2, 1])
ax4.hist(scatter1a, bins=collect(0:0.1:25), color="grey", alpha=0.5, density=true, label="ice")
ax4.hist(scatter1b, bins=collect(0:0.1:25), color="orange", alpha=0.5, density=true, label="pond")
ax4.set_xlabel(L"\overline{w'T'}~\mathrm{[W~m^{-2}]}")
ax4.set_ylabel("PDF")
ax4.legend()
ax4.grid(alpha=0.3)

ax5 = fig.add_subplot(gs[2, 2])
ax5.hist(scatter2a, bins=collect(0:0.1:25), color="grey", alpha=0.5, density=true, label="ice")
ax5.hist(scatter2b, bins=collect(0:0.1:25), color="orange", alpha=0.5, density=true, label="pond")
ax5.set_xlabel(L"\overline{w'q'}~\mathrm{[W~m^{-2}]}")
ax5.legend()
ax5.grid(alpha=0.3)

ax6 = fig.add_subplot(gs[2, 3])
ax6.hist(scatter3a, bins=collect(0:0.1:25), color="grey", alpha=0.5, density=true, label="ice")
ax6.hist(scatter3b, bins=collect(0:0.1:25), color="orange", alpha=0.5, density=true, label="pond")
ax6.set_xlabel(L"\overline{w'T'}~\mathrm{[W~m^{-2}]}")
ax6.legend()
ax6.grid(alpha=0.3)

# --- Row 3: Histograms of ratios ---
ax7 = fig.add_subplot(gs[3, 1])
ax7.hist(ratio1, bins=collect(0:0.01:2.0), color="black", alpha=0.7, density=true)
ax7.axvline(1, color="red", linestyle="--", alpha=0.3, label="1:1")
ax7.set_xlabel(L"\left(\overline{w'T'}_{pond} ~/~\overline{w'T'}_{ice}\right)_{IRG}")
ax7.set_ylabel("PDF")
ax7.grid(alpha=0.3)

ax8 = fig.add_subplot(gs[3, 2])
ax8.hist(ratio2, bins=collect(-1.0:0.01:2.0), color="blue", alpha=0.7, density=true)
ax8.axvline(1, color="red", linestyle="--", alpha=0.3, label="1:1")
ax8.set_xlabel(L"\left(\overline{w'q'}_{pond} ~/~ \overline{w'q'}_{ice}\right)_{CSAT}")
ax8.grid(alpha=0.3)

ax9 = fig.add_subplot(gs[3, 3])
ax9.hist(ratio3, bins=collect(-0.5:0.01:1.5), color="black", alpha=0.7, density=true)
ax9.axvline(1, color="red", linestyle="--", alpha=0.3, label="1:1")
ax9.set_xlabel(L"\left(\overline{w'T'}_{pond} ~/~\overline{w'T'}_{ice}\right)_{CSAT}")
ax9.grid(alpha=0.3)

PyPlot.tight_layout()
##
output_folder = "/home/haugened/Documents/data/CONTRASTS/plots/correlation/"
PyPlot.savefig(joinpath(output_folder, "1a.pdf"), bbox_inches="tight")

#calculating correlations
using NaNStatistics
cor1 = nancor(fx1.wT, fx3.wT)
cor2 = nancor(fx1.wq, fx3.wq)
cor3 = nancor(fx2.wT, fx4.wT)

#calculating mean bias
scatter1a_mean = nanmean(scatter1a)
scatter1b_mean = nanmean(scatter1b)
scatter2a_mean = nanmean(scatter2a)
scatter2b_mean = nanmean(scatter2b)
scatter3a_mean = nanmean(scatter3a)
scatter3b_mean = nanmean(scatter3b)

mad_irg_wT = scatter1a_mean - scatter1b_mean
mad_irg_wq = scatter2a_mean - scatter2b_mean
mad_csat_wT = scatter3a_mean - scatter3b_mean
######################################################