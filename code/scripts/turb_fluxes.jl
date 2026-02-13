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
pt1 = ax.plot(fx1.time[1:20*600:end], fx1.wT[1:20*600:end] .* (ρ_air * c_p), label = instr_labels[1])
pt2 = ax.plot(fx2.time[1:20*60:end], fx2.wq[1:20*60:end].* (L_v * 1e-3),     label = instr_labels[2])
pt3 = ax.plot(fx3.time[1:20*60:end], fx3.wT[1:20*60:end] .* (ρ_air * c_p),   label = instr_labels[3])
pt4 = ax.plot(fx4.time[1:20*60:end], fx4.wq[1:20*60:end].* (L_v * 1e-3),     label = instr_labels[4])
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
labels = instr_labels
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
output_folder_wT = output_folder = "/home/haugened/Documents/data/CONTRASTS/plots/wT_wq/"
#PyPlot.savefig(joinpath(output_folder_wT, "1b_2.pdf"), bbox_inches="tight")
##
output_folder = "/home/haugened/Documents/data/CONTRASTS/plots/wT_wq/"
PyPlot.savefig(joinpath(output_folder, "1a.pdf"), bbox_inches="tight")
##
# Momentum fluxes and friction velocity plot
fig, (ax1, ax2) = PyPlot.subplots(2, 1, figsize=(10, 8), sharex=true)

# Define colors for consistent legend
colors = ["C0", "C1", "C2", "C3"]

# Define plotting step size
step = 20*60  # every 1min

# Upper subplot - Vertical momentum fluxes
ax1.set_title("Vertical Momentum Fluxes and Friction Velocity")
uw1 = ax1.plot(fx1.time[1:step:end], fx1.uw[1:step:end], color=colors[1])
uw2 = ax1.plot(fx2.time[1:step:end], fx2.uw[1:step:end], color=colors[2])
uw3 = ax1.plot(fx3.time[1:step:end], fx3.uw[1:step:end], color=colors[3])
uw4 = ax1.plot(fx4.time[1:step:end], fx4.uw[1:step:end], color=colors[4])
ax1.set_ylabel(L"\overline{u'w'} ~\mathrm{[m^2~s^{-2}]}")
ax1.grid()

# Lower subplot - Friction velocity
us1 = ax2.plot(fx1.time[1:step:end], fx1.u_star[1:step:end], color=colors[1])
us2 = ax2.plot(fx2.time[1:step:end], fx2.u_star[1:step:end], color=colors[2])
us3 = ax2.plot(fx3.time[1:step:end], fx3.u_star[1:step:end], color=colors[3])
us4 = ax2.plot(fx4.time[1:step:end], fx4.u_star[1:step:end], color=colors[4])
ax2.set_ylabel(L"u_*~\mathrm{[m~s^{-1}]}")
ax2.set_xlabel("Time")
ax2.xaxis_date()
ax2.grid()

# Create a single legend for the entire figure
handles = [uw1[1], uw2[1], uw3[1], uw4[1]]
#labels same as above
ax1.legend(handles, labels)

PyPlot.tight_layout()
output_folder_uw = output_folder = "/home/haugened/Documents/data/CONTRASTS/plots/uw_ustar/"
#PyPlot.savefig(joinpath(output_folder_uw, "1a.pdf"), bbox_inches="tight")
##
##########################################################
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
###         Scatter plot fluxes CONTRASTS          ###
######################################################
mpl_scatter_density = pyimport("mpl_scatter_density")

# Scatter data preparation, such that "ice" is always the "a" part, effect for heat and momentum
var1a = fx3
var1b = fx1
var3a = fx4
var3b = fx2
sc_lbl_idx = [3, 1, 4, 2] #indexing for labels

##################
#Heat fluxes

scatter1a = var1a.wT[.!isnan.(fx1.wT) .&& .!isnan.(fx3.wT)]
scatter1b = var1b.wT[.!isnan.(fx1.wT) .&& .!isnan.(fx3.wT)]

scatter1a .*= ρ_air .* c_p
scatter1b .*= ρ_air .* c_p

scatter2b = var1b.wq[.!isnan.(fx1.wq) .&& .!isnan.(fx3.wq)]
scatter2a = var1a.wq[.!isnan.(fx1.wq) .&& .!isnan.(fx3.wq)]

scatter2a .*= L_v .* 1e-3
scatter2b .*= L_v .* 1e-3

scatter3b = var3b.wT[.!isnan.(fx2.wT) .&& .!isnan.(fx4.wT)]
scatter3a = var3a.wT[.!isnan.(fx2.wT) .&& .!isnan.(fx4.wT)]

scatter3a .*= ρ_air .* c_p
scatter3b .*= ρ_air .* c_p

ratio1 = scatter1b ./ scatter1a
ratio2 = scatter2b ./ scatter2a
ratio3 = scatter3b ./ scatter3a

##
fig = PyPlot.figure(figsize=(12, 9.5))
fig.suptitle("1a - 09.07.2025 to 10.07.2025 16:00")
gs = gridspec.GridSpec(3, 3, height_ratios=(2, 1, 1))

# --- Row 1: Scatter density ---
ax1 = fig.add_subplot(gs[1, 1], projection="scatter_density")
ax1.set_title("1m Sensible Heat Fluxes")
ax1.axhline(0, color="grey", alpha=0.6)
ax1.axvline(0, color="grey", alpha=0.6)
ax1.plot([-25,25], [-25,25], color="grey", alpha=0.3)
ax1.scatter_density(scatter1a, scatter1b, color="black", vmin=00, vmax=1000)
lim_1 = (0, 25)
ax1.set_xlim(lim_1)
ax1.set_ylim(lim_1)
#ax1.set_yticks(collect(-20:10:20))
ax1.set_xlabel("\$\\overline{w'T'}_{$(srf_type[sc_lbl_idx[1]])}~\\mathrm{[W~m^{-2}]}\$")
ax1.set_ylabel("\$\\overline{w'T'}_{$(srf_type[sc_lbl_idx[2]])}~\\mathrm{[W~m^{-2}]}\$")
ax1.grid()
ax1.set_aspect("equal")

ax2 = fig.add_subplot(gs[1, 2], projection="scatter_density")
ax2.set_title("1m Latent Heat Fluxes")
ax2.axhline(0, color="grey")
ax2.axvline(0, color="grey")
ax2.plot([-25,25], [-25,25], color="grey", alpha=0.3)
ax2.scatter_density(scatter2a, scatter2b, color="blue", vmin=0, vmax=1000)
lim_2 = (0, 25)
ax2.set_xlim(lim_2)
ax2.set_ylim(lim_2)
#ax2.set_yticks(collect(-10:5:10))
ax2.set_xlabel("\$\\overline{w'q'}_{$(srf_type[sc_lbl_idx[1]])}~\\mathrm{[W~m^{-2}]}\$")
ax2.set_ylabel("\$\\overline{w'q'}_{$(srf_type[sc_lbl_idx[2]])}~\\mathrm{[W~m^{-2}]}\$")
ax2.grid()
ax2.set_aspect("equal")

ax3 = fig.add_subplot(gs[1, 3], projection="scatter_density")
ax3.plot([-25,30], [-25,30], color="grey", alpha=0.3)
ax3.scatter_density(scatter3a, scatter3b, color="black", vmin=0, vmax=1000)
ax3.set_title("2m Sensible Heat Fluxes")
ax3.axhline(0, color="grey")
ax3.axvline(0, color="grey")
lim_3 = (0, 25)
ax3.set_xlim(lim_3)
ax3.set_ylim(lim_3)
#ax3.set_yticks(collect(-20:10:10))
ax3.set_xlabel("\$\\overline{w'T'}_{$(srf_type[sc_lbl_idx[3]])}~\\mathrm{[W~m^{-2}]}\$")
ax3.set_ylabel("\$\\overline{w'T'}_{$(srf_type[sc_lbl_idx[4]])}~\\mathrm{[W~m^{-2}]}\$")
ax3.grid()
ax3.set_aspect("equal")

# --- Row 2: Overlaid histograms of individual fluxes ---
ax4 = fig.add_subplot(gs[2, 1])
ax4.hist(scatter1a, bins=collect(0:0.1:25), color="grey", alpha=0.5, density=true, label=instr_labels[sc_lbl_idx[1]])
ax4.hist(scatter1b, bins=collect(0:0.1:25), color="orange", alpha=0.5, density=true, label=instr_labels[sc_lbl_idx[2]])
ax4.set_xlabel(L"\overline{w'T'}~\mathrm{[W~m^{-2}]}")
ax4.set_ylabel("PDF")
ax4.set_xlim(ax1.get_xlim())
ax4.legend()
ax4.grid(alpha=0.3)

ax5 = fig.add_subplot(gs[2, 2])
ax5.hist(scatter2a, bins=collect(0:0.1:25), color="grey", alpha=0.5, density=true, label=instr_labels[sc_lbl_idx[1]])
ax5.hist(scatter2b, bins=collect(0:0.1:25), color="orange", alpha=0.5, density=true, label=instr_labels[sc_lbl_idx[2]])
ax5.set_xlabel(L"\overline{w'q'}~\mathrm{[W~m^{-2}]}")
ax5.set_xlim(ax2.get_xlim())
ax5.legend()
ax5.grid(alpha=0.3)

ax6 = fig.add_subplot(gs[2, 3])
ax6.hist(scatter3a, bins=collect(0:0.1:25), color="grey", alpha=0.5, density=true, label=instr_labels[sc_lbl_idx[3]])
ax6.hist(scatter3b, bins=collect(0:0.1:25), color="orange", alpha=0.5, density=true, label=instr_labels[sc_lbl_idx[4]])
ax6.set_xlabel(L"\overline{w'T'}~\mathrm{[W~m^{-2}]}")
ax6.set_xlim(ax3.get_xlim())
ax6.legend()
ax6.grid(alpha=0.3)

# --- Row 3: Histograms of ratios ---
ax7 = fig.add_subplot(gs[3, 1])
ax7.hist(ratio1, bins=collect(0:0.01:2.0), color="black", alpha=0.7, density=true)
ax7.axvline(1, color="red", linestyle="--", alpha=0.3, label="1:1")
ax7.set_xlabel("\$\\left(\\overline{w'T'}_{$(srf_type[sc_lbl_idx[2]])} ~/~\\overline{w'T'}_{$(srf_type[sc_lbl_idx[1]])}\\right)_{IRG}\$")
ax7.set_ylabel("PDF")
ax7.grid(alpha=0.3)

ax8 = fig.add_subplot(gs[3, 2])
ax8.hist(ratio2, bins=collect(-1.0:0.01:2.0), color="blue", alpha=0.7, density=true)
ax8.axvline(1, color="red", linestyle="--", alpha=0.3, label="1:1")
ax8.set_xlabel("\$\\left(\\overline{w'q'}_{$(srf_type[sc_lbl_idx[2]])} ~/~ \\overline{w'q'}_{$(srf_type[sc_lbl_idx[1]])}\\right)_{CSAT}\$")
ax8.grid(alpha=0.3)

ax9 = fig.add_subplot(gs[3, 3])
ax9.hist(ratio3, bins=collect(-0.5:0.01:1.5), color="black", alpha=0.7, density=true)
ax9.axvline(1, color="red", linestyle="--", alpha=0.3, label="1:1")
ax9.set_xlabel("\$\\left(\\overline{w'T'}_{$(srf_type[sc_lbl_idx[4]])} ~/~\\overline{w'T'}_{$(srf_type[sc_lbl_idx[3]])}\\right)_{CSAT}\$")
ax9.grid(alpha=0.3)

PyPlot.tight_layout()
##
output_folder = "/home/haugened/Documents/data/CONTRASTS/plots/correlation_heat/"
#PyPlot.savefig(joinpath(output_folder, "1a.pdf"), bbox_inches="tight")

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
##
##############################
#Scatter plot momentum fluxes CONTRASTS
# Momentum flux (u_star) scatter data preparation
# Columns 1-2: Use var1a/var1b and var3a/var3b from heat flux setup (tower comparisons)
scatter_mom1a = var1a.u_star[.!isnan.(var1a.u_star) .&& .!isnan.(var1b.u_star)]
scatter_mom1b = var1b.u_star[.!isnan.(var1a.u_star) .&& .!isnan.(var1b.u_star)]

scatter_mom2a = var3a.u_star[.!isnan.(var3a.u_star) .&& .!isnan.(var3b.u_star)]
scatter_mom2b = var3b.u_star[.!isnan.(var3a.u_star) .&& .!isnan.(var3b.u_star)]

# Columns 3-4: Direct comparison of heights at same location
scatter_mom3a = fx1.u_star[.!isnan.(fx1.u_star) .&& .!isnan.(fx2.u_star)]
scatter_mom3b = fx2.u_star[.!isnan.(fx1.u_star) .&& .!isnan.(fx2.u_star)]

scatter_mom4a = fx3.u_star[.!isnan.(fx3.u_star) .&& .!isnan.(fx4.u_star)]
scatter_mom4b = fx4.u_star[.!isnan.(fx3.u_star) .&& .!isnan.(fx4.u_star)]

# Calculate ratios (y/x for each column)
ratio_mom1 = scatter_mom1b ./ scatter_mom1a
ratio_mom2 = scatter_mom2b ./ scatter_mom2a
ratio_mom3 = scatter_mom3b ./ scatter_mom3a
ratio_mom4 = scatter_mom4b ./ scatter_mom4a

fig_mom = PyPlot.figure(figsize=(16, 9.5))
fig_mom.suptitle("1b-2 - Friction Velocity")
gs_mom = gridspec.GridSpec(3, 4, height_ratios=(2, 1, 1))

# --- Row 1: Scatter density plots ---
# Column 1: 1m comparison (tower comparison using var1a/var1b)
ax_m1 = fig_mom.add_subplot(gs_mom[1, 1], projection="scatter_density")
ax_m1.set_title("1m")
ax_m1.plot([0,1], [0,1], color="grey", alpha=0.3)
ax_m1.scatter_density(scatter_mom1a, scatter_mom1b, color="red", vmin=0, vmax=1000)
ax_m1.set_xlabel("\$u_{*,$(srf_type[sc_lbl_idx[1]])}~\\mathrm{[m~s^{-1}]}\$")
ax_m1.set_ylabel("\$u_{*,$(srf_type[sc_lbl_idx[2]])}~\\mathrm{[m~s^{-1}]}\$")
ax_m1.grid()
lim_m1 = (0.3, 0.8)
ax_m1.set_xlim(lim_m1)
ax_m1.set_ylim(lim_m1)
ax_m1.set_aspect("equal")

# Column 2: 2m comparison (tower comparison using var3a/var3b)
ax_m2 = fig_mom.add_subplot(gs_mom[1, 2], projection="scatter_density")
ax_m2.set_title("2m")
ax_m2.plot([0,1], [0,1], color="grey", alpha=0.3)
ax_m2.scatter_density(scatter_mom2a, scatter_mom2b, color="red", vmin=0, vmax=1000)
ax_m2.set_xlabel("\$u_{*,$(srf_type[sc_lbl_idx[3]])}~\\mathrm{[m~s^{-1}]}\$")
ax_m2.set_ylabel("\$u_{*,$(srf_type[sc_lbl_idx[4]])}~\\mathrm{[m~s^{-1}]}\$")
ax_m2.grid()
lim_m2 = (0.4, 0.9)
ax_m2.set_xlim(lim_m2)
ax_m2.set_ylim(lim_m2)
ax_m2.set_aspect("equal")

# Column 3: Height comparison at tower 1 (fx1 vs fx2)
ax_m3 = fig_mom.add_subplot(gs_mom[1, 3], projection="scatter_density")
ax_m3.set_title("$(srf_type[1]) (tower 1)")
ax_m3.plot([0,1], [0,1], color="grey", alpha=0.3)
ax_m3.scatter_density(scatter_mom3a, scatter_mom3b, color="red", vmin=0, vmax=1000)
ax_m3.set_xlabel("\$u_{*,$(heights[1])m}~\\mathrm{[m~s^{-1}]}\$")
ax_m3.set_ylabel("\$u_{*,$(heights[2])m}~\\mathrm{[m~s^{-1}]}\$")
ax_m3.grid()
lim_m3 = (0.3, 0.9)
ax_m3.set_xlim(lim_m3)
ax_m3.set_ylim(lim_m3)
ax_m3.set_aspect("equal")

# Column 4: Height comparison at tower 2 (fx3 vs fx4)
ax_m4 = fig_mom.add_subplot(gs_mom[1, 4], projection="scatter_density")
ax_m4.set_title("$(srf_type[3]) (tower 2)")
ax_m4.plot([0,1], [0,1], color="grey", alpha=0.3)
ax_m4.scatter_density(scatter_mom4a, scatter_mom4b, color="red", vmin=0, vmax=1000)
ax_m4.set_xlabel("\$u_{*,$(heights[3])m}~\\mathrm{[m~s^{-1}]}\$")
ax_m4.set_ylabel("\$u_{*,$(heights[4])m}~\\mathrm{[m~s^{-1}]}\$")
ax_m4.grid()
lim_m4 = (0.3, 0.9)
ax_m4.set_xlim(lim_m4)
ax_m4.set_ylim(lim_m4)
ax_m4.set_aspect("equal")

# --- Row 2: Overlaid histograms of individual u_star ---
ax_m5 = fig_mom.add_subplot(gs_mom[2, 1])
ax_m5.hist(scatter_mom1a, bins=50, color="grey", alpha=0.5, density=true, label=instr_labels[sc_lbl_idx[1]])
ax_m5.hist(scatter_mom1b, bins=50, color="orange", alpha=0.5, density=true, label=instr_labels[sc_lbl_idx[2]])
ax_m5.set_xlabel(L"u_*~\mathrm{[m~s^{-1}]}")
ax_m5.set_ylabel("PDF")
ax_m5.set_xlim(ax_m1.get_xlim())
ax_m5.legend()
ax_m5.grid(alpha=0.3)

ax_m6 = fig_mom.add_subplot(gs_mom[2, 2])
ax_m6.hist(scatter_mom2a, bins=50, color="grey", alpha=0.5, density=true, label=instr_labels[sc_lbl_idx[3]])
ax_m6.hist(scatter_mom2b, bins=50, color="orange", alpha=0.5, density=true, label=instr_labels[sc_lbl_idx[4]])
ax_m6.set_xlabel(L"u_*~\mathrm{[m~s^{-1}]}")
ax_m6.set_xlim(ax_m2.get_xlim())
ax_m6.legend()
ax_m6.grid(alpha=0.3)

ax_m7 = fig_mom.add_subplot(gs_mom[2, 3])
ax_m7.hist(scatter_mom3a, bins=50, color="grey", alpha=0.5, density=true, label=instr_labels[1])
ax_m7.hist(scatter_mom3b, bins=50, color="orange", alpha=0.5, density=true, label=instr_labels[2])
ax_m7.set_xlabel(L"u_*~\mathrm{[m~s^{-1}]}")
ax_m7.set_xlim(ax_m3.get_xlim())
ax_m7.legend()
ax_m7.grid(alpha=0.3)

ax_m8 = fig_mom.add_subplot(gs_mom[2, 4])
ax_m8.hist(scatter_mom4a, bins=50, color="grey", alpha=0.5, density=true, label=instr_labels[3])
ax_m8.hist(scatter_mom4b, bins=50, color="orange", alpha=0.5, density=true, label=instr_labels[4])
ax_m8.set_xlabel(L"u_*~\mathrm{[m~s^{-1}]}")
ax_m8.set_xlim(ax_m4.get_xlim())
ax_m8.legend()
ax_m8.grid(alpha=0.3)

# --- Row 3: Histograms of ratios ---
ax_m9 = fig_mom.add_subplot(gs_mom[3, 1])
ax_m9.hist(ratio_mom1, bins=collect(0.9:0.001:1.1), color="red", alpha=0.7, density=true)
ax_m9.axvline(1, color="black", linestyle="--", alpha=0.3, label="1:1")
ax_m9.set_xlabel("\$\\left(u_{*,$(srf_type[sc_lbl_idx[2]])} ~/~ u_{*,$(srf_type[sc_lbl_idx[1]])}\\right)_{1m}\$")
ax_m9.set_ylabel("PDF")
ax_m9.grid(alpha=0.3)

ax_m10 = fig_mom.add_subplot(gs_mom[3, 2])
ax_m10.hist(ratio_mom2, bins=collect(0.9:0.001:1.1), color="red", alpha=0.7, density=true)
ax_m10.axvline(1, color="black", linestyle="--", alpha=0.3, label="1:1")
ax_m10.set_xlabel("\$\\left(u_{*,$(srf_type[sc_lbl_idx[4]])} ~/~ u_{*,$(srf_type[sc_lbl_idx[3]])}\\right)_{2m}\$")
ax_m10.grid(alpha=0.3)

ax_m11 = fig_mom.add_subplot(gs_mom[3, 3])
ax_m11.hist(ratio_mom3, bins=collect(0.9:0.001:1.3), color="red", alpha=0.7, density=true)
ax_m11.axvline(1, color="black", linestyle="--", alpha=0.3, label="1:1")
ax_m11.set_xlabel("\$\\left(u_{*,$(heights[2])m} ~/~ u_{*,$(heights[1])m}\\right)_{$(srf_type[1])}\$")
ax_m11.grid(alpha=0.3)

ax_m12 = fig_mom.add_subplot(gs_mom[3, 4])
ax_m12.hist(ratio_mom4, bins=collect(0.9:0.001:1.2), color="red", alpha=0.7, density=true)
ax_m12.axvline(1, color="black", linestyle="--", alpha=0.3, label="1:1")
ax_m12.set_xlabel("\$\\left(u_{*,$(heights[4])m} ~/~ u_{*,$(heights[3])m}\\right)_{$(srf_type[3])}\$")
ax_m12.grid(alpha=0.3)

PyPlot.tight_layout()
##
output_folder_mom = "/home/haugened/Documents/data/CONTRASTS/plots/correlation_momentum/"
#PyPlot.savefig(joinpath(output_folder_mom, "1b_2.pdf"), bbox_inches="tight")

#calculating correlations for momentum
cor_mom1 = nancor(var1a.u_star, var1b.u_star)
cor_mom2 = nancor(var3a.u_star, var3b.u_star)
cor_mom3 = nancor(fx1.u_star, fx2.u_star)
cor_mom4 = nancor(fx3.u_star, fx4.u_star)

#calculating mean bias for momentum
scatter_mom1a_mean = nanmean(scatter_mom1a)
scatter_mom1b_mean = nanmean(scatter_mom1b)
scatter_mom2a_mean = nanmean(scatter_mom2a)
scatter_mom2b_mean = nanmean(scatter_mom2b)
scatter_mom3a_mean = nanmean(scatter_mom3a)
scatter_mom3b_mean = nanmean(scatter_mom3b)
scatter_mom4a_mean = nanmean(scatter_mom4a)
scatter_mom4b_mean = nanmean(scatter_mom4b)

mad_mom1 = scatter_mom1b_mean - scatter_mom1a_mean
mad_mom2 = scatter_mom2b_mean - scatter_mom2a_mean
mad_mom3 = scatter_mom3b_mean - scatter_mom3a_mean
mad_mom4 = scatter_mom4b_mean - scatter_mom4a_mean
##
######################################################