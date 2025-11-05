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
datapath = "/home/haugened/Documents/data/"
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
turb.missing2nan!(evaldf5)
turb.missing2nan!(evaldf6)

#Reynolds averaging times (determine usind MRD)
ra1 = Millisecond(2^11*timestep) #Second(330)
ra2 = Millisecond(2^11*timestep) #Second(330)
ra3 = Millisecond(2^11*timestep) #Second(320)
ra4 = Millisecond(2^11*timestep) #Second(210)
ra5 = Millisecond(2^9*timestep) #Second(300)
ra6 = Millisecond(2^10*timestep) #Second(55)

fx1_raw = turb.turbflux(evaldf1, ra1)
fx2_raw = turb.turbflux(evaldf2, ra2)
fx3_raw = turb.turbflux(evaldf3, ra3)
fx4_raw = turb.turbflux(evaldf4, ra4)
fx5_raw = turb.turbflux(evaldf5, ra5)
fx6_raw = turb.turbflux(evaldf6, ra6)

#averaging
fx1 = turb.avgflux(fx1_raw, Second(600))
fx2 = turb.avgflux(fx2_raw, Second(600))
fx3 = turb.avgflux(fx3_raw, Second(600))
fx4 = turb.avgflux(fx4_raw, Second(600))
fx5 = turb.avgflux(fx5_raw, Second(600))
fx6 = turb.avgflux(fx6_raw, Second(600))

#Obukhov-length
L1 = turb.obukhov(fx1_raw, Minute(30))
L2 = turb.obukhov(fx2_raw, Minute(10))
L3 = turb.obukhov(fx3_raw, Minute(10))
L4 = turb.obukhov(fx4_raw, Minute(10))
try
    L5 = turb.obukhov(fx5_raw, Minute(30))
catch e
end
L6 = turb.obukhov(fx6_raw, Minute(10))

######################################
#=
#save to NetCDF4
outpath = "/home/haugened/Documents/data/fluxes"

#condense to 1min data
fx1_exp = fx1[1:20*60:end,:]
fx2_exp = fx2[1:20*60:end,:]
fx3_exp = fx3[1:20*60:end,:]
fx4_exp = fx4[1:20*60:end,:]
fx5_exp = fx5[1:20*60:end,:]
fx6_exp = fx6[1:20*60:end,:]
turb.saveturbasnetcdf(fx1_exp, joinpath(outpath, "t1irg_fx.nc"))
turb.saveturbasnetcdf(fx2_exp, joinpath(outpath, "t2irg_fx.nc"))
turb.saveturbasnetcdf(fx3_exp, joinpath(outpath, "t2lcsat_fx.nc"))
turb.saveturbasnetcdf(fx4_exp, joinpath(outpath, "t2ucsat_fx.nc"))
turb.saveturbasnetcdf(fx5_exp, joinpath(outpath, "kaijo_fx.nc"))
turb.saveturbasnetcdf(fx6_exp, joinpath(outpath, "tjk_fx.nc"))
=#

######################################
#calculate OSHD-'model'-fluxes
#=
tjksonicmean = gen.movingaverage(sqrt.(evaldf6.u .^ 2 .+ evaldf6.v .^ 2), 20 * 600)

timestotake = fill(false, length(evaldf6.time))
for ix in axes(tjkmeteodata.time, 1)
    tmp = abs.(evaldf6.time .- tjkmeteodata.time[ix])
    minim = minimum(tmp)
    apos = findfirst(x -> x == minim, tmp)
    timestotake[apos] = true
end

tjksonicmeantotake = round.(tjksonicmean[timestotake]; digits=2)

(h_oshd, RiB, dT) = turb.OSHD_SHF(tjkmeteodata.tair_CS215_avg .+ 273.15, fill(273.15, length(tjksonicmeantotake)), tjksonicmeantotake, 0.005, 3.7, 4.4)
(h_oshd_lwrout, RiB_LW, dT_LW) = turb.OSHD_SHF(tjkmeteodata.tair_CS215_avg .+ 273.15, (tjkmeteodata.LWR_out ./ (5.67e-8)) .^ (1 / 4), tjksonicmeantotake, 0.005, 3.7, 4.4)

#plot OSHD flux comparison
fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.plot(tjkmeteodata.time, h_oshd, label="FSM (273K)")#T_{srf}=273.15 K")
#ax.plot(tjkmeteodata.time, h_oshd_lwrout, label="FSM (LWR)") #T_{srf}=(LWR_{out}/σ)^(1/4)")
ax.plot(fx6.time[1:20*600:end], fx6.wT[1:20*600:end] .* (c_p * ρ_air), label="TJK")
ax.plot(fx2.time[1:20*600:end], fx2.wT[1:20*600:end] .* (c_p * ρ_air), label="T2IRG")
#ax.plot(fx3.time[1:20*600:end], fx3.wT[1:20*600:end].*(c_p*ρ_air), label="T2LCSAT")
#ax.plot(fx4.time[1:20*600:end], fx4.wT[1:20*600:end].*(c_p*ρ_air), label="T2UCSAT")
#ax.plot(fx5.time[1:20*600:end], fx5.wT[1:20*600:end].*(c_p*ρ_air), label="KAIJO")
ax.set_ylabel(L"H~\mathrm{[W~m^{-2}]}")
ax.grid()
ax.legend()

#=
#create additional data for day and night
daystart = Time(08, 00, 00)
dayend = Time(18, 00, 00)
(fx1day, fx1night) = turb.splitdaynight(fx1, daystart, dayend)
(fx2day, fx2night) = turb.splitdaynight(fx2, daystart, dayend)
(fx3day, fx3night) = turb.splitdaynight(fx3, daystart, dayend)
(fx4day, fx4night) = turb.splitdaynight(fx4, daystart, dayend)
(evaldf1day, evaldf1night) = turb.splitdaynight(evaldf1, daystart, dayend)
(evaldf2day, evaldf2night) = turb.splitdaynight(evaldf2, daystart, dayend)
(evaldf3day, evaldf3night) = turb.splitdaynight(evaldf3, daystart, dayend)
(evaldf4day, evaldf4night) = turb.splitdaynight(evaldf4, daystart, dayend)
=#
##
=#

######################################################
###                    PLOTS                       ###
######################################################
##
#flux time series
fig = PyPlot.figure()
ax = fig.add_subplot(111)
#ax.set_title("Turbulent fluxes")
#pt1 = ax.plot(fx1.time[1:20*600:end], fx1.wT[1:20*600:end], label="T1IRG")
pt2 = ax.plot(fx2.time[1:20*60:end], fx2.wT[1:20*60:end], label="T2IRG")
pt3 = ax.plot(fx3.time[1:20*60:end], fx3.wT[1:20*60:end], label="T2LCSAT")
pt4 = ax.plot(fx4.time[1:20*60:end], fx4.wT[1:20*60:end], label="T2UCSAT")
pt5 = ax.plot(fx5.time[1:20*60:end], fx5.wT[1:20*60:end], label="KAIJO")
#pt6 = ax.plot(fx6.time[1:20*60:end], fx6.wT[1:20*60:end], label="TJK")
#pt2wq = ax.plot(fx2.time[1:20*60:end], fx2.wq[1:20*60:end], label="w'q'(T2IRG)")
#pt1wq = ax.plot(fx1.time[1:20*600:end], fx1.wq[1:20*600:end], label="w'q'(T1IRG)")
#ax.set_ylabel(L"u_\ast~\mathrm{[m~s^{-1}]}")
#ax.set_ylabel(L"\overline{u'w'}~\mathrm{[m^2~s^{-2}]}")
#ax.set_ylabel(L"\overline{w'q'}~\mathrm{[m~g~s^{-1}~m^{-3}]}")
#ax.set_xlabel("31.05.2021")
ax.set_ylabel(L"\overline{w'T_S'}~\mathrm{[K~m~s^{-1}]}")
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
##
######################################################
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
######################################################
#compare fluxes with modeled in components
##
step = 20 * 600
cmap = PyPlot.get_cmap("tab10")
fig = PyPlot.figure()
gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 1])
ax3 = fig.add_subplot(gs[3, 1])
ax3.set_ylabel(L"Ri_B")
ax3.plot(tjkmeteodata.time, RiB_LW, color=cmap(0))
#minorlocator = pydates.MinuteLocator(interval=1)
#ax5.xaxis.set_minor_locator(minorlocator)
ax3.grid()
ax1 = fig.add_subplot(gs[1, 1], sharex=ax3)
ax1.set_ylabel(L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}")
ax1.plot(tjkmeteodata.time, h_oshd_lwrout, color=cmap(0), label=L"FSM")
ax1.plot(fx6.time[1:step:end], fx6.wT[1:step:end] .* (ρ_air * c_p), color=cmap(1), label="TJK")
ax1.plot(fx2.time[1:step:end], fx2.wT[1:step:end] .* (ρ_air * c_p), color=cmap(2), label="T2IRG")
#ax1.plot(fx5.time[1:step:end], fx5.wT[1:step:end].*(ρ_air*c_p), color=cmap(3), label="KAIJO")
ax1.legend()
ax1.grid()
ax1.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(gs[2, 1], sharex=ax3)
ax2.set_ylabel(L"\Delta T = T_S - T_A~\mathrm{[K]}")
ax2.plot(tjkmeteodata.time, dT_LW, color=cmap(0))
ax2.grid()
ax2.tick_params(axis="x", labelbottom=false)
PyPlot.show()
##
######################################################
#4-window fluxes
fig = PyPlot.figure()
gs = gridspec.GridSpec(4, 1)
axl = fig.add_subplot(gs[4, 1])
axl.grid()
axl.plot(fx2.time[1:1200:end], fx2.wT[1:1200:end])
axl.set_ylabel(L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}")
ax1 = fig.add_subplot(gs[1, 1], sharex=axl)
ax1.grid()
ax1.plot(fx6.time[1:1200:end], fx6.wT[1:1200:end])
ax1.set_ylabel(L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}")
ax1.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(gs[2, 1], sharex=axl)
ax2.grid()
ax2.plot(fx4.time[1:1200:end], fx4.wT[1:1200:end])
ax2.set_ylabel(L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}")
ax2.tick_params(axis="x", labelbottom=false)
ax3 = fig.add_subplot(gs[3, 1], sharex=axl)
ax3.grid()
ax3.plot(fx3.time[1:1200:end], fx3.wT[1:1200:end])
ax3.set_ylabel(L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}")
ax3.tick_params(axis="x", labelbottom=false)


######################################################
#fluxes at different scales

scale1 = Second(600)
scale2 = Second(10)

fx1_s1 = turb.avgflux(turb.turbflux(evaldf1, scale1), scale1)
fx2_s1 = turb.avgflux(turb.turbflux(evaldf2, scale1), scale1)
fx3_s1 = turb.avgflux(turb.turbflux(evaldf3, scale1), scale1)
fx4_s1 = turb.avgflux(turb.turbflux(evaldf4, scale1), scale1)
fx5_s1 = turb.avgflux(turb.turbflux(evaldf5, scale1), scale1)
fx6_s1 = turb.avgflux(turb.turbflux(evaldf6, scale1), scale1)

fx1_s2 = turb.avgflux(turb.turbflux(evaldf1, scale2), scale2)
fx2_s2 = turb.avgflux(turb.turbflux(evaldf2, scale2), scale2)
fx3_s2 = turb.avgflux(turb.turbflux(evaldf3, scale2), scale2)
fx4_s2 = turb.avgflux(turb.turbflux(evaldf4, scale2), scale2)
fx5_s2 = turb.avgflux(turb.turbflux(evaldf5, scale2), scale2)
fx6_s2 = turb.avgflux(turb.turbflux(evaldf6, scale2), scale2)
##
fig = PyPlot.figure()
gs = gridspec.GridSpec(2, 1)
ax1 = fig.add_subplot(gs[1, 1])
ax1.set_ylabel(L"\overline{w'T'}_{10min-10s}~\mathrm{[K~m~s^{-1}]}")
ax1.plot(fx1_s1.time[1:10:end], fx1_s1.wT[1:10:end] .- gen.movingaverage(fx1_s2.wT, round(Int, 20 * scale1 / scale2))[1:10:end], label="T1IRG")
ax1.plot(fx2_s1.time[1:10:end], fx2_s1.wT[1:10:end] .- gen.movingaverage(fx2_s2.wT, round(Int, 20 * scale1 / scale2))[1:10:end], label="T2IRG")
ax1.plot(fx3_s1.time[1:10:end], fx3_s1.wT[1:10:end] .- gen.movingaverage(fx3_s2.wT, round(Int, 20 * scale1 / scale2))[1:10:end], label="T2UCSAT")
ax1.plot(fx4_s1.time[1:10:end], fx4_s1.wT[1:10:end] .- gen.movingaverage(fx4_s2.wT, round(Int, 20 * scale1 / scale2))[1:10:end], label="T2LCSAT")
ax1.plot(fx5_s1.time[1:10:end], fx5_s1.wT[1:10:end] .- gen.movingaverage(fx5_s2.wT, round(Int, 20 * scale1 / scale2))[1:10:end], label="Kaijo")
ax1.plot(fx6_s1.time[1:10:end], fx6_s1.wT[1:10:end] .- gen.movingaverage(fx6_s2.wT, round(Int, 20 * scale1 / scale2))[1:10:end], label="TJK")
ax1.grid()
ax1.legend()
ax2 = fig.add_subplot(gs[2, 1], sharex=ax1)
ax2.set_ylabel(L"\overline{w'T'}_{10s}~\mathrm{[K~m~s^{-1}]}")
ax2.plot(fx1_s2.time[1:10:end], fx1_s2.wT[1:10:end], label="T1IRG")
ax2.plot(fx2_s2.time[1:10:end], fx2_s2.wT[1:10:end], label="T2IRG")
ax2.plot(fx3_s2.time[1:10:end], fx3_s2.wT[1:10:end], label="T2LCSAT")
ax2.plot(fx4_s2.time[1:10:end], fx4_s2.wT[1:10:end], label="T2UCSAT")
ax2.plot(fx5_s2.time[1:10:end], fx5_s2.wT[1:10:end], label="Kaijo")
ax2.plot(fx6_s2.time[1:10:end], fx6_s2.wT[1:10:end], label="TJK")
ax2.grid()
#ax2.legend()
##

######################################################
#comparison u'w' and v'w'
plotstep = 120*20
fig = PyPlot.figure()
gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
fig.suptitle("Comparison u'w' and v'w' - TJK")
ax1 = fig.add_subplot(gs[1, 1])
ax1.set_ylabel(L"\overline{u'w'}~[\mathrm{m^2~s^{-2}}]")
ax1.plot(fx6.time[1:plotstep:end], fx6.uw[1:plotstep:end])
ax1.tick_params(axis="x", labelbottom=false)
ax1.grid()
ax1.legend()
ax2 = fig.add_subplot(gs[2, 1], sharex=ax1)
ax2.plot(fx6.time[1:plotstep:end], fx6.vw[1:plotstep:end])
ax2.set_ylabel(L"\overline{v'w'}~[\mathrm{m^2~s^{-2}}]")
#date_format = pydates.DateFormatter("%H:%M")
#ax2.xaxis.set_major_formatter(date_format)
ax2.grid()
ax3 = fig.add_subplot(gs[3,1], sharex=ax1)
ax3.set_ylabel(L"\overline{v'w'}/\overline{u'w'}")
ax3.plot(fx6.time[1:plotstep:end], fx6.vw[1:plotstep:end]./fx6.uw[1:plotstep:end])
ax3.grid()
ax3.set_ylim(-1,1)
##

######################################################
#plot special u'w' and v'w' with DR per period
fx_dr_raw = turb.turbfluxdrperperiod(evaldf6, Second(100), Minute(30))
fx_dr = turb.avgflux(fx_dr_raw, Minute(30))
plotstep = 120*20
fig = PyPlot.figure()
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
fig.suptitle("Comparison u'w' and v'w' - TJK")
ax1 = fig.add_subplot(gs[1, 1])
ax1.set_ylabel(L"\overline{u'w'}~[\mathrm{m^2~s^{-2}}]")
ax1.plot(fx_dr.time[1:plotstep:end], fx_dr.uw[1:plotstep:end])
ax1.tick_params(axis="x", labelbottom=false)
ax1.grid()
ax1.set_ylim(-0.1,0.05)
#ax1.legend()
ax2 = fig.add_subplot(gs[2, 1], sharex=ax1)
ax2.plot(fx_dr.time[1:plotstep:end], fx_dr.vw[1:plotstep:end])
ax2.set_ylabel(L"\overline{v'w'}~[\mathrm{m^2~s^{-2}}]")
#date_format = pydates.DateFormatter("%H:%M")
#ax2.xaxis.set_major_formatter(date_format)
ax2.grid()
ax2.set_ylim(-0.1,0.05)

######################################################
#plot special u'w'-profile with DR for each period
fx2_dr_raw = turb.turbfluxdrperperiod(evaldf2, Second(100), Minute(30))
fx2_dr = turb.avgflux(fx2_dr_raw, Minute(30))
fx3_dr_raw = turb.turbfluxdrperperiod(evaldf3, Second(100), Minute(30))
fx3_dr = turb.avgflux(fx3_dr_raw, Minute(30))
fx4_dr_raw = turb.turbfluxdrperperiod(evaldf4, Second(100), Minute(30))
fx4_dr = turb.avgflux(fx4_dr_raw, Minute(30))
fx6_dr_raw = turb.turbfluxdrperperiod(evaldf6, Second(100), Minute(30))
fx6_dr = turb.avgflux(fx6_dr_raw, Minute(30))

plotstep = 120*20
fig = PyPlot.figure()
gs = gridspec.GridSpec(1, 1)
fig.suptitle("u'w' - vertical profile")
ax1 = fig.add_subplot(gs[1, 1])
ax1.set_ylabel(L"\overline{u'w'}~[\mathrm{m^2~s^{-2}}]")
ax1.plot(fx2_dr.time[1:plotstep:end], fx2_dr.uw[1:plotstep:end], label="T2IRG", alpha=0.5)
ax1.plot(fx3_dr.time[1:plotstep:end], fx3_dr.uw[1:plotstep:end], label="T2LCSAT", alpha=0.5)
ax1.plot(fx4_dr.time[1:plotstep:end], fx4_dr.uw[1:plotstep:end], label="T2UCSAT", alpha=0.5)
ax1.grid()
ax1.set_ylim(-0.1,0.05)
ax1.legend()

fig = PyPlot.figure()
gs = gridspec.GridSpec(1, 1)
fig.suptitle("u'T' - vertical profile")
ax1 = fig.add_subplot(gs[1, 1])
ax1.set_ylabel(L"\overline{u'T'}~[\mathrm{K~m~s^{-1}}]")
ax1.plot(fx2_dr.time[1:plotstep:end], fx2_dr.uT[1:plotstep:end], label="T2IRG", alpha=0.5)
#ax1.plot(fx3_dr.time[1:plotstep:end], fx3_dr.uw[1:plotstep:end], label="T2LCSAT", alpha=0.5)
ax1.plot(fx4_dr.time[1:plotstep:end], fx4_dr.uT[1:plotstep:end], label="T2UCSAT", alpha=0.5)
ax1.grid()
ax1.set_ylim(-0.1,0.05)
ax1.legend()

t2ucsat_abovejet = zeros(Bool, size(fx4_dr, 1))
t2irg_abovejet = zeros(Bool, size(fx2_dr, 1))
t2ucsat_belowjet = zeros(Bool, size(fx4_dr, 1))
t2irg_belowjet = zeros(Bool, size(fx2_dr, 1))

t2ucsat_abovejet = fx4_dr.uw .> 0
t2irg_abovejet = fx2_dr.uw .> 0
t2ucsat_belowjet = fx4_dr.uw .< 0
t2irg_belowjet = fx2_dr.uw .< 0

t2ucsat_uT_abovejet = fx4_dr.uT .< 0
t2irg_uT_abovejet = fx2_dr.uT .< 0
t2ucsat_uT_belowjet = fx4_dr.uT .> 0
t2irg_uT_belowjet = fx2_dr.uT .> 0

jetbetween_t2ucsat_t2irg = t2ucsat_abovejet .&& t2irg_belowjet
jetbetween_t2ucsat_t2irg_uT = t2ucsat_uT_abovejet .&& t2irg_uT_belowjet

category = zeros(Int64, length(t2ucsat_uT_abovejet))
category[t2ucsat_belowjet] .= 1
category[jetbetween_t2ucsat_t2irg] .= 2
category[t2irg_abovejet] .= 3

#create winddirection classes
windclass = fill(NaN, length(tjkmeteodata.wind_mean_vector_direction))
windclass[310 .<= tjkmeteodata.wind_mean_vector_direction.||tjkmeteodata.wind_mean_vector_direction.<=20] .= 1
windclass[130 .<= tjkmeteodata.wind_mean_vector_direction .<= 200] .= 2

fig = PyPlot.figure()
gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
ax = fig.add_subplot(gs[1,1])
#ax.plot(fx4_dr.time[1:plotstep:end], jetbetween_t2ucsat_t2irg[1:plotstep:end], ".", markersize=0.001)
ax.plot(fx2_dr.time[1:plotstep:end], fx2_dr.uw[1:plotstep:end], label="T2IRG", alpha=0.7)
ax.plot(fx3_dr.time[1:plotstep:end], fx3_dr.uw[1:plotstep:end], label="T2LCSAT", alpha=0.7)
ax.plot(fx4_dr.time[1:plotstep:end], fx4_dr.uw[1:plotstep:end], label="T2UCSAT", alpha=0.7)
ax.plot(fx6_dr.time[1:plotstep:end], fx6_dr.uw[1:plotstep:end], label="TJK", alpha=0.5)
#ax.plot(fx4_dr.time[1:plotstep:end], jetbetween_t2ucsat_t2irg_uT[1:plotstep:end].-0.02, ".")
#ax.pcolormesh(pydates.date2num.(fx4_dr.time[1:plotstep:end]), collect(range(0.95, 1.05, 2)), repeat(transpose(category[1:plotstep:end]), 2), shading="auto", cmap="Dark2", alpha=1, vmin=0, vmax=9)
#ax.pcolormesh(pydates.date2num.(fx4_dr.time[1:plotstep:end]), collect(range(0.9, 0.95, 2)), repeat(transpose(jetbetween_t2ucsat_t2irg_uT[1:plotstep:end]), 2), shading="auto", cmap="Set1", alpha=0.5, vmin=0, vmax=9)
#ax.set_ylim(0.95, 1.05)
ax.grid()
ax.legend()
#ax.tick_params(axis="y", labelleft=false)
ax.tick_params(axis="x", labelbottom=false)
ax.set_ylabel(L"\overline{u'w'}~\mathrm{[m^2~s^{-2}]}")
#ax.set_ylabel("LLJ height")
ax2 = fig.add_subplot(gs[2,1], sharex=ax)
ax2.plot(pydates.date2num.(tjkmeteodata.time), tjkmeteodata.wind_mean_scalar, color="black")
ax2.pcolormesh(pydates.date2num.(tjkmeteodata.time), collect(range(-10, 10, 2)), repeat(transpose(windclass), 2), shading="auto", cmap="Set1", vmin=1, vmax=9, alpha=0.5)
#ax2.tick_params(axis="y", labelleft=false)
#ax2.set_ylim(0.95,1.05)
ax2.grid()
ax2.set_ylim(0,10)
ax2.set_ylabel(L"u~\mathrm{[m~s^{-1}]}")

######################################################
#comparison latent and sensible heat fluxes
PyPlot.pygui(true)
fig = PyPlot.figure()
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
fig.suptitle("Comparison latent - sensible heat fluxes")
ax1 = fig.add_subplot(gs[1, 1])
ax1.set_ylabel(L"\overline{w'\theta'}~[\mathrm{W~m^{-2}}]")
#ax1.plot(fx1.time, fx1.wT .* (ρ_air * c_p), label="T1IRG")
ax1.plot(fx2.time, fx2.wT .* (ρ_air * c_p), label="T2IRG")
ax1.plot(fx5.time, fx5.wT .* (ρ_air * c_p), label="KAIJO", alpha=0.5)
ax1.tick_params(axis="x", labelbottom=false)
ax1.grid()
ax1.legend()
ax2 = fig.add_subplot(gs[2, 1], sharex=ax1)
#ax2.plot(fx1.time, fx1.wq .* (L_v * 1e-3), label="T1IRG")
ax2.plot(fx2.time, fx2.wq .* (L_v * 1e-3), label="T2IRG")
ax2.set_ylabel(L"\overline{w'q'}~[\mathrm{W~m^{-2}}]")
date_format = pydates.DateFormatter("%H:%M")
ax2.xaxis.set_major_formatter(date_format)
ax2.set_xlabel("31.05.2021")
ax2.grid()

######################################################
##
#scatterplot winddir Tdash
ecdf = evaldf3
ecname = "T2LCSAT"

scatterx = mod.(rad2deg.(atan.(-ecdf.v, ecdf.u)), 360)
scattery = ecdf.T .- gen.movingaverage(ecdf.T, 1200)

sortpermu = sortperm(scatterx)
scatterxsort = scatterx[sortpermu]
scatterysort = scattery[sortpermu]

medianbinwidth = 2
steps = collect(0:medianbinwidth:360)
medianofbins = zeros(Float64, length(steps) - 1)

for i in 2:length(steps)
    to_median = scatterysort[steps[i-1].<=scatterxsort.<steps[i]]
    if length(to_median) != 0
        medianofbins[i-1] = median(to_median)
    else
        medianofbins[i-1] = NaN
    end
end

cmap = PyPlot.get_cmap("tab10")
fig = PyPlot.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(ecname)
ax1.set_xlabel("winddir (360=up-valley, 270=west towards screen)")
ax1.set_ylabel(L"T-\overline{T}_{\mathrm{1min}}~\mathrm{[K]}")
ax1.scatter(scatterx, scattery, s=1)
ax1.plot(steps[1:end-1] .+ medianbinwidth / 2, medianofbins, label=L"2^\circ\mathrm{-bin~median}", color=cmap(1))
ax1.grid()
ax1.legend()
##
######################################################
##
#scatterplot winddir fluxraw
ecdf = evaldf6
fxrawdf = fx6_raw
ecname = "TJK"

scatterx = mod.(rad2deg.(atan.(-ecdf.v, ecdf.u)), 360)
scattery = fxrawdf.wT

sortpermu = sortperm(scatterx)
scatterxsort = scatterx[sortpermu]
scatterysort = scattery[sortpermu]

medianbinwidth = 2
steps = collect(0:medianbinwidth:360)
medianofbins = zeros(Float64, length(steps) - 1)

for i in 2:length(steps)
    to_median = scatterysort[steps[i-1].<=scatterxsort.<steps[i]]
    if length(to_median) != 0
        medianofbins[i-1] = median(to_median)
    else
        medianofbins[i-1] = NaN
    end
end

cmap = PyPlot.get_cmap("tab10")
fig = PyPlot.figure()
ax1 = fig.add_subplot(111)
ax1.set_title(ecname)
ax1.set_xlabel("winddir (360=up-valley, 270=west towards screen)")
ax1.set_ylabel(L"\overline{w'T'}_{\mathrm{raw}}~\mathrm{[K~m~s^{-1}]}")
ax1.scatter(scatterx, scattery, s=1)
ax1.plot(steps[1:end-1] .+ medianbinwidth / 2, medianofbins, label=L"2^\circ\mathrm{-bin~median}", color=cmap(1))
ax1.grid()
ax1.legend()
##
######################################################
#scatterplot with mean u and mean T
meanu2 = gen.movingaverage(evaldf2.u, round(Int, Millisecond(Second(100)) / timestep))
meanT = gen.movingaverage(evaldf2.T, round(Int, Millisecond(Second(100)) / timestep))
dirafterdr = 180 .- mod.(atan.(-evaldf2.v, evaldf2.u), 180)

PyPlot.pygui(true)
fig = PyPlot.figure(figsize=(13, 6))
gs = gridspec.GridSpec(1, 3, width_ratios=[15, 15, 1])
fig.suptitle(string(evalstart, " - ", evalend, " T2IRG"))
ax1 = fig.add_subplot(gs[1, 1])
ax1.grid()
ax1.set_ylabel(L"\overline{w'\theta'}~ [\mathrm{K~m~s^{-1}}]")
ax1.set_xlabel(L"\overline{u}~[\mathrm{m~s^{-1}}]")
#ax1.set_xlim([-4,6])
#ax1.set_ylim([-0.1,0.15])
#ax1.set_ylim([-0.125,0.075])
sp = ax1.scatter(meanu2[1:500:end], fx2.wT[1:500:end], s=0.5, c=pydates.date2num.(fx4.time[1:500:end]), cmap="turbo")
ax2 = fig.add_subplot(gs[1, 2], sharey=ax1)
#ax2.set_xlim([3,16])
#ax2.set_ylim([-0.1,0.15])
ax2.grid()
ax2.set_xlabel(L"\overline{T}~[\mathrm{K}]")
ax2.tick_params(axis="y", labelleft=false)
#ax2.set_ylabel(L"\overline{w'\theta'}~ [\mathrm{K~m~s^{-1}}]")
sp2 = ax2.scatter(meanT[1:500:end], fx2.wT[1:500:end], s=0.5, c=pydates.date2num.(fx4.time[1:500:end]), cmap="turbo")
ax3 = fig.add_subplot(gs[1, 3])
cb2 = PyPlot.colorbar(sp, cax=ax3)
cb2.set_label("date/time")
loc = pydates.AutoDateLocator()
cb2.ax.yaxis.set_major_locator(loc)
date_format = pydates.DateFormatter("%d.%m %H")
cb2.ax.yaxis.set_major_formatter(date_format)
#cb2.ax.yaxis.set_major_formatter(pydates.ConciseDateFormatter(loc))
##
######################################################
#scatterplot sensible and latent heat flux
PyPlot.pygui(true)
fig = PyPlot.figure()
gs = gridspec.GridSpec(1, 2, width_ratios=[30, 1])
fig.suptitle(string(evalstart, " - ", evalend, " T2IRG"))
ax1 = fig.add_subplot(gs[1, 1])
ax1.grid()
ax1.set_ylabel(L"\overline{w'\theta'}~[\mathrm{W~m^{-2}}]")
ax1.set_xlabel(L"\overline{w'q'}~[\mathrm{W~m^{-2}}]")
sp = ax1.scatter(fx2.wq[1:500:end] .* (L_v * 1e-3), fx2.wT[1:500:end] .* (ρ_air * c_p), s=0.5, c=pydates.date2num.(fx4.time[1:500:end]), cmap="turbo")
ax1.set_xlim([-50, 200])
ax1.set_ylim([-100, 200])
ax2 = fig.add_subplot(gs[1, 2])
cb2 = PyPlot.colorbar(sp, cax=ax2)
cb2.set_label("date/time")
loc = pydates.AutoDateLocator()
cb2.ax.yaxis.set_major_locator(loc)
date_format = pydates.DateFormatter("%d.%m %H")
cb2.ax.yaxis.set_major_formatter(date_format)
#cb2.ax.yaxis.set_major_formatter(pydates.ConciseDateFormatter(loc))
##
######################################################
#scatterplot latent heat fluxes IRG1 and IRG2
PyPlot.pygui(true)
fig = PyPlot.figure()
gs = gridspec.GridSpec(1, 2, width_ratios=[30, 1])
fig.suptitle(string(evalstart, " - ", evalend))
ax1 = fig.add_subplot(gs[1, 1])
ax1.grid()
ax1.set_ylabel(L"\overline{w'q'}_{T2IRG}~[\mathrm{W~m^{-2}}]")
ax1.set_xlabel(L"\overline{w'q'}_{T1IRG}~[\mathrm{W~m^{-2}}]")
sp = ax1.scatter(fx1.wq[1:500:end] .* (L_v * 1e-3), fx2.wq[1:500:end] .* (L_v * 1e-3), s=0.5, c=pydates.date2num.(fx4.time[1:500:end]), cmap="turbo") #meanu2[1:500:end]
#ax1.set_xlim([-50,200])
#ax1.set_ylim([-100,200])
ax2 = fig.add_subplot(gs[1, 2])
cb2 = PyPlot.colorbar(sp, cax=ax2)
cb2.set_label("date/time")#L"\overline{u}~\mathrm{[m~s^{-1}]}")
loc = pydates.AutoDateLocator()
cb2.ax.yaxis.set_major_locator(loc)
date_format = pydates.DateFormatter("%d.%m %H")
cb2.ax.yaxis.set_major_formatter(date_format)
#cb2.ax.yaxis.set_major_formatter(pydates.ConciseDateFormatter(loc))
##
######################################################
#scatterplot turbulent part between T1IRG and T2IRG
PyPlot.pygui(true)
fig = PyPlot.figure()
gs = gridspec.GridSpec(1, 2, width_ratios=[30, 1])
fig.suptitle(string(evalstart, " - ", evalend))
ax1 = fig.add_subplot(gs[1, 1])
ax1.grid()
ax1.set_ylabel(L"u'_{T2IRG}~[\mathrm{m~s^{-1}}]")
ax1.set_xlabel(L"u'_{T1IRG}~[\mathrm{m~s^{-1}}]")
avgidx = round(Int, Millisecond(ra1) / Millisecond(timestep))
sp = ax1.scatter(gen.movingaverage(evaldf1.u - gen.movingaverage(evaldf1.u, avgidx), avgidx)[1:500:end],
    gen.movingaverage(evaldf2.u - gen.movingaverage(evaldf2.u, avgidx), avgidx)[1:500:end], s=0.5, c=meanu2[1:500:end], cmap="turbo") #fx2.wT[1:500:end]
#ax1.set_xlim([-50,200])
#ax1.set_ylim([-100,200])
ax2 = fig.add_subplot(gs[1, 2])
cb2 = PyPlot.colorbar(sp, cax=ax2)
cb2.set_label(L"u~\mathrm{[m~s^{-1}]}")#\overline{w'T'}_{T2IRG}~\mathrm{[K~m~s^{-1}]}")
#cb2.ax.yaxis.set_major_formatter(pydates.ConciseDateFormatter(loc))
##
######################################################
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

######################################################
#check for different situations
##
#upvalley wind
windnorth = BitVector(undef, size(evaldf1, 1))

#parameters for statistics using the TJK station
upvalleywind_leftlim = 305 #degree
upvalleywind_rightlim = 10 #degree
downvalleywind_leftlim = 135#degree
downvalleywind_rightlim = 200#degree

firstentry = turb.findnearest(evaldf1.time[1], tjkmeteodata2.time)

starttime = evaldf1.time[firstentry[1]]
endtime = starttime + Minute(5)
@showprogress for idx in 1:size(tjkmeteodata2, 1)
    isnorth = upvalleywind_leftlim .<= tjkmeteodata2.wind_mean_vector_direction[idx] ||
              tjkmeteodata2.wind_mean_vector_direction[idx] .<= upvalleywind_rightlim
    windnorth[starttime.<=evaldf1.time.<endtime] .= isnorth
    starttime = endtime
    endtime = endtime + Minute(10)
end
##

T1IRGpos = fx1.wT .> 0   #positive flux at T1IRG
T2IRGneg = fx2.wT .< 0   #negative flux at T2IRG
T2LCSATpos = fx3.wT .> 0 #positive flux at T2LCSAT
T2UCSATpos = fx4.wT .> 0 #positive flux at T2UCSAT
T2CSATspos = T2UCSATpos .&& T2LCSATpos

#kaijodataexists
kaijodata = BitVector(undef, size(evaldf1, 1))
kaijodata .= false
#firstentry = turb.findnearest(evaldf1.time[1], evaldf5.time)
curridx = findall(x -> x == evaldf5.time[1], evaldf1.time)[1]
kaijodata[curridx] = true
@showprogress for idx in 2:size(evaldf5, 1)
    if evaldf5.time[idx] - evaldf5.time[idx-1] == Millisecond(50)
        curridx += 1
    else
        curridx = findfirst(x -> x == evaldf5.time[idx], evaldf1.time)[1]
    end
    kaijodata[curridx] = true
end

#plot
PyPlot.pygui(true)
fig = PyPlot.figure()
gs = gridspec.GridSpec(5, 1)
ax1 = fig.add_subplot(gs[1, 1])
ax1.plot(evaldf1.time[1:20:end], windnorth[1:20:end] .- 0.5, ".", markersize=0.5)
ax1.set_ylim([0.4, 0.6])
ax1.set_ylabel("upvalley wind")
ax1.tick_params(axis="x", labelbottom=false)
ax1.tick_params(axis="y", labelleft=false)
ax2 = fig.add_subplot(gs[2, 1], sharex=ax1)
ax2.plot(evaldf1.time[1:20:end], T1IRGpos[1:20:end] .- 0.5, ".", markersize=0.5)
ax2.set_ylim([0.4, 0.6])
ax2.set_ylabel("T1IRG +")
ax2.tick_params(axis="x", labelbottom=false)
ax2.tick_params(axis="y", labelleft=false)
ax3 = fig.add_subplot(gs[3, 1], sharex=ax1)
ax3.plot(evaldf2.time[1:20:end], T2IRGneg[1:20:end] .- 0.5, ".", markersize=0.5)
ax3.set_ylim([0.4, 0.6])
ax3.set_ylabel("T2IRG -")
ax3.tick_params(axis="x", labelbottom=false)
ax3.tick_params(axis="y", labelleft=false)
ax4 = fig.add_subplot(gs[4, 1], sharex=ax1)
ax4.plot(evaldf3.time[1:20:end], T2CSATspos[1:20:end] .- 0.5, ".", markersize=0.5)
ax4.set_ylim([0.4, 0.6])
ax4.set_ylabel("T2CSATs +")
ax4.tick_params(axis="x", labelbottom=false)
ax4.tick_params(axis="y", labelleft=false)
ax5 = fig.add_subplot(gs[5, 1], sharex=ax1)
ax5.plot(evaldf1.time[1:20:end], kaijodata[1:20:end] .- 0.5, ".", markersize=0.5)
ax5.set_ylim([0.4, 0.6])
ax5.set_ylabel("Kaijo")
ax5.tick_params(axis="y", labelleft=false)
ax5.set_xlabel("time")
fig.autofmt_xdate()

######################################################
#plot mean wind directions
avgl = 12000
m1u = gen.movingaverage(evaldf1.u, avgl)
m1v = gen.movingaverage(evaldf1.v, avgl)
m2u = gen.movingaverage(evaldf2.u, avgl)
m2v = gen.movingaverage(evaldf2.v, avgl)
m3u = gen.movingaverage(evaldf3.u, avgl)
m3v = gen.movingaverage(evaldf3.v, avgl)
m4u = gen.movingaverage(evaldf4.u, avgl)
m4v = gen.movingaverage(evaldf4.v, avgl)
m5u = gen.movingaverage(evaldf5.u, avgl)
m5v = gen.movingaverage(evaldf5.v, avgl)

dir1 = turb.winddir.(m1u, m1v)
dir2 = turb.winddir.(m2u, m2v)
dir3 = turb.winddir.(m3u, m3v)
dir4 = turb.winddir.(m4u, m4v)
dir5 = turb.winddir.(m5u, m5v)

tjkmean = mean(tjkmeteodata2.wind_mean_vector_direction)
sh1 = 0 #tjkmean - mean(dir1)
sh2 = 0 #tjkmean - mean(dir2)
sh3 = 0 #tjkmean - mean(dir3)
sh4 = 0 #tjkmean - mean(dir4)
sh5 = 0 #tjkmean - mean(dir5)

PyPlot.pygui(true)
wfig = PyPlot.figure()
wax = wfig.add_subplot(111)
wax.set_title("Sonics wind direction")
wax.set_xlabel("time")
wax.set_ylabel(L"dir~\mathrm{[^\circ C]}")
wax.plot(tjkmeteodata2.time, tjkmeteodata2.wind_mean_vector_direction, label="TJK")
wax.plot(evaldf1.time, mod.(dir1 .+ sh1, 360), label="T1IRG", alpha=0.7)
wax.plot(evaldf2.time, mod.(dir2 .+ sh2, 360), label="T2IRG", alpha=0.7)
wax.plot(evaldf3.time, mod.(dir3 .+ sh3, 360), label="T2LCSAT", alpha=0.7)
wax.plot(evaldf4.time, mod.(dir4 .+ sh4, 360), label="T2UCSAT", alpha=0.7)
wax.plot(evaldf5.time, mod.(dir5 .+ sh5, 360), label="KAIJO", alpha=0.7)
wax.legend()
wax.grid()
wfig.autofmt_xdate()
##

PyPlot.pygui(true)
wsfig = PyPlot.figure()
wsax = wsfig.add_subplot(111)
wsax.set_title("Windspeed - 10min avg")
wsax.set_ylabel(L"\sqrt{u^2+v^2}~\mathrm{[m~s^{-1}]}")
wsax.set_xlabel("10.06.2021")
#wsax.plot(tjkmeteodata2.time, tjkmeteodata2.wind_mean_scalar, alpha=0.5, label="TJK")
wsax.plot(evaldf1.time, gen.movingaverage(sqrt.(evaldf1.u .^ 2 + evaldf1.v .^ 2), 12000), label="T1IRG")
wsax.plot(evaldf2.time, gen.movingaverage(sqrt.(evaldf2.u .^ 2 + evaldf2.v .^ 2), 12000), label="T2IRG")
wsax.plot(evaldf3.time, gen.movingaverage(sqrt.(evaldf3.u .^ 2 + evaldf3.v .^ 2), 12000), label="T2LCSAT")
wsax.plot(evaldf4.time, gen.movingaverage(sqrt.(evaldf4.u .^ 2 + evaldf4.v .^ 2), 12000), label="T2UCSAT")
wsax.plot(evaldf5.time, gen.movingaverage(sqrt.(evaldf5.u .^ 2 + evaldf5.v .^ 2), 12000), label="KAIJO")
wsax.legend()
wsax.grid()
date_format = pydates.DateFormatter("%H:%M")
wsax.xaxis.set_major_formatter(date_format)
fig.autofmt_xdate()
##

######################################################
###               OBUKHOV-LENGTH                   ###
######################################################
#multi-panel plot to investigate Obhukov-length
turb.missing2nan!(evaldf4)
PyPlot.pygui(true)
fig = PyPlot.figure()
gs = gridspec.GridSpec(4, 1)
axl = fig.add_subplot(gs[4, 1])
axl.plot(evaldf4.time, gen.movingaverage(evaldf4.T, 600 * 20))
axl.set_ylabel(L"T_v~\mathrm{[^\circ C]}")
axl.set_xlabel(L"31.05.2021")
axl.grid()
ax1 = fig.add_subplot(gs[1, 1], sharex=axl)
ax1.plot(L4.time, L4.L)
ax1.set_ylabel(L"L~\mathrm{[m]}")
ax1.grid()
ax1.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(gs[2, 1], sharex=axl)
ax2.plot(fx4.time, fx4.u_star .^ 3)
ax2.set_ylabel(L"u_\ast^3~\mathrm{[m^3~s^{-3}]}")
ax2.grid()
ax2.tick_params(axis="x", labelbottom=false)
ax3 = fig.add_subplot(gs[3, 1], sharex=axl)
ax3.plot(fx4.time, fx4.wT)
ax3.set_ylabel(L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}")
ax3.grid()
ax3.tick_params(axis="x", labelbottom=false)
date_format = pydates.DateFormatter("%H:%M")
axl.xaxis.set_major_formatter(date_format)

##
######################################################
L1[Time(08, 00, 00).>Time.(L1.time).||Time(18, 00, 00).<Time.(L1.time), :L] .= NaN
L2[Time(08, 00, 00).>Time.(L2.time).||Time(18, 00, 00).<Time.(L2.time), :L] .= NaN
L3[Time(08, 00, 00).>Time.(L3.time).||Time(18, 00, 00).<Time.(L3.time), :L] .= NaN
L4[Time(08, 00, 00).>Time.(L4.time).||Time(18, 00, 00).<Time.(L4.time), :L] .= NaN
L5[Time(08, 00, 00).>Time.(L5.time).||Time(18, 00, 00).<Time.(L5.time), :L] .= NaN
L6[Time(08, 00, 00).>Time.(L6.time).||Time(18, 00, 00).<Time.(L6.time), :L] .= NaN

L1[Date.(L1.time) .== Date(2021,05,23), :L] .= NaN
L2[Date.(L2.time) .== Date(2021,05,23), :L] .= NaN
L3[Date.(L3.time) .== Date(2021,05,23), :L] .= NaN
L4[Date.(L4.time) .== Date(2021,05,23), :L] .= NaN
L5[Date.(L5.time) .== Date(2021,05,23), :L] .= NaN
L6[Date.(L6.time) .== Date(2021,05,23), :L] .= NaN

L1[DateTime(2021,05,22,08,00,00) .<= L1.time .<= DateTime(2021,05,22,12,15,00), :L] .= NaN
L2[DateTime(2021,05,22,08,00,00) .<= L2.time .<= DateTime(2021,05,22,12,15,00), :L] .= NaN
L3[DateTime(2021,05,22,08,00,00) .<= L3.time .<= DateTime(2021,05,22,12,15,00), :L] .= NaN
L4[DateTime(2021,05,22,08,00,00) .<= L4.time .<= DateTime(2021,05,22,12,15,00), :L] .= NaN
L5[DateTime(2021,05,22,08,00,00) .<= L5.time .<= DateTime(2021,05,22,12,15,00), :L] .= NaN
L6[DateTime(2021,05,22,08,00,00) .<= L6.time .<= DateTime(2021,05,22,12,15,00), :L] .= NaN

L1[DateTime(2021,05,24,14,00,00) .<= L1.time .<= DateTime(2021,05,24,18,15,00), :L] .= NaN
L2[DateTime(2021,05,24,14,00,00) .<= L2.time .<= DateTime(2021,05,24,18,15,00), :L] .= NaN
L3[DateTime(2021,05,24,14,00,00) .<= L3.time .<= DateTime(2021,05,24,18,15,00), :L] .= NaN
L4[DateTime(2021,05,24,14,00,00) .<= L4.time .<= DateTime(2021,05,24,18,15,00), :L] .= NaN
L5[DateTime(2021,05,24,14,00,00) .<= L5.time .<= DateTime(2021,05,24,18,15,00), :L] .= NaN
L6[DateTime(2021,05,24,14,00,00) .<= L6.time .<= DateTime(2021,05,24,18,15,00), :L] .= NaN

L6[DateTime(2021,05,25,08,00,00) .<= L6.time .<= DateTime(2021,05,25,08,30,00), :L] .= NaN


nrdays = Day(L2.time[end] - L2.time[1])
zLq1 = fill(NaN, Dates.value(nrdays), 6)
zLmed = fill(NaN, Dates.value(nrdays), 6)
zLq3 = fill(NaN, Dates.value(nrdays), 6)
startdate = Date(L2.time[1])
zLdates = collect(DateTime(startdate)+Hour(12):Day(1):DateTime(startdate + nrdays - Day(1))+Hour(12))

for i in 1:size(zLq1, 1)
    date = startdate + Day(i - 1)
    #(zLq1[i, 1], zLmed[i, 1], zLq3[i, 1]) = quantile(filter(!isnan, L1[Date.(L1.time).==date,:L]), [0.25, 0.5, 0.75])
    try
        (zLq1[i, 2], zLmed[i, 2], zLq3[i, 2]) = quantile(filter(!isnan, 1 ./ L2[Date.(L2.time).==date, :L]), [0.25, 0.5, 0.75])
    catch e
    end
    try
        (zLq1[i, 3], zLmed[i, 3], zLq3[i, 3]) = quantile(filter(!isnan, 2.2 ./ L3[Date.(L3.time).==date,:L]), [0.25, 0.5, 0.75])
    catch e
    end
    try
        (zLq1[i, 4], zLmed[i, 4], zLq3[i, 4]) = quantile(filter(!isnan, 2.9 ./ L4[Date.(L4.time).==date,:L]), [0.25, 0.5, 0.75])       
    catch e
    end
    #(zLq1[i, 5], zLmed[i, 5], zLq3[i, 5]) = quantile(filter(!isnan, L5[Date.(L5.time).==date,:L]), [0.25, 0.5, 0.75])
    try
        (zLq1[i, 6], zLmed[i, 6], zLq3[i, 6]) = quantile(filter(!isnan, 5.0 ./ L6[Date.(L6.time).==date, :L]), [0.25, 0.5, 0.75])
    catch e
    end
end

##
#stability plot (Obukhov-length)
limmin = -5
limmax = 5
PyPlot.pygui(true)
fig = PyPlot.figure(figsize=(10,6))
gs = gridspec.GridSpec(1, 1)
axl = fig.add_subplot(gs[1, 1])
#axl.set_title("T2IRG")
#axl.plot(L2.time[1:100:end], 1.2 ./ L2.L[1:100:end], alpha=0.3)
#axl.plot(L6.time[1:1200:end], 5 ./ L6.L[1:1200:end], alpha=0.3)
axl.set_ylabel(L"\frac{z}{L}}", fontsize=16)
axl.set_ylim(limmin, limmax)
axl.grid()
axl.errorbar(zLdates, zLmed[:, 2], yerr=transpose(abs.(zLmed[:, 2] .- [zLq1[:, 2] zLq3[:, 2]])), marker="s", ms=4, capsize=2, label="T2IRG")
#axl.errorbar(zLdates, zLmed[:, 3], yerr=transpose(abs.(zLmed[:, 3] .- [zLq1[:, 3] zLq3[:, 3]])), marker="s", ms=4, capsize=2, label="T2LCSAT")
axl.errorbar(zLdates, zLmed[:, 4], yerr=transpose(abs.(zLmed[:, 4] .- [zLq1[:, 4] zLq3[:, 4]])), marker="s", ms=4, capsize=2, label="T2UCSAT")
axl.errorbar(zLdates, zLmed[:, 6], yerr=transpose(abs.(zLmed[:, 6] .- [zLq1[:, 6] zLq3[:, 6]])), marker="s", ms=4, capsize=2, label="TJK")
axl.legend()
majorloc = pydates.DayLocator()
minorloc = pydates.HourLocator(interval=12)
majorformatter = pydates.DateFormatter("%d.%m")
minorformatter = pydates.DateFormatter("%H")
axl.xaxis.set_major_locator(majorloc)
axl.xaxis.set_minor_locator(minorloc)
axl.xaxis.set_major_formatter(majorformatter)
for label in axl.xaxis.get_ticklabels()[1:2:end]
    label.set_visible(false)
end

#=
ax1 = fig.add_subplot(gs[1, 1], sharex=axl)
ax1.set_title("TJK")
ax1.set_ylabel(L"\frac{z}{L}}")
ax1.grid()
ax1.set_ylim(limmin, limmax)
ax1.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(gs[2, 1], sharex=axl)
ax2.set_title("T2UCSAT")
ax2.plot(L4.time[1:1200:end], 5 ./ L4.L[1:1200:end])
ax2.set_ylabel(L"\frac{z}{L}}")
ax2.grid()
ax2.set_ylim(limmin, limmax)
ax2.tick_params(axis="x", labelbottom=false)
ax3 = fig.add_subplot(gs[3, 1], sharex=axl)
ax3.set_title("T2LCSAT")
ax3.plot(L3.time[1:1200:end], 5 ./ L3.L[1:1200:end])
ax3.set_ylabel(L"\frac{z}{L}}")
ax3.grid()
ax3.set_ylim(limmin, limmax)
ax3.tick_params(axis="x", labelbottom=false)
=#
##
h2 = 1.1
h3 = 2.0
h4 = 3.0
h6 = 5.0

medtimestart = DateTime(2021, 05, 31, 12, 00, 00)
medtimeend = DateTime(2021, 05, 31, 16, 00, 00)
L2totake = L2[medtimestart.<=L2.time.<=medtimeend, :]
L3totake = L3[medtimestart.<=L3.time.<=medtimeend, :]
L4totake = L4[medtimestart.<=L4.time.<=medtimeend, :]
L6totake = L6[medtimestart.<=L6.time.<=medtimeend, :]

(zL2q1, zL2med, zL2q3) = quantile(h2 ./ L2totake.L, [0.25, 0.5, 0.75])
(zL3q1, zL3med, zL3q3) = quantile(h3 ./ L3totake.L, [0.25, 0.5, 0.75])
(zL4q1, zL4med, zL4q3) = quantile(h4 ./ L4totake.L, [0.25, 0.5, 0.75])
(zL6q1, zL6med, zL6q3) = quantile(h6 ./ L6totake.L, [0.25, 0.5, 0.75])

#stability-profile plot
fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.errorbar(zL2med, 1.1, xerr=transpose(abs.(zL2med .- [zL2q1 zL2q3])), fmt="o")
ax.errorbar(zL3med, 2.0, xerr=transpose(abs.(zL3med .- [zL3q1 zL3q3])), fmt="o")
ax.errorbar(zL4med, 3.0, xerr=transpose(abs.(zL4med .- [zL4q1 zL4q3])), fmt="o")
#ax.errorbar(zL6med, 5.0, xerr=transpose(abs.(zL6med .- [zL6q1 zL6q3])), fmt="o")

##

######################################################
###            SONIC TEMPERATURE OFFSET            ###
######################################################
#sonic temperature offset
#avgT1 = mean(evaldf1.T)
#avgT2 = mean(evaldf2.T)

#T2minusT1 = avgT2 - avgT1

#plot offsets
PyPlot.pygui(true)
fig2 = PyPlot.figure()
ax2 = fig2.add_subplot(111)
ax2.set_title(L"\mathrm{T_{sonic}~offset}")
ax2.set_xlabel("10.06.2021")
ax2.set_ylabel(L"T_{sonic}~\mathrm{[^\circ C]}")
ax2.grid()
ax2.plot(evaldf1.time, gen.movingaverage(evaldf1.T, 12000), label="T1IRG orig")
ax2.plot(evaldf2.time, gen.movingaverage(evaldf2.T, 12000), label="T2IRG")
ax2.plot(evaldf3.time, gen.movingaverage(evaldf3.T, 12000), label="T2LCSAT")
ax2.plot(evaldf4.time, gen.movingaverage(evaldf4.T, 12000), label="T2UCSAT")
ax2.plot(evaldf5.time, gen.movingaverage(evaldf5.T, 12000), label="KAIJO")
#ax2.plot(evaldf1.time, gen.movingaverage(evaldf1.T, 12000).+T2minusT1, label="T1IRG corr")
ax2.legend()
date_format = pydates.DateFormatter("%H:%M")
ax2.xaxis.set_major_formatter(date_format)

#evaldf1mod = evaldf1
#evaldf1mod.T = evaldf1.T .+ T2minusT1
##
#is it moisture related?
#calculate q
q1 = evaldf1.h2o .* (1e-3 / ρ_air)
q2 = evaldf2.h2o .* (1e-3 / ρ_air)

#plot q
qfig = PyPlot.figure()
qax = qfig.add_subplot(111)
qax.set_xlabel("10.06.2021")
qax.set_ylabel(L"q~\mathrm{\left[\cdot 10^{-3}\frac{(g~m^{-3})_{H_2O}}{(g~m^{-3})_{air}}\right]}")
qax.plot(evaldf1.time, gen.movingaverage(q1, 12000) .* 1e3, label="T1IRG")
qax.plot(evaldf2.time, gen.movingaverage(q2, 12000) .* 1e3, label="T2IRG")
qax.grid()
qax.legend()
date_format = pydates.DateFormatter("%H:%M")
qax.xaxis.set_major_formatter(date_format)

#correct sonic temperature for moisture effects
T1_moistcorr = evaldf1.T .* (0.51 * q1 .+ 1) .^ (-1)
T2_moistcorr = evaldf2.T .* (0.51 * q2 .+ 1) .^ (-1)

#plot corrected/not corrected
fig2 = PyPlot.figure()
ax2 = fig2.add_subplot(111)
ax2.set_title(L"\mathrm{T_{sonic}~moisture~correction}")
ax2.set_xlabel("10.06.2021")
ax2.set_ylabel(L"T_{sonic}~\mathrm{[^\circ C]}")
ax2.grid()
ax2.plot(evaldf1.time, gen.movingaverage(evaldf1.T, 12000), label="T1IRG orig")
ax2.plot(evaldf2.time, gen.movingaverage(evaldf2.T, 12000), label="T2IRG orig")
ax2.plot(evaldf1.time, gen.movingaverage(T1_moistcorr, 12000), label="T1IRG corr")
ax2.plot(evaldf2.time, gen.movingaverage(T2_moistcorr, 12000), label="T2IRG corr")
ax2.legend()
date_format = pydates.DateFormatter("%H:%M")
ax2.xaxis.set_major_formatter(date_format)

#-> too small

avgT1 = gen.movingaverage(evaldf1night.T, 12000)
avgT2 = gen.movingaverage(evaldf2night.T, 12000)

#check if data over night gives hints
PyPlot.pygui(true)
fig3 = PyPlot.figure()
ax3 = fig3.add_subplot(111)
ax3.set_title(L"\mathrm{T_{sonic}~nights}")
ax3.set_ylabel(L"T_{sonic}~\mathrm{[^\circ C]}")
ax3.grid()
ax3.plot(evaldf1night.time[1:500:end], avgT1[1:500:end], label="T1IRG")
ax3.plot(evaldf2night.time[1:500:end], avgT2[1:500:end], label="T2IRG")
ax3.legend()
date_format = pydates.DateFormatter("%d.%m %H:%M")
ax3.xaxis.set_major_formatter(date_format)
fig3.autofmt_xdate()
@show(mean(evaldf1night.T) - mean(evaldf2night.T))

#check for wind
spd1 = sqrt.(evaldf1night.u .^ 2 + evaldf1night.v .^ 2)
spd2 = sqrt.(evaldf2night.u .^ 2 + evaldf2night.v .^ 2)

hist1fig = PyPlot.figure()
hax = hist1fig.add_subplot(111)
hax.set_title("night wind speeds - 8.6.-11.6.")
hax.set_xlabel(L"\sqrt{u^2+v^2}~\mathrm{[m~s^{-1}]}")
hax.grid()
hax.hist(spd1, bins=collect(0:0.05:6), density=true, alpha=0.5, label="T1IRG")
hax.hist(spd2, bins=collect(0:0.05:6), density=true, alpha=0.5, label="T2IRG")
hax.legend()

spd_thresh = 3 #m/s threshold to be considered as "wind mixed night"

#average
spd1_avg = gen.movingaverage(spd1, 12000)
spd2_avg = gen.movingaverage(spd2, 12000)

#select temperature values for wind above threshold
cond = spd1_avg .>= spd_thresh .&& spd2_avg .>= spd_thresh

#histogram for windy nights
windyhistfig = PyPlot.figure()
whax = windyhistfig.add_subplot(111)
whax.set_title("8.6.-11.6. - nights with spd>3m/s")
#whax.set_xlabel(L"\left(\overline{T_{T1IRG}}\right)_{10min}-\left(\overline{T_{T2IRG}}\right)_{10min}")
whax.set_xlabel(L"T_{T1IRG}-T_{T2IRG}")
whax.grid()
#whax.hist(evaldf1night.T - evaldf2night.T, bins=collect(-2:0.01:2), density=true, alpha=0.5, label=L"20~\mathrm{Hz}")
#whax.hist(avgT1 - avgT2, bins=collect(-2:0.01:2), density=true, alpha=0.5, label=L"10~\mathrm{min~avg.}")
whax.hist((avgT1-avgT2)[cond], bins=collect(-0.5:0.001:-0.3), density=true, label=L"10~\mathrm{min~avg.}")
whax.legend()
##
######################################################
###                  ADVECTION                     ###
######################################################
@warn("Rotating of the wind direction (maybe) needs to be done!")

#from Mott20(HEFEX): HA = -\frac{ΔT}{Δy}̅v with ̅v=(v₁+v₂)/2; ΔT = T₁-T₂

evaldfadv1 = evaldf1
evaldfadv2 = evaldf2
#correction due to moisture
#calculate q
q1 = evaldfadv1.h2o .* (1e-3 / ρ_air)
q2 = evaldfadv2.h2o .* (1e-3 / ρ_air)
#correct sonic temperature for moisture effects
evaldfadv1.T = evaldf1.T .* (0.51 * q1 .+ 1) .^ (-1)
evaldfadv2.T = evaldf2.T .* (0.51 * q2 .+ 1) .^ (-1)
#=
#apply modification to T1IRG
evaldf1mod = evaldf1
evaldf1mod.T = evaldf1mod.T .+ 0.58053369496
=#
irgadv = turb.advect(evaldfadv1, evaldfadv2, "T", "u", ra1, 30, 6.15)

irgadv.adv = gen.movingaverage(irgadv.adv, round(Int, Millisecond(ra1) / timestep))

#set advection to NaN if wind from south
#irgadv.adv[irgadv.meanwind.>0].=NaN
#irgadv.Δq[irgadv.meanwind.>0].=NaN

#Plot advection time series
##
PyPlot.pygui(true)
fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.set_title("Advection of sensible heat between T1IRG and T2IRG")
ptirg = ax.plot(irgadv.time[1:end], irgadv.adv[1:end], label="moisture corr")
ax.set_ylabel(L"HA~\mathrm{[K~m^{-1}~s^{-1}]}")
#ax.set_ylim([-0.25,0.25])
ax.xaxis_date()
#majorlocator = pydates.HourLocator([0,12])#interval=6)
#minorlocator = pydates.HourLocator(interval=1)
#ax.xaxis.set_major_locator(majorlocator)
#ax.xaxis.set_minor_locator(minorlocator)
date_format = pydates.DateFormatter("%H:%M")
ax.xaxis.set_major_formatter(date_format)
ax.set_xlabel("10.06.2021")
#fig.autofmt_xdate()
#ax.set_xlim(DateTime(2021,05,21,07,00,00), DateTime(2021,05,24,22,00,00))
ax.grid()
ax.legend()

#Plot to understand advection
##
PyPlot.pygui(true)
fig = PyPlot.figure()
ax = fig.add_subplot(311)
ax.set_title("Advection of sensible heat between T1IRG and T2IRG")
ptirg = ax.plot(irgadv.time[1:500:end], irgadv.adv[1:500:end], label="T1->T2 (IRG)")
ax.set_ylabel(L"HA~\mathrm{[K~m^{-1}~s^{-1}]}")
ax.legend()
ax.grid()
ax.xaxis_date()
ax.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(312, sharex=ax)
ptmw = ax2.plot(irgadv.time[1:500:end], irgadv.meanwind[1:500:end], label="mean wind")
ax2.set_ylabel(L"̅\overline{u}~\mathrm{[m~s^{-1}]}")
ax2.legend()
ax2.grid()
ax2.tick_params(axis="x", labelbottom=false)
ax3 = fig.add_subplot(313, sharex=ax)
ptmw = ax3.plot(irgadv.time[1:500:end], irgadv.Δq[1:500:end], label=L"\Delta q")
ax3.set_ylabel(L"\Delta T = T_1-T_2~\mathrm{[K]}")
ax3.legend()
ax3.grid()
#majorlocator = pydates.HourLocator([0,12])#interval=6)
#minorlocator = pydates.HourLocator(interval=1)
#ax3.xaxis.set_major_locator(majorlocator)
#ax3.xaxis.set_minor_locator(minorlocator)
date_format = pydates.DateFormatter("%H:%M")
ax3.xaxis.set_major_formatter(date_format)
#fig.autofmt_xdate()
ax3.set_xlabel("10.06.2021")
##
##
PyPlot.pygui(true)
fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.set_title("Advection of sensible heat between T1IRG and T2IRG")
ptirg = ax.scatter(irgadv.time[1:500:end], irgadv.adv[1:500:end], s=0.5, c=irgadv.Δq[1:500:end], cmap="turbo", label="T1->T2 (IRG)")
cbar = PyPlot.colorbar(ptirg, ax=ax)
cbar.set_label(L"\Delta T = T_1-T_2~ \mathrm{[K]}")
ax.set_ylabel(L"HA~\mathrm{[K~m^{-1}~s^{-1}]}")
ax.grid()
ax.xaxis_date()
majorlocator = pydates.HourLocator([0, 12])#interval=6)
minorlocator = pydates.HourLocator(interval=1)
ax.xaxis.set_major_locator(majorlocator)
ax.xaxis.set_minor_locator(minorlocator)
date_format = pydates.DateFormatter("%d.%m %H:%M")
ax.xaxis.set_major_formatter(date_format)
fig.autofmt_xdate()
##
######################################################
###                  PLOTTING                      ###
######################################################
#Plot advection and turb flux
##
PyPlot.pygui(true)
fig = PyPlot.figure()
ax = fig.add_subplot(211)
#ax.set_title("Advection of sensible heat between T1IRG and T2IRG")
ptirg = ax.plot(irgadv.time[1:500:end], irgadv.adv[1:500:end], label="T1->T2 (IRG)")
ax.set_ylabel(L"HA~\mathrm{[K~m^{-1}~s^{-1}]}")
#ax.set_ylim([-0.25,0.25])
ax.xaxis_date()
#ax.set_xlim(DateTime(2021,05,21,07,00,00), DateTime(2021,05,24,22,00,00))
ax.grid()
ax.legend()
ax2 = fig.add_subplot(212, sharex=ax)
ax2.set_title("")
ax2.plot(fx2.time[1:500:end], fx2.wT[1:500:end], label="T2IRG")
ax2.set_ylabel(L"\overline{w'\theta'}~\mathrm{[K~m~s^{-1}]}")
#ax2.set_ylim([-0.25,0.25])
majorlocator = pydates.HourLocator([0, 12])#interval=6)
minorlocator = pydates.HourLocator(interval=1)
ax2.xaxis.set_major_locator(majorlocator)
ax2.xaxis.set_minor_locator(minorlocator)
date_format = pydates.DateFormatter("%d.%m %H:%M")
ax2.xaxis.set_major_formatter(date_format)
fig.autofmt_xdate()
ax2.grid()
ax2.legend()
##
#scatter plot of HA and w'T'
PyPlot.pygui(true)
scatfig = PyPlot.figure()
ax = scatfig.add_subplot(111)
ax.set_title(string(evaldf1.time[1], " - ", evaldf1.time[end], " only up-valley"))
sc = ax.scatter(irgadv.adv, fx2.wT, c=irgadv.meanwind, cmap="turbo", s=1)
cb = PyPlot.colorbar(sc, ax=ax)
cb.set_label(L"\overline{u}~\mathrm{[m~s^{-1}]}")
ax.set_xlabel(L"HA_{T_1 \rightarrow T_2}~\mathrm{[K~m^{-1}~s^{-1}]}")
ax.set_ylabel(L"\overline{w'\theta'}_{T2IRG}~\mathrm{[K~m~s^{-1}]}")
ax.grid()

######################################################