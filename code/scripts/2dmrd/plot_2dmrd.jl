######################################################
###                   PLOT 2D-MRD                  ###
###            author: Michi Haugeneder            ###
######################################################
#=
Load data calculated with '2dmrd_distributed.jl' and
plot.
=#
using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings, DataFrames
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
cramericm = pyimport("cmcrameri.cm")

importdir = joinpath(@__DIR__, "..", "..")
datapath = "/home/haugened/Documents/data/"
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

######################################################
#variables

readfilename = joinpath(datapath, "2dmrd", "db_tot", "raw", "t2lcsat_uu_norm.nc")
#Reynolds averaging time for flux
#ra = Second(100)
norm = false #decide if data shall be normalized (true) or not (false)
exportnormed = false #write normed mrd_data to .netcdf-file

#cmpfile1 = joinpath(datapath, "2dmrd", "31051100_to_08060900_T2UCSAT_block1h.nc")
#cmpfile2 = joinpath(datapath, "2dmrd", "31051100_to_08060900_T2IRG_block1h.nc")

######################################################

println()
println("-----------S-T-A-R-T-------------")

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Millisecond(2^10*timestep),
    "KAIJO" => Millisecond(2^9*timestep))

######################################################
###                LOAD AND PLOT                   ###
######################################################
(time_middle, x, mrd_data, meanwind, evalstart, evalend,
    cols_for_mrd, nrpoints, blocklen, shift, sourcefile) =
    MRD.read2dmrdfromnetcdf(readfilename)

figtitle = MRD.create2dmrdtitle(sourcefile, evalstart, evalend)
#=
evalstart = DateTime(2021,05,29,00,00,00)
evalend = DateTime(2021,06,06,00,00,00)
mrd_data = mrd_data[1:end-49,12167:23686]
x = x[1:end-49,12167:23686]
time_middle = time_middle[12167:23686]
=#
fxavgperiod = blocklen * (evalend - evalstart)

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[evalstart.<=tjkmeteo[:, 1].<=evalend, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 1
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 2

#load data and calc fluxes
inst = MRD.instrumentnamefromsourcefile(sourcefile)
fluxfile = inst2file[inst]
ra = inst2ra[inst]

fluxdf = turb.readturbasnetcdf(fluxfile, evalstart, evalend)
#show statistics about missing data
turb.printmissstats(fluxdf)
#interpolate missing
fluxdf = turb.interpolatemissing(fluxdf)
#double rotation
if inst == "TJK"
    turb.drdf!(fluxdf, periodwise=false)
else
    turb.drdf!(fluxdf)
end
turb.missing2nan!(fluxdf)
wndspd = gen.movingaverage(sqrt.(fluxdf.u .^ 2 .+ fluxdf.v .^ 2 .+ fluxdf.w .^ 2), 20 * 3600)
fx1_raw = turb.turbflux(fluxdf, ra)
fx1 = turb.avgflux(fx1_raw, fxavgperiod)
fx1_raw = nothing
#Obukhov length
#L1 = turb.obukhov(fx1, blocklen * (evalend - evalstart))

#normalization
normfct = zeros(Float64, size(mrd_data, 2))
if norm
    #calculate sum of dyadic scales up to 2^11 for 2DMRD
    if inst == "Kaijo"
        scalestofind = Millisecond(50) .* 2 .^ collect(1:9)
    elseif inst == "TJK"
        scalestofind = Millisecond(50) .* 2 .^ collect(1:10)
    else
        scalestofind = Millisecond(50) .* 2 .^ collect(1:11)
    end
    idxscales = zeros(Int64, length(scalestofind))
    for ix in 1:length(scalestofind)
        idxscales[ix] = findall(x->x==scalestofind[ix], x[:,7])[1]
    end
    sumuptoscale = zeros(Float64, size(x, 2))
    for ix in 1:size(x, 2)
        sumuptoscale[ix] = sum(mrd_data[idxscales, ix])
    end

    #find correct column in flux file
    fx_str = string(cols_for_mrd[1], cols_for_mrd[2])
    fx_str_ret = string(cols_for_mrd[2], cols_for_mrd[1])
    if in(fx_str, names(fx1))
    elseif in(fx_str_ret, names(fx1))
        fx_str = fx_str_ret
    else
        @warn("Norming: Fluxes not found in fluxfile! Applying no norm!")
    end
 
    #take reference flux values for the times needed
    fxref = zeros(Float64, length(time_middle))
    idxfirst = findfirst(x->x==time_middle[1], fx1.time)
    fxref[1] = fx1[idxfirst, fx_str]
    shiftidcs = shift * (blocklen*(evalend-evalstart))/timestep
    for ix in 2:length(fxref)
        fxref[ix] = fx1[round(Int, idxfirst + (ix-1)*shiftidcs), fx_str]
    end

    #calculate normalization factor
    normfct = sumuptoscale ./ fxref

    #set limits to prohibit outliers
    normfct[normfct .> 5.0] .= 5.0
    normfct[normfct .< -5.0] .= -5.0

    #apply normalization
    for ix in 1:size(mrd_data, 2)
        mrd_data[:,ix] ./= normfct[ix]
    end
end

if exportnormed && norm
    @info("Writing to file...")
    MRD.write2dmrddatatonetcdf(x, mrd_data, time_middle, meanwind, evalstart, evalend, cols_for_mrd,
    nrpoints, blocklen, shift, sourcefile, string(readfilename[1:end-3], "_norm.nc"))
end

#-----------------------------------------------
#temporal y-scale

multfact = 1000
cmtouse = cramericm.vik
cmin = -5 #colorbar limits
cmax = -cmin
if "w" in cols_for_mrd
    if "u" in cols_for_mrd
        labeltext = L"$C_{uw} [10^{-3}~\mathrm{m^2~s^{-2}}]$"
        labelbottom = L"$\overline{uw} [\mathrm{m^2~s^{-2}}]$"
        flux = fx1.uw
        cmin = -5
        cmax = -cmin
        fluxscale = (-0.1, 0.2)
        rightscale = (-2.5,2)
    elseif "v" in cols_for_mrd
        labeltext = L"$C_{vw} [10^{-3}~\mathrm{m^2~s^{-2}}]$"
        labelbottom = L"$\overline{vw} [\mathrm{m^2~s^{-2}}]$"
        flux = fx1.vw
        cmin = -5
        cmax = -cmin
        fluxscale = (-0.1, 0.2)
        rightscale = (-2.5,2)
    elseif "T" in cols_for_mrd
        labeltext = L"$C_{wT_S} [10^{-3}~\mathrm{K~m~s^{-1}}]$"
        labelbottom = L"$\overline{wT_s} [\mathrm{K~m~s^{-1}}]$"
        flux = fx1.wT
        cmin = -5
        cmax = -cmin
        fluxscale = (-0.1, 0.2)
        rightscale = (-2.5,2)
    elseif cols_for_mrd == ["w", "w"]
        labeltext = L"$C_{ww} [10^{-3}~\mathrm{m^2}~s^{-2}]$"
        labelbottom =  L"$\overline{ww} [\mathrm{m^2}~s^{-2}]$"
        flux = fx1.ww
        cmtouse = "plasma"
        cmin = 0
        cmax = 10
        fluxscale = (0, 0.3)
        rightscale = (0,10)
    else
        @warn("neither 'u' nor 'v' nor 'T' in cols_for_mrd!")
    end
elseif cols_for_mrd == ["T"]
    labeltext = L"$C_{T} [10^{-3}~\mathrm{K}]$"
    labelbottom = L"$\overline{wT_s} [\mathrm{K~m~s^{-1}}]$"
    flux = fx1.wT
    cmin = -5
    cmax = -cmin
    fluxscale = (-0.1, 0.2)
    rightscale = (-2.5,2)
elseif cols_for_mrd == ["T", "T"]
    labeltext = L"$C_{T_s^2} [10^{-2}~\mathrm{K^2}]$"
    labelbottom = L"$\overline{(T_s)^2} [\mathrm{K^{2}}]$"
    flux = fx1.TT
    multfact=100
    cmtouse = cramericm.batlow
    cmin=0
    cmax=6
    fluxscale = (-0.1, 0.2)
    rightscale = (-2.5,2)
elseif cols_for_mrd == ["u", "u"]
    labeltext = L"$C_{uu} [10^{-3}~\mathrm{m^2}~s^{-2}]$"
    labelbottom =  L"$\overline{uu} [\mathrm{m^2}~s^{-2}]$"
    flux = fx1.uu
    cmtouse = "plasma"
    cmin = 0
    cmax = 50
    fluxscale = (0, 1)
    rightscale = (0,50)
elseif cols_for_mrd == ["tke"]
    @error("use part of script further below to plot 2D-MRD of TKE")
else
    @error("no w in cols_for_mrd!")
end

mrd_avg = zeros(Float64, size(mrd_data, 1))
mrd_q1 = zeros(Float64, size(mrd_data, 1))
mrd_q3 = zeros(Float64, size(mrd_data, 1))
x_avg = Vector{Millisecond}(undef, size(mrd_data, 1))

for i in axes(x, 1)
    (mrd_q1[i], mrd_avg[i], mrd_q3[i]) = quantile(filter(!isnan, mrd_data[i, :]), [0.25, 0.5, 0.75])
    x_avg[i] = Millisecond.(median(filter(x -> (x .!= 999999 && x .!= -9999), Dates.value.(x[i, :]))))
end

firstbiggerra = findfirst(x -> x > ra, Millisecond.(Dates.value.(x_avg)))

#plot temporal scale (y)

fig = PyPlot.figure(figsize=(12, 7.5))
gs = gridspec.GridSpec(3, 2, width_ratios=[7, 1], height_ratios=[1, 8, 1], hspace=0.05, wspace=0.05)
ax4 = fig.add_subplot(gs[1, 1])
ax4.set_title(figtitle)
ax4.plot(fluxdf.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
yl = ax4.get_ylim()
ax4.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="Set1", alpha=0.5, vmin=1, vmax=9)
ax4.set_ylabel(L"| \mathbf{u}|~\mathrm{[m~s^{-1}]}")
ax4.tick_params(axis="x", labelbottom=false)
ax4.grid()
ax4.set_ylim(0,8)
ax = fig.add_subplot(gs[2, 1], sharex=ax4)
mrdplot = ax.pcolormesh(time_middle, Dates.value.(x) ./ 10^3, mrd_data .* multfact, cmap=cmtouse, vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax.set_ylabel("time scale [s]")
ax.xaxis_date()
ax.set_ylim([0.1, 2e3])
ax.axhline(Dates.value(Millisecond(ra))/1e3, ls="--", color="black", alpha=0.3)
ax.fill_betweenx((Dates.value(Millisecond(ra))/1e3, last(ax.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(gs[2, 2], sharey=ax)
ax2.set_xlabel(labeltext)
ax2.grid()
(mrd1d,) = ax2.plot(mrd_avg .* 1000, Dates.value.(x_avg) ./ 1000)
xlim = ax2.get_xlim()
mrd1d_q = ax2.fill_betweenx(Dates.value.(x_avg) ./ 1000, mrd_q1 .* 1000, mrd_q3 .* 1000, alpha=0.3)
(mrd_white,) = ax2.plot(mrd_avg[firstbiggerra:end] .* 1000, Dates.value.(x_avg[firstbiggerra:end]) ./ 1000, color="white", alpha=0.5)
ax2.fill_betweenx(Dates.value.(x_avg[firstbiggerra:end]) ./ 1000, -100, 100, color="white", alpha=0.3)
ax2.axhline(Dates.value(ra), ls="--", color="black", alpha=0.3)
ax2.tick_params(axis="y", labelleft=false)
ax2.set_xlim(rightscale)
cbar = fig.colorbar(mrdplot, ax=ax2)#, orientation = "horizontal")
cbar.set_label(labeltext)
ax3 = fig.add_subplot(gs[3, 1], sharex=ax4)
#ax3.axhline(color="black", alpha=0.6)
ax3.plot(fx1.time[1:1200:end], flux[1:1200:end], color="black")
ax3.set_ylim(fluxscale)
ax3.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="Set1", alpha=0.5, vmin=1, vmax=9)
ax3.set_ylabel(labelbottom)
ax3.grid()
##
#=
#-----------------------------------------------
#spatial y-scale

x_spatial = fill(NaN, size(x, 1), size(x, 2))
meanwind_tmp = replace(meanwind, NaN => 1.0)
for ix in axes(x_spatial, 1)
    x_spatial[ix, :] = (Dates.value.(x[ix, :]) ./ 10^3) .* meanwind_tmp
end

ra_spatial_min = Dates.value(Second(ra)) .* minimum(meanwind_tmp)
firstbiggerra_sp = findfirst(x -> x >= ra_spatial_min, x_spatial)[1]

min_x_spatial = minimum(filter(!isnan, x_spatial))
max_x_spatial = maximum(filter(!isnan, x_spatial))
x_spatial_avg = collect(10 .^ range(log10(min_x_spatial), log10(max_x_spatial), nrpoints))

mrd_spatial_avg = zeros(Float64, size(x_spatial_avg, 1))
mrd_spatial_q1 = zeros(Float64, size(x_spatial_avg, 1))
mrd_spatial_q3 = zeros(Float64, size(x_spatial_avg, 1))
(mrd_spatial_q1[1], mrd_spatial_avg[1], mrd_spatial_q3[1]) =
    quantile(filter(!isnan, mrd_data[x_spatial.==min_x_spatial]), [0.25, 0.5, 0.75])

for i in axes(x_spatial_avg, 1)[2:end]
    try
        (mrd_spatial_q1[i], mrd_spatial_avg[i], mrd_spatial_q3[i]) =
            quantile(filter(!isnan, mrd_data[x_spatial_avg[i-1].<x_spatial.<=x_spatial_avg[i]]), [0.25, 0.5, 0.75])
    catch e
    end
end

#plot
cmin = -5 #colorbar limits
cmax = -cmin
fig = PyPlot.figure(figsize=(12, 7.5))
gs = gridspec.GridSpec(3, 2, width_ratios=[7, 1], height_ratios=[1, 8, 1], hspace=0.05, wspace=0.05)
ax4 = fig.add_subplot(gs[1, 1])
ax4.set_title(figtitle)
ax4.plot(fluxdf.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
ax4.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="Set1", alpha=0.5, vmin=1, vmax=9)
ax4.set_ylabel(L"\sqrt{u^2+v^2+w^2}~\mathrm{[m~s^{-1}]}")
ax4.tick_params(axis="x", labelbottom=false)
ax4.grid()
ax = fig.add_subplot(gs[2, 1], sharex=ax4)
mrdplot = ax.pcolormesh(time_middle, x_spatial, mrd_data .* 1000, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax.set_ylabel("avg. scale [m]")
ax.xaxis_date()
ax.set_ylim([0.1, 3e4])
ax.axhline(ra_spatial_min, ls="--", color="black", alpha=0.3)
ax.fill_betweenx((ra_spatial_min, last(ax.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(gs[2, 2], sharey=ax)
ax2.set_xlabel(L"$C_{w'T'} [\cdot 10^{-3}~\mathrm{K~m~s^{-1}}]$")
ax2.grid()
ax2.plot(mrd_spatial_avg .* 1000, x_spatial_avg)
xlim = ax2.get_xlim()
ax2.fill_betweenx(x_spatial_avg, mrd_spatial_q1 .* 1000, mrd_spatial_q3 .* 1000, alpha=0.3)
ax2.plot(mrd_spatial_avg[firstbiggerra_sp:end] .* 1000, x_spatial_avg[firstbiggerra_sp:end], color="white", alpha=0.5)
ax2.fill_betweenx(x_spatial_avg[firstbiggerra_sp:end], mrd_spatial_q1[firstbiggerra_sp:end] .* 1000, mrd_spatial_q3[firstbiggerra_sp:end] .* 1000, color="white", alpha=0.3)
ax2.axhline(ra_spatial_min, ls="--", color="black", alpha=0.3)
ax2.tick_params(axis="y", labelleft=false)
ax2.set_xlim(xlim)
cbar = fig.colorbar(mrdplot, ax=ax2)#, orientation = "horizontal")
cbar.set_label(L"$C_{w'T'} [\cdot 10^{-3}~\mathrm{K~m~s^{-1}}]$")
ax3 = fig.add_subplot(gs[3, 1], sharex=ax4)
#ax3.axhline(color="black", alpha=0.6)
ax3.plot(fx1.time[1:1200:end], fx1.wT[1:1200:end], color="black")
ax3.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="Set1", alpha=0.5, vmin=1, vmax=9)
ax3.set_ylabel(L"$\overline{w'T'}~\mathrm{[\mathrm{K~m~s^{-1}}]}$")
ax3.grid()
=#
######################################################
###        TAKE (SUB-)INTERVAL FROM 2D-MRD         ###
######################################################
timemiddle_start = DateTime(2021, 04, 01, 00, 00, 00)
timemiddle_end = DateTime(2021, 07, 01, 00, 00, 00)

mrd_data_sub = copy(mrd_data[:, timemiddle_start.<=time_middle.<=timemiddle_end])
time_middle_sub = copy(time_middle[timemiddle_start.<=time_middle.<=timemiddle_end])
x_sub = copy(x[:, timemiddle_start.<=time_middle.<=timemiddle_end])

######################################################
###   CATEGORIZE WT-DECOMP IN STABILITY CLASSES    ###
######################################################

#----------------------------------------------------
#take only specific wind directions if wanted
disallowmissing!(tjkdata)
windclass_2dmrd = fill(NaN, size(time_middle_sub, 1))
turb.extendtofinertimeseries!(windclass_2dmrd, time_middle_sub, windclass, tjkdata.time)

condition = windclass_2dmrd .== 1
mrd_data_sub[:, .!condition] .= NaN

#----------------------------------------------------

mrd_data_todo = mrd_data_sub

stabclass = zeros(Float64, size(mrd_data_todo, 2))
sumupto = Second(50)
sumuptoidx = findfirst(x -> x > sumupto, x[:, 1])
for i in 1:size(mrd_data_todo, 2)
    stabclass[i] = sum(mrd_data_todo[1:sumuptoidx, i])
end

occurpos = count(x -> x .>= 0, stabclass)/size(mrd_data_todo, 2)
occurneg = count(x -> x .< 0, stabclass)/size(mrd_data_todo, 2)

posmrd_data = mrd_data_sub[:, stabclass.>0]
negmrd_data = mrd_data_sub[:, stabclass.<=0]

posmrd_med = zeros(Float64, size(posmrd_data, 1))
posmrd_q1 = similar(posmrd_med)
posmrd_q3 = similar(posmrd_med)
negmrd_med = zeros(Float64, size(negmrd_data, 1))
negmrd_q1 = similar(negmrd_med)
negmrd_q3 = similar(negmrd_med)

for i in 1:size(posmrd_data, 1)
    (posmrd_q1[i], posmrd_med[i], posmrd_q3[i]) = quantile(filter(!isnan, posmrd_data[i, :]), [0.25, 0.5, 0.75])
    (negmrd_q1[i], negmrd_med[i], negmrd_q3[i]) = quantile(filter(!isnan, negmrd_data[i, :]), [0.25, 0.5, 0.75])
end

fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.set_title("T2LCSAT down valley wind")
ax.plot(Dates.value.(x[:, 1]) ./ 1000, posmrd_med .* 1e3, color="red", label="unstable+neutral")
ax.fill_between(Dates.value.(x[:, 1]) ./ 1000, posmrd_q1 .* 1e3, posmrd_q3 .* 1e3, color="red", alpha=0.5)
ax.plot(Dates.value.(x[:, 1]) ./ 1000, negmrd_med .* 1e3, color="blue", label="stable")
ax.fill_between(Dates.value.(x[:, 1]) ./ 1000, negmrd_q1 .* 1e3, negmrd_q3 .* 1e3, color="blue", alpha=0.5)
PyPlot.xscale("log")
ax.grid()
ax.set_xlabel("time scale [s]")
ax.set_ylabel(L"C_{wT_s} [10^{-3}~\mathrm{K~m~s^{-1}}]")
#ax.legend()

@show occurneg
@show occurpos

######################################################
###   CATEGORIZE TKE-DECOMP IN WIND DIR CLASSES    ###
######################################################
timemiddle_start = DateTime(2021, 05, 30, 00, 00, 00)
timemiddle_end = DateTime(2021, 07, 01, 00, 00, 00)

mrd_data_up = copy(mrd_data[:, timemiddle_start.<=time_middle.<=timemiddle_end])
mrd_data_down = copy(mrd_data[:, timemiddle_start.<=time_middle.<=timemiddle_end])
time_middle_sub = copy(time_middle[timemiddle_start.<=time_middle.<=timemiddle_end])
x_sub = copy(x[:, timemiddle_start.<=time_middle.<=timemiddle_end])

#----------------------------------------------------
#take only specific wind directions if wanted
disallowmissing!(tjkdata)
windclass_2dmrd = fill(NaN, size(time_middle_sub, 1))
turb.extendtofinertimeseries!(windclass_2dmrd, time_middle_sub, windclass, tjkdata.time)

upwind = windclass_2dmrd .== 1   #1 = up valley wind; 2 = down valley wind
downwind = windclass_2dmrd .== 2   #1 = up valley wind; 2 = down valley wind
mrd_data_up[:, .!upwind] .= NaN
mrd_data_down[:, .!downwind] .= NaN

#----------------------------------------------------

lin_mrd_up_med = zeros(Float64, size(mrd_data_up, 1))
lin_mrd_up_q1 = similar(lin_mrd_up_med)
lin_mrd_up_q3 = similar(lin_mrd_up_med)
lin_mrd_down_med = zeros(Float64, size(mrd_data_down, 1))
lin_mrd_down_q1 = similar(lin_mrd_down_med)
lin_mrd_down_q3 = similar(lin_mrd_down_med)

for i in 1:size(mrd_data_up, 1)
    (lin_mrd_up_q1[i], lin_mrd_up_med[i], lin_mrd_up_q3[i]) = quantile(filter(!isnan, mrd_data_up[i, :]), [0.25, 0.5, 0.75])
    (lin_mrd_down_q1[i], lin_mrd_down_med[i], lin_mrd_down_q3[i]) = quantile(filter(!isnan, mrd_data_down[i, :]), [0.25, 0.5, 0.75])
end

fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.set_title("T2IRG")
ax.plot(Dates.value.(x[:, 1]) ./ 1000, lin_mrd_down_med .* 1e2, color="blue")
ax.fill_between(Dates.value.(x[:, 1]) ./ 1000, lin_mrd_down_q1 .* 1e2, lin_mrd_down_q3 .* 1e2, color="blue", alpha=0.5)
ax.plot(Dates.value.(x[:, 1]) ./ 1000, lin_mrd_up_med .* 1e2, color="red")
ax.fill_between(Dates.value.(x[:, 1]) ./ 1000, lin_mrd_up_q1 .* 1e2, lin_mrd_up_q3 .* 1e2, color="red", alpha=0.5)
PyPlot.xscale("log")
ax.grid()
ax.set_xlabel("time scale [s]")
ax.set_ylabel(L"$C_{e}~[10^{-2}~\mathrm{J~kg^{-1}}]$")
#ax.legend()

######################################################
###               COMPARE 2DMRDs                   ###
######################################################
#=
(cmp1_time_middle, cmp1_x, cmp1_mrd_data, cmp1_meanwind,
cmp1_evalstart, cmp1_evalend, cmp1_cols_for_mrd,
cmp1_nrpoints, cmp1_blocklen, cmp1_shift, cmp1_sourcefile) = 
MRD.read2dmrdfromnetcdf(cmpfile1)

(cmp2_time_middle, cmp2_x, cmp2_mrd_data, cmp2_meanwind,
cmp2_evalstart, cmp2_evalend, cmp2_cols_for_mrd,
cmp2_nrpoints, cmp2_blocklen, cmp2_shift, cmp2_sourcefile) = 
MRD.read2dmrdfromnetcdf(cmpfile2)

#calculate factor
cmp_mrd = sign.(cmp2_mrd_data ./ cmp1_mrd_data)

inst1 = MRD.instrumentnamefromsourcefile(cmp1_sourcefile)
inst2 = MRD.instrumentnamefromsourcefile(cmp2_sourcefile)

#plot
cmin = -0.5 #colorbar limits
cmax = -cmin
fig = PyPlot.figure()
gs = gridspec.GridSpec(1, 1)#, width_ratios=[5, 1])
ax = fig.add_subplot(gs[1, 1])
ax.set_title(string("sign.(", inst2, " ./ ", inst1, ")"))
mrdplot = ax.pcolormesh(cmp1_time_middle, Dates.value.(cmp1_x) ./ 10^3, cmp_mrd, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax.set_ylabel("avg. time [s]")
ax.xaxis_date()
ax.set_ylim([0.1, 1.5e4])
cbar = fig.colorbar(mrdplot, ax=ax)#, orientation = "horizontal")
cbar.set_label(L"$C_{w'T'} ./ C_{w'T'}$")
#axgb.tick_params(axis="y", labelleft=false)
=#
######################################################
###               4 PANEL 2DMRDs                   ###
######################################################
fluxtype = "TT"
fluxfilename = fluxtype

if fluxtype == "wT" || fluxtype == "tke"
    fluxfilename = string(fluxfilename, "_norm.nc")
else
    fluxfilename = string(fluxfilename, "_norm.nc")
end

panelfile1 = joinpath(datapath, "2dmrd", "db_tot", string("tjk_", fluxfilename))
panelfile2 = joinpath(datapath, "2dmrd", "db_tot", string("t2ucsat_", fluxfilename))
panelfile3 = joinpath(datapath, "2dmrd", "db_tot", string("t2lcsat_", fluxfilename))
panelfile4 = joinpath(datapath, "2dmrd", "db_tot", string("t2irg_", fluxfilename))

#=
panelfile1 = joinpath(datapath, "2dmrd", "db_tot", string("t2ucsat_", fluxtype, "_norm.nc"))
panelfile2 = joinpath(datapath, "2dmrd", "db_tot", string("t2lcsat_", fluxtype, "_norm.nc"))
panelfile3 = joinpath(datapath, "2dmrd", "db_tot", string("t2irg_", fluxtype, "_norm.nc"))
panelfile4 = joinpath(datapath, "2dmrd", "db_tot", string("kaijo_", fluxtype, "_norm.nc"))
=#

(panel1_time_middle, panel1_x, panel1_mrd_data, panel1_meanwind,
    panel1_evalstart, panel1_evalend, panel1_cols_for_mrd,
    panel1_nrpoints, panel1_blocklen, panel1_shift, panel1_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile1)

(panel2_time_middle, panel2_x, panel2_mrd_data, panel2_meanwind,
    panel2_evalstart, panel2_evalend, panel2_cols_for_mrd,
    panel2_nrpoints, panel2_blocklen, panel2_shift, panel2_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile2)

(panel3_time_middle, panel3_x, panel3_mrd_data, panel3_meanwind,
    panel3_evalstart, panel3_evalend, panel3_cols_for_mrd,
    panel3_nrpoints, panel3_blocklen, panel3_shift, panel3_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile3)

(panel4_time_middle, panel4_x, panel4_mrd_data, panel4_meanwind,
    panel4_evalstart, panel4_evalend, panel4_cols_for_mrd,
    panel4_nrpoints, panel4_blocklen, panel4_shift, panel4_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile4)

panelinst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
panelinst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
panelinst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
panelinst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

panel1ra = inst2ra[panelinst1]
panel2ra = inst2ra[panelinst2]
panel3ra = inst2ra[panelinst3]
panel4ra = inst2ra[panelinst4]

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[panel1_evalstart.<=tjkmeteo[:, 1].<=panel2_evalend, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 1
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 2

#load data and calc fluxes
inst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
inst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
inst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
inst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)
#=
panel1_x_avg = Vector{Millisecond}(undef, size(panel1_mrd_data, 1))
panel2_x_avg = Vector{Millisecond}(undef, size(panel2_mrd_data, 1))
panel3_x_avg = Vector{Millisecond}(undef, size(panel3_mrd_data, 1))
panel4_x_avg = Vector{Millisecond}(undef, size(panel4_mrd_data, 1))

for i in axes(panel1_x, 1)
    panel1_x_avg[i] = Millisecond.(median(filter(x -> (x .!= 999999 && x .!= -9999), Dates.value.(panel1_x[i, :]))))
end
for i in axes(panel2_x, 1)
    panel2_x_avg[i] = Millisecond.(median(filter(x -> (x .!= 999999 && x .!= -9999), Dates.value.(panel2_x[i, :]))))
end
for i in axes(panel3_x, 1)
    panel3_x_avg[i] = Millisecond.(median(filter(x -> (x .!= 999999 && x .!= -9999), Dates.value.(panel1_x[i, :]))))
end
for i in axes(panel4_x, 1)
    panel4_x_avg[i] = Millisecond.(median(filter(x -> (x .!= 999999 && x .!= -9999), Dates.value.(panel1_x[i, :]))))
end

firstbiggerra1 = size(panel1_x, 1)#findfirst(x -> x > panel1ra, Millisecond.(Dates.value.(panel1_x_avg)))
firstbiggerra2 = size(panel2_x, 1)#findfirst(x -> x > panel2ra, Millisecond.(Dates.value.(panel2_x_avg)))
firstbiggerra3 = size(panel3_x, 1)#findfirst(x -> x > panel3ra, Millisecond.(Dates.value.(panel3_x_avg)))
firstbiggerra4 = size(panel4_x, 1)#findfirst(x -> x > panel4ra, Millisecond.(Dates.value.(panel4_x_avg)))

panel1_sum = zeros(Float64, size(panel1_time_middle, 1))
panel2_sum = zeros(Float64, size(panel2_time_middle, 1))
panel3_sum = zeros(Float64, size(panel3_time_middle, 1))
panel4_sum = zeros(Float64, size(panel4_time_middle, 1))

for i in axes(panel1_sum, 1)
    panel1_sum[i] = sum(filter(!isnan, panel1_mrd_data[1:firstbiggerra1, i]))
end
for i in axes(panel2_sum, 1)
    panel2_sum[i] = sum(filter(!isnan, panel2_mrd_data[1:firstbiggerra2, i]))
end
for i in axes(panel3_sum, 1)
    panel3_sum[i] = sum(filter(!isnan, panel3_mrd_data[1:firstbiggerra3, i]))
end
for i in axes(panel4_sum, 1)
    panel4_sum[i] = sum(filter(!isnan, panel4_mrd_data[1:firstbiggerra4, i]))
end
=#
#fluxfile1 = inst2file[inst1]
fluxfile2 = inst2file[inst2]
#fluxfile3 = inst2file[inst3]
#fluxfile4 = inst2file[inst4]
#ra1 = inst2ra[inst1]
ra2 = inst2ra[inst2]
#ra3 = inst2ra[inst3]
#ra4 = inst2ra[inst4]

fx2df = turb.readturbasnetcdf(fluxfile2, panel2_evalstart, panel2_evalend)
#show statistics about missing data
turb.printmissstats(fx2df)
#interpolate missing
fx2df = turb.interpolatemissing(fx2df)
#double rotation
if inst2 == "TJK"
    turb.drdf!(fx2, periodwise=false)
else
    turb.drdf!(fx2df)
end
turb.missing2nan!(fx2df)
wndspd = gen.movingaverage(sqrt.(fx2df.u .^ 2 .+ fx2df.v .^ 2 .+ fx2df.w .^ 2), 20 * 3600)

fluxpath = "/home/haugened/Documents/data/fluxes/"
fx1file = joinpath(fluxpath, string(lowercase(inst1), "_fx.nc"))
fx2file = joinpath(fluxpath, string(lowercase(inst2), "_fx.nc"))
fx3file = joinpath(fluxpath, string(lowercase(inst3), "_fx.nc"))
fx4file = joinpath(fluxpath, string(lowercase(inst4), "_fx.nc"))

fx1 = turb.readturbasnetcdf(fx1file, panel1_evalstart, panel1_evalend)
fx2 = turb.readturbasnetcdf(fx2file, panel2_evalstart, panel2_evalend)
fx3 = turb.readturbasnetcdf(fx3file, panel3_evalstart, panel3_evalend)
fx4 = turb.readturbasnetcdf(fx4file, panel4_evalstart, panel4_evalend)

multfactor = 1

if fluxtype == "wT"
    cmin = -5 #colorbar limits
    cmax = -cmin
    multfactor = 1000
elseif fluxtype == "tke"
    cmin = 0
    cmax = 15
    multfactor = 100
elseif fluxtype == "TT"
    cmin = 0
    cmax = 6
    multfactor = 100
else
    @error("fluxtype not correct")
end

##
#plot
fig = PyPlot.figure(figsize=(11,8))
gs = gridspec.GridSpec(6, 2, height_ratios=[5, 15, 15, 15, 15, 8], width_ratios=[50,1])
date_format = pydates.DateFormatter("%d.%m.")
majorloc = pydates.DayLocator(interval=1)
minorloc = pydates.HourLocator([6,12,18])
ax = fig.add_subplot(gs[1, 1])
ax.plot(fx2df.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
yl = ax.get_ylim()
ax.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="Set1", alpha=0.5, vmin=1, vmax=9)
ax.set_ylabel(L"| \mathbf{u}|_{3\mathrm{m}}~\mathrm{[m~s^{-1}]}")
ax.tick_params(axis="x", labelbottom=false)
ax.grid()
ax.set_ylim(0,5.5)
ax.set_xlim(DateTime(2021,05,29,00,00,00), DateTime(2021,06,07,00,00,00))
fig.text(0.825, 0.73, "5m", rotation=-90, fontweight="bold")
ax1 = fig.add_subplot(gs[2, 1], sharex=ax)
#ax1.set_title("5m meteo station")
mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* multfactor, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax1.xaxis_date()
ax1.set_ylim([0.1, 2e3])
ax1.tick_params(axis="x", labelbottom=false)
ax1.axhline(Dates.value(Millisecond(panel1ra))./1000, ls="--", color="black", alpha=0.3)
ax1.fill_betweenx((Dates.value(Millisecond(panel1ra))./1000, last(ax1.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.57, "3m", rotation=-90, fontweight="bold")
ax2 = fig.add_subplot(gs[3, 1], sharex=ax, sharey=ax1)
#ax2.set_title("3m T2")
mrdplot2 = ax2.pcolormesh(panel2_time_middle, Dates.value.(panel2_x) ./ 1000, panel2_mrd_data .* multfactor, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
ax2.tick_params(axis="x", labelbottom=false)
#ax2.tick_params(axis="y", labelleft=false)
ax2.axhline(Dates.value(Millisecond(panel2ra))./1000, ls="--", color="black", alpha=0.3)
ax2.fill_betweenx((Dates.value(Millisecond(panel2ra))./1000, last(ax2.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.42, "2m", rotation=-90, fontweight="bold")
ax3 = fig.add_subplot(gs[4, 1], sharex=ax, sharey=ax1)
#ax3.set_title("2m T2")
mrdplot3 = ax3.pcolormesh(panel3_time_middle, Dates.value.(panel3_x) ./ 1000, panel3_mrd_data .* multfactor, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
ax3.axhline(Dates.value(Millisecond(panel3ra))./1000, ls="--", color="black", alpha=0.3)
ax3.fill_betweenx((Dates.value(Millisecond(panel3ra))./1000, last(ax3.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax3.tick_params(axis="x", labelbottom=false)
ax4 = fig.add_subplot(gs[5, 1], sharex=ax, sharey=ax1)
#ax4.set_title("1m T2")
mrdplot4 = ax4.pcolormesh(panel4_time_middle, Dates.value.(panel4_x) ./ 1000, panel4_mrd_data .* multfactor, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
#ax4.tick_params(axis="y", labelleft=false)
ax4.axhline(Dates.value(Millisecond(panel4ra))./1000, ls="--", color="black", alpha=0.3)
ax4.fill_betweenx((Dates.value(Millisecond(panel4ra))./1000, last(ax4.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax4.tick_params(axis="x", labelbottom=false)
fig.text(0.825, 0.26, "1m", rotation=-90, fontweight="bold")
ax5 = fig.add_subplot(py"$(gs)[2:-2, 1]")
cbar = fig.colorbar(mrdplot, cax=ax5, orientation="vertical")
if fluxtype == "wT"
    cbar.set_label(L"$C_{wT_s} [10^{-3}~\mathrm{K~m~s^{-1}}]$")
elseif fluxtype == "tke"
    cbar.set_label(L"$C_{e}~[\mathrm{10^{-2}~m^2~s^{-2}}]$")
elseif fluxtype == "TT"
    cbar.set_label(L"$C_{T_s^2}~[\mathrm{10^{-2}~K^2}]$")
end
fig.text(0.05, 0.45, "time scale [s]", rotation="vertical")
ax6 = fig.add_subplot(gs[6, 1], sharex = ax)
if fluxtype == "wT"
    ax6.plot(fx2.time, fx2.wT, label="3m")
    ax6.plot(fx4.time, fx4.wT, alpha=0.9, label="1m")
    ax6.set_ylabel(L"$\overline{w'T_s'}~[\mathrm{K~m~s^{-1}}]$")
    ax6.set_ylim(-0.05, 0.2)
elseif fluxtype == "tke"
    ax6.plot(fx2.time, fx2.tke, label="3m")
    ax6.plot(fx4.time, fx4.tke, alpha=0.9, label="1m")
    ax6.set_ylabel(L"$e~\mathrm{[m^2~s^{-2}]}$")
    ax6.set_ylim(0,1.2)
elseif fluxtype == "TT"
    ax6.plot(fx2.time, fx2.TT, label="3m")
    ax6.plot(fx4.time, fx4.TT, label="1m", alpha=0.8)
    ax6.set_ylabel(L"$\overline{\left( T_s' \right)^2}~\mathrm{[K^2]}$")
    ax6.set_ylim(0,0.8)
end
ax6.legend(loc="center left", bbox_to_anchor=(1,0.5))
ax6.grid()
ax6.xaxis.set_major_locator(majorloc)
ax6.xaxis.set_minor_locator(minorloc)
ax6.xaxis.set_major_formatter(date_format)
##

######################################################
###          4 PANEL 2DMRDs WITH KAIJO             ###
######################################################
fluxtype = "TT"
fluxfilename = fluxtype

if fluxtype == "wT" || fluxtype == "tke"
    fluxfilename = string(fluxfilename, "_norm.nc")
else
    fluxfilename = string(fluxfilename, "_norm.nc")
end

panelfile1 = joinpath(datapath, "2dmrd", "db_tot", string("t2ucsat_", fluxfilename))
panelfile2 = joinpath(datapath, "2dmrd", "db_tot", string("t2lcsat_", fluxfilename))
panelfile3 = joinpath(datapath, "2dmrd", "db_tot", string("t2irg_", fluxfilename))
panelfile4 = joinpath(datapath, "2dmrd", "db_tot", string("kaijo_", fluxfilename))

(panel1_time_middle, panel1_x, panel1_mrd_data, panel1_meanwind,
    panel1_evalstart, panel1_evalend, panel1_cols_for_mrd,
    panel1_nrpoints, panel1_blocklen, panel1_shift, panel1_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile1)

(panel2_time_middle, panel2_x, panel2_mrd_data, panel2_meanwind,
    panel2_evalstart, panel2_evalend, panel2_cols_for_mrd,
    panel2_nrpoints, panel2_blocklen, panel2_shift, panel2_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile2)

(panel3_time_middle, panel3_x, panel3_mrd_data, panel3_meanwind,
    panel3_evalstart, panel3_evalend, panel3_cols_for_mrd,
    panel3_nrpoints, panel3_blocklen, panel3_shift, panel3_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile3)

(panel4_time_middle, panel4_x, panel4_mrd_data, panel4_meanwind,
    panel4_evalstart, panel4_evalend, panel4_cols_for_mrd,
    panel4_nrpoints, panel4_blocklen, panel4_shift, panel4_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile4)

panelinst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
panelinst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
panelinst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
panelinst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

panel1ra = inst2ra[panelinst1]
panel2ra = inst2ra[panelinst2]
panel3ra = inst2ra[panelinst3]
panel4ra = inst2ra[panelinst4]

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[panel1_evalstart.<=tjkmeteo[:, 1].<=panel2_evalend, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 1
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 2

#load data and calc fluxes
inst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
inst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
inst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
inst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)
#=
panel1_x_avg = Vector{Millisecond}(undef, size(panel1_mrd_data, 1))
panel2_x_avg = Vector{Millisecond}(undef, size(panel2_mrd_data, 1))
panel3_x_avg = Vector{Millisecond}(undef, size(panel3_mrd_data, 1))
panel4_x_avg = Vector{Millisecond}(undef, size(panel4_mrd_data, 1))

for i in axes(panel1_x, 1)
    panel1_x_avg[i] = Millisecond.(median(filter(x -> (x .!= 999999 && x .!= -9999), Dates.value.(panel1_x[i, :]))))
end
for i in axes(panel2_x, 1)
    panel2_x_avg[i] = Millisecond.(median(filter(x -> (x .!= 999999 && x .!= -9999), Dates.value.(panel2_x[i, :]))))
end
for i in axes(panel3_x, 1)
    panel3_x_avg[i] = Millisecond.(median(filter(x -> (x .!= 999999 && x .!= -9999), Dates.value.(panel1_x[i, :]))))
end
for i in axes(panel4_x, 1)
    panel4_x_avg[i] = Millisecond.(median(filter(x -> (x .!= 999999 && x .!= -9999), Dates.value.(panel1_x[i, :]))))
end

firstbiggerra1 = size(panel1_x, 1)#findfirst(x -> x > panel1ra, Millisecond.(Dates.value.(panel1_x_avg)))
firstbiggerra2 = size(panel2_x, 1)#findfirst(x -> x > panel2ra, Millisecond.(Dates.value.(panel2_x_avg)))
firstbiggerra3 = size(panel3_x, 1)#findfirst(x -> x > panel3ra, Millisecond.(Dates.value.(panel3_x_avg)))
firstbiggerra4 = size(panel4_x, 1)#findfirst(x -> x > panel4ra, Millisecond.(Dates.value.(panel4_x_avg)))

panel1_sum = zeros(Float64, size(panel1_time_middle, 1))
panel2_sum = zeros(Float64, size(panel2_time_middle, 1))
panel3_sum = zeros(Float64, size(panel3_time_middle, 1))
panel4_sum = zeros(Float64, size(panel4_time_middle, 1))

for i in axes(panel1_sum, 1)
    panel1_sum[i] = sum(filter(!isnan, panel1_mrd_data[1:firstbiggerra1, i]))
end
for i in axes(panel2_sum, 1)
    panel2_sum[i] = sum(filter(!isnan, panel2_mrd_data[1:firstbiggerra2, i]))
end
for i in axes(panel3_sum, 1)
    panel3_sum[i] = sum(filter(!isnan, panel3_mrd_data[1:firstbiggerra3, i]))
end
for i in axes(panel4_sum, 1)
    panel4_sum[i] = sum(filter(!isnan, panel4_mrd_data[1:firstbiggerra4, i]))
end
=#
#fluxfile1 = inst2file[inst1]
fluxfile2 = inst2file[inst2]
#fluxfile3 = inst2file[inst3]
#fluxfile4 = inst2file[inst4]
#ra1 = inst2ra[inst1]
ra2 = inst2ra[inst2]
#ra3 = inst2ra[inst3]
#ra4 = inst2ra[inst4]

fx2df = turb.readturbasnetcdf(fluxfile2, panel2_evalstart, panel2_evalend)
#show statistics about missing data
turb.printmissstats(fx2df)
#interpolate missing
fx2df = turb.interpolatemissing(fx2df)
#double rotation
if inst2 == "TJK"
    turb.drdf!(fx2, periodwise=false)
else
    turb.drdf!(fx2df)
end
turb.missing2nan!(fx2df)
wndspd = gen.movingaverage(sqrt.(fx2df.u .^ 2 .+ fx2df.v .^ 2 .+ fx2df.w .^ 2), 20 * 3600)

fluxpath = "/home/haugened/Documents/data/fluxes/"
fx1file = joinpath(fluxpath, string(lowercase(inst1), "_fx.nc"))
fx2file = joinpath(fluxpath, string(lowercase(inst2), "_fx.nc"))
fx3file = joinpath(fluxpath, string(lowercase(inst3), "_fx.nc"))
fx4file = joinpath(fluxpath, string(lowercase(inst4), "_fx.nc"))

fx1 = turb.readturbasnetcdf(fx1file, panel1_evalstart, panel1_evalend)
fx2 = turb.readturbasnetcdf(fx2file, panel2_evalstart, panel2_evalend)
fx3 = turb.readturbasnetcdf(fx3file, panel3_evalstart, panel3_evalend)
fx4 = turb.readturbasnetcdf(fx4file, panel4_evalstart, panel4_evalend)

multfactor = 1

if fluxtype == "wT"
    cmin = -5 #colorbar limits
    cmax = -cmin
    multfactor = 1000
elseif fluxtype == "tke"
    cmin = 0
    cmax = 15
    multfactor = 100
elseif fluxtype == "TT"
    cmin = 0
    cmax = 6
    multfactor = 100
else
    @error("fluxtype not correct")
end

##
#plot
fig = PyPlot.figure(figsize=(11,8))
gs = gridspec.GridSpec(6, 2, height_ratios=[5, 15, 15, 15, 15, 8], width_ratios=[50,1])
date_format = pydates.DateFormatter("%d.%m. %H:%M")
majorloc = pydates.HourLocator([0,6,12,18])
minorloc = pydates.HourLocator(interval=1)
ax = fig.add_subplot(gs[1, 1])
ax.plot(fx2df.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
yl = ax.get_ylim()
ax.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="Set1", alpha=0.5, vmin=1, vmax=9)
ax.set_ylabel(L"| \mathbf{u}|_{3\mathrm{m}}~\mathrm{[m~s^{-1}]}")
ax.tick_params(axis="x", labelbottom=false)
ax.grid()
ax.set_ylim(0,5.5)
ax.set_xlim(DateTime(2021,06,01,06,00,00), DateTime(2021,06,02,12,00,00))
fig.text(0.825, 0.73, "3m", rotation=-90, fontweight="bold")
ax1 = fig.add_subplot(gs[2, 1], sharex=ax)
#ax1.set_title("5m meteo station")
mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* multfactor, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax1.xaxis_date()
ax1.set_ylim([0.1, 2e3])
ax1.tick_params(axis="x", labelbottom=false)
ax1.axhline(Dates.value(Millisecond(panel1ra))./1000, ls="--", color="black", alpha=0.3)
ax1.fill_betweenx((Dates.value(Millisecond(panel1ra))./1000, last(ax1.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.57, "2m", rotation=-90, fontweight="bold")
ax2 = fig.add_subplot(gs[3, 1], sharex=ax, sharey=ax1)
#ax2.set_title("3m T2")
mrdplot2 = ax2.pcolormesh(panel2_time_middle, Dates.value.(panel2_x) ./ 1000, panel2_mrd_data .* multfactor, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
ax2.tick_params(axis="x", labelbottom=false)
#ax2.tick_params(axis="y", labelleft=false)
ax2.axhline(Dates.value(Millisecond(panel2ra))./1000, ls="--", color="black", alpha=0.3)
ax2.fill_betweenx((Dates.value(Millisecond(panel2ra))./1000, last(ax2.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.42, "1m", rotation=-90, fontweight="bold")
ax3 = fig.add_subplot(gs[4, 1], sharex=ax, sharey=ax1)
#ax3.set_title("2m T2")
mrdplot3 = ax3.pcolormesh(panel3_time_middle, Dates.value.(panel3_x) ./ 1000, panel3_mrd_data .* multfactor, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
ax3.axhline(Dates.value(Millisecond(panel3ra))./1000, ls="--", color="black", alpha=0.3)
ax3.fill_betweenx((Dates.value(Millisecond(panel3ra))./1000, last(ax3.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax3.tick_params(axis="x", labelbottom=false)
ax4 = fig.add_subplot(gs[5, 1], sharex=ax, sharey=ax1)
#ax4.set_title("1m T2")
mrdplot4 = ax4.pcolormesh(panel4_time_middle, Dates.value.(panel4_x) ./ 1000, panel4_mrd_data .* multfactor, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
#ax4.tick_params(axis="y", labelleft=false)
ax4.axhline(Dates.value(Millisecond(panel4ra))./1000, ls="--", color="black", alpha=0.3)
ax4.fill_betweenx((Dates.value(Millisecond(panel4ra))./1000, last(ax4.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax4.tick_params(axis="x", labelbottom=false)
fig.text(0.825, 0.26, "0.3m", rotation=-90, fontweight="bold")
ax5 = fig.add_subplot(py"$(gs)[2:-2, 1]")
cbar = fig.colorbar(mrdplot, cax=ax5, orientation="vertical")
if fluxtype == "wT"
    cbar.set_label(L"$C_{wT_s} [10^{-3}~\mathrm{K~m~s^{-1}}]$")
elseif fluxtype == "tke"
    cbar.set_label(L"$C_{e}~[\mathrm{10^{-2}~m^2~s^{-2}}]$")
elseif fluxtype == "TT"
    cbar.set_label(L"$C_{T_s^2}~[\mathrm{10^{-2}~K^2}]$")
end
fig.text(0.05, 0.45, "time scale [s]", rotation="vertical")
ax6 = fig.add_subplot(gs[6, 1], sharex = ax)
if fluxtype == "wT"
    ax6.plot(fx2.time, fx2.wT, label="2m")
    ax6.plot(fx4.time, fx4.wT, alpha=0.9, label="0.3m")
    ax6.set_ylabel(L"$\overline{w'T_s'}~[\mathrm{K~m~s^{-1}}]$")
    ax6.set_ylim(-0.05, 0.2)
elseif fluxtype == "tke"
    ax6.plot(fx2.time, fx2.tke, label="2m")
    ax6.plot(fx4.time, fx4.tke, alpha=0.9, label="0.3m")
    ax6.set_ylabel(L"$e~\mathrm{[m^2~s^{-2}]}$")
    ax6.set_ylim(0,1.2)
elseif fluxtype == "TT"
    ax6.plot(fx2.time, fx2.TT, label="2m")
    ax6.plot(fx4.time, fx4.TT, alpha=0.9, label="0.3m")
    ax6.set_ylabel(L"$\overline{\left( T_s' \right)^2}~\mathrm{[K^2]}$")
    ax6.set_ylim(0,0.8)
end
ax6.legend(loc="center left", bbox_to_anchor=(1,0.5))
ax6.grid()
ax6.xaxis.set_major_locator(majorloc)
ax6.xaxis.set_minor_locator(minorloc)
ax6.xaxis.set_major_formatter(date_format)
##
######################################################
###                    TKE PLOT                    ###
######################################################
#temporal y-scale

mrd_avg = zeros(Float64, size(mrd_data, 1))
mrd_q1 = zeros(Float64, size(mrd_data, 1))
mrd_q3 = zeros(Float64, size(mrd_data, 1))
x_avg = Vector{Millisecond}(undef, size(mrd_data, 1))

for i in axes(x, 1)
    (mrd_q1[i], mrd_avg[i], mrd_q3[i]) = quantile(filter(!isnan, mrd_data[i, :]), [0.25, 0.5, 0.75])
    x_avg[i] = Millisecond.(median(filter(x -> (x .!= 999999 && x .!= -9999), Dates.value.(x[i, :]))))
end

firstbiggerra = findfirst(x -> x > ra, Millisecond.(Dates.value.(x_avg)))

tke_total = zeros(Float64, size(mrd_data, 2))
for i in 1:length(tke_total)
    tke_total[i] = sum(mrd_data[1:firstbiggerra-1, i])
end

#plot temporal scale (y)
cmin = 0 #colorbar limits
cmax = 15
fig = PyPlot.figure(figsize=(12, 7.5))
gs = gridspec.GridSpec(3, 2, width_ratios=[7, 1], height_ratios=[1, 8, 1], hspace=0.05, wspace=0.05)
ax4 = fig.add_subplot(gs[1, 1])
ax4.set_title(figtitle)
#ax4.plot(tjkdata.time, tjkdata.wind_mean_vector_direction, color="black")
ax4.plot(fluxdf.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
ax4.set_ylim(0, 6)
yl = ax4.get_ylim()
ax4.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="Set1", alpha=0.5, vmin=1, vmax=9)
ax4.set_ylabel(L"| \mathbf{u}|~\mathrm{[m~s^{-1}]}")
#ax4.set_ylabel("winddir [°]")
ax4.tick_params(axis="x", labelbottom=false)
ax4.grid()
ax4.set_xlim(DateTime(2021,05,27,00,00,00), DateTime(2021,06,11,00,00,00))
ax = fig.add_subplot(gs[2, 1], sharex=ax4)
mrdplot = ax.pcolormesh(time_middle, Dates.value.(x) ./ 10^3, mrd_data .*100, cmap=cramericm.batlow, vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax.set_ylabel("time scale [s]")
ax.xaxis_date()
ax.set_ylim([0.1, 2e3])
ax.axhline(Dates.value(Millisecond(ra))./1e3, ls="--", color="black", alpha=0.3)
ax.fill_betweenx((Dates.value(Millisecond(ra))./1e3, last(ax.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(gs[2, 2], sharey=ax)
ax2.set_xlabel(L"$C_{e} [\mathrm{J~kg^{-1}}]$")
ax2.grid()
(mrd1d,) = ax2.plot(mrd_avg, Dates.value.(x_avg) ./ 1000)
mrd1d_q = ax2.fill_betweenx(Dates.value.(x_avg) ./ 1000, mrd_q1, mrd_q3, alpha=0.3)
xlim = ax2.get_xlim()
(mrd_white,) = ax2.plot(mrd_avg[firstbiggerra:end], Dates.value.(x_avg[firstbiggerra:end]) ./ 1000, color="white", alpha=0.5)
ax2.fill_betweenx(Dates.value.(x_avg[firstbiggerra:end]) ./ 1000, -100, 100, color="white", alpha=0.3)
ax2.axhline(Dates.value(ra), ls="--", color="black", alpha=0.3)
ax2.tick_params(axis="y", labelleft=false)
ax2.set_xlim(xlim)
cbar = fig.colorbar(mrdplot, ax=ax2)#, orientation = "horizontal")
cbar.set_label(L"$C_{e} [10^{-2}~\mathrm{J~kg^{-1}}]$")
ax3 = fig.add_subplot(gs[3, 1], sharex=ax4)
#ax3.axhline(color="black", alpha=0.6)
#ax3.plot(pydates.date2num.(time_middle), tke_total, color="black")
ax3.plot(fx1.time[1:1200:end], fx1.tke[1:1200:end], color="black")
#ax3.set_ylim((0, 30))
ax3.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(0, 1.5 * maximum(filter(!isnan, tke_total)), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="Set1", alpha=0.5, vmin=1, vmax=9)
ax3.set_ylabel(L"$e~\mathrm{[J~kg^{-1}]}$")
ax3.grid()
ax3.set_ylim(0,1.5)
#axbt = fig.add_axes([0.8,0.05,0.1,0.075])
#bupd = widgets.Button(axbt, "Update 1D-MRD")
#mrd1d_q = bupd.on_clicked(update_1dmrd)

######################################################
###                 U STAR PLOT                    ###
######################################################
#temporal y-scale
#=
mrd_avg = zeros(Float64, size(mrd_data, 1))
mrd_q1 = zeros(Float64, size(mrd_data, 1))
mrd_q3 = zeros(Float64, size(mrd_data, 1))
x_avg = Vector{Millisecond}(undef, size(mrd_data, 1))

for i in axes(x, 1)
    (mrd_q1[i], mrd_avg[i], mrd_q3[i]) = quantile(filter(!isnan, mrd_data[i, :]), [0.25, 0.5, 0.75])
    x_avg[i] = Millisecond.(median(filter(x -> (x .!= 999999 && x .!= -9999), Dates.value.(x[i, :]))))
end

firstbiggerra = findfirst(x -> x > ra, Millisecond.(Dates.value.(x_avg)))

ustar_total = zeros(Float64, size(mrd_data, 2))
for i in 1:length(ustar_total)
    ustar_total[i] = sum(mrd_data[1:firstbiggerra-1, i])
end

#plot temporal scale (y)
cmin = 0 #colorbar limits
cmax = 0.3
fig = PyPlot.figure(figsize=(12, 7.5))
gs = gridspec.GridSpec(3, 2, width_ratios=[7, 1], height_ratios=[1, 8, 1], hspace=0.05, wspace=0.05)
ax4 = fig.add_subplot(gs[1, 1])
ax4.set_title(figtitle)
ax4.plot(fluxdf.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
ax4.set_ylim(0, 6)
yl = ax4.get_ylim()
ax4.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="Set1", alpha=0.5, vmin=1, vmax=9)
ax4.set_ylabel(L"\sqrt{u^2+v^2+w^2}~\mathrm{[m~s^{-1}]}")
ax4.tick_params(axis="x", labelbottom=false)
ax4.grid()
ax = fig.add_subplot(gs[2, 1], sharex=ax4)
mrdplot = ax.pcolormesh(time_middle, Dates.value.(x) ./ 10^3, mrd_data, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax.set_ylabel("avg. time [s]")
ax.xaxis_date()
ax.set_ylim([0.1, 2e3])
ax.axhline(Dates.value(Second(ra)), ls="--", color="black", alpha=0.3)
ax.fill_betweenx((Dates.value(Second(ra)), last(ax.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(gs[2, 2], sharey=ax)
ax2.set_xlabel(L"$C_{u_\ast} [\mathrm{m~s^{-1}}]$")
ax2.grid()
(mrd1d,) = ax2.plot(mrd_avg, Dates.value.(x_avg) ./ 1000)
mrd1d_q = ax2.fill_betweenx(Dates.value.(x_avg) ./ 1000, mrd_q1, mrd_q3, alpha=0.3)
xlim = ax2.get_xlim()
(mrd_white,) = ax2.plot(mrd_avg[firstbiggerra:end], Dates.value.(x_avg[firstbiggerra:end]) ./ 1000, color="white", alpha=0.5)
ax2.fill_betweenx(Dates.value.(x_avg[firstbiggerra:end]) ./ 1000, -100, 100, color="white", alpha=0.3)
ax2.axhline(Dates.value(ra), ls="--", color="black", alpha=0.3)
ax2.tick_params(axis="y", labelleft=false)
ax2.set_xlim(xlim)
cbar = fig.colorbar(mrdplot, ax=ax2)#, orientation = "horizontal")
cbar.set_label(L"$C_{u_\ast} [\mathrm{m~s^{-1}}]$")
ax3 = fig.add_subplot(gs[3, 1], sharex=ax4)
#ax3.axhline(color="black", alpha=0.6)
#ax3.plot(pydates.date2num.(time_middle), ustar_total ./ 40, color="green")
ax3.plot(fx1.time[1:1200:end], fx1.u_star[1:1200:end], color="black")
x3ylim = ax3.get_ylim()
ax3.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(0, 1.5 * maximum(filter(!isnan, ustar_total)), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="Set1", alpha=0.5, vmin=1, vmax=9)
ax3.set_ylabel(L"$u_\ast~\mathrm{[m~s^{-1}]}$")
ax3.grid()
ax3.set_ylim(x3ylim)
#axbt = fig.add_axes([0.8,0.05,0.1,0.075])
#bupd = widgets.Button(axbt, "Update 1D-MRD")
#mrd1d_q = bupd.on_clicked(update_1dmrd)

#=
function update_1dmrd(mrd_q1=mrd_q1, mrd_q3=mrd_q3, ax=ax, ax2=ax2, mrd1d=mrd1d, mrd1d_q=mrd1d_q, mrd_white=mrd_white, x_avg=x_avg, mrd_data=mrd_data, time_middle=time_middle, firstbiggerra=firstbiggerra)
    exc = ax._get_view()
    #@show exc

    xleft = exc[1]
    xright = exc[2]
    ylow = exc[3]
    yhigh = exc[4]

    time_totake = xleft .<= pydates.date2num.(time_middle) .<= xright
    mrd_totake = mrd_data[:, time_totake]

    mrd_avg_plot = zeros(Float64, size(mrd_totake, 1))
    mrd_q1_plot = zeros(Float64, size(mrd_totake, 1))
    mrd_q3_plot = zeros(Float64, size(mrd_totake, 1))

    for i in axes(mrd_avg_plot, 1)
        (mrd_q1_plot[i], mrd_avg_plot[i], mrd_q3_plot[i]) = quantile(filter(!isnan, mrd_totake[i, :]), [0.25, 0.5, 0.75])
    end

    mrd1d.set_xdata(mrd_avg_plot .* 1000)
    c = mrd1d.get_color()
    mrd_white.set_xdata(mrd_avg_plot[firstbiggerra:end] .* 1000)
    mrd1d_q.set_visible(false)
    mrd1d_q = ax2.fill_betweenx(Dates.value.(x_avg) ./ 1000, mrd_q1_plot .* 1000, mrd_q3_plot .* 1000, color=c, alpha=0.3)

    PyPlot.draw()
    return mrd1d_q
end
=#
=#
######################################################
###               PLOT LINEAR MRDS                 ###
######################################################
fig = PyPlot.figure()
ax = fig.add_subplot(111)
for i in 1:30:size(posmrd_data, 2)
    ax.plot(Dates.value.(x[:, 12]) ./ 1e3, posmrd_data[:, i], color="black", alpha=0.5)
end
PyPlot.xscale("log")
