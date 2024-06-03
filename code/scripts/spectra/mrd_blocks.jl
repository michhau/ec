######################################################
### MULTI_RESOLUTION FLUX DECOMPOSITION FOR BLOCKS ###
###            author: Michi Haugeneder            ###
######################################################
#=
Calculate 'normal' and quasi-continuos MRD. Load data with
load_data.jl
Calculate statistics for given blocks and save figures
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter, Suppressor
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")

importdir = joinpath(@__DIR__, "..", "..")
datapath = "/home/haugened/Documents/data/"
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .gen
PyPlot.pygui(true)

timestep = Millisecond(50)

#parameter for MRD (length = 2^M*timestep)
M = 14
println("Length MRD blocks: ", 2^M * timestep)
MRDshift = round(Int, 0.1 * 2^M) #shift between two blocks of MRD
MRDoutfolder = "/home/haugened/Pictures/MRDs/"

#parameters for statistics using the TJK station
upvalleywind_leftlim = 305 #degree
upvalleywind_rightlim = 10 #degree
downvalleywind_leftlim = 135#degree
downvalleywind_rightlim = 200#degree

#select data and measurement period to be evaluated
evaldf1 = t1irg
evaldf2 = t2irgdb
evaldf3 = t2lcsatdb
evaldf4 = t2ucsatdb
evaldf5 = kaijodf
evaldf6 = tjkmeteodata

#blocklength; first included, last excluded
bllength = Hour(2)

blockstarttime = starttimeperiod - bllength
blockendtime = starttimeperiod
iblock = 0

#iterate over blocks
while blockstarttime <= endtimeperiod
    iblock += 1
    blockstarttime += bllength
    blockendtime += bllength

    #cut data
    blockdata1 = evaldf1[blockstarttime.<=evaldf1.time.<blockendtime, :]
    blockdata2 = evaldf2[blockstarttime.<=evaldf2.time.<blockendtime, :]
    blockdata3 = evaldf3[blockstarttime.<=evaldf3.time.<blockendtime, :]
    blockdata4 = evaldf4[blockstarttime.<=evaldf4.time.<blockendtime, :]
    # blockdata5 = evaldf5[blockstarttime .<= evaldf5.time .< blockendtime,:]
    blockdata6 = evaldf6[blockstarttime.<=evaldf6.time.<blockendtime, :]

    if size(blockdata1, 1) > 2^M && size(blockdata2, 1) > 2^M && size(blockdata3, 1) > 2^M && size(blockdata4, 1) > 2^M
        @suppress begin
            #double rotation
            turb.drdf!(blockdata1)
            turb.drdf!(blockdata2)
            turb.drdf!(blockdata3)
            turb.drdf!(blockdata4)
            #turb.drdf!(blockdata5)

            #calculate wind direction statistics using the TJK station
            upvalleywind_counts = count(x -> upvalleywind_leftlim .<= x || x .<= upvalleywind_rightlim, blockdata6.wind_mean_vector_direction)
            downvalleywind_counts = count(x -> downvalleywind_leftlim .<= x .<= downvalleywind_rightlim, blockdata6.wind_mean_vector_direction)
            total_counts = size(blockdata6, 1)
            upwindper = round(upvalleywind_counts / total_counts * 100, digits=1)
            downwindper = round(downvalleywind_counts / total_counts * 100, digits=1)
            statsstring = string("up-valley: ", upwindper, "%\n", "down-valley: ", downwindper, "%")

            #MRD
            #Mulitresolution Flux Decomposition
            blocklength = ceil(2^M * timestep, Dates.Second)
            (mrd_x_1, mrd_data_1, time_middle_1) = turb.completemrd(blockdata1, "w", "T", M, MRDshift)
            (mrd_time_1, mrd_D_median_1, mrd_D_min_1, mrd_D_max_1, mrd_D_quant1_1, mrd_D_quant3_1) =
                turb.mrdpp(mrd_x_1, mrd_data_1)
            (mrd_x_2, mrd_data_2, time_middle_2) = turb.completemrd(blockdata2, "w", "T", M, MRDshift)
            (mrd_time_2, mrd_D_median_2, mrd_D_min_2, mrd_D_max_2, mrd_D_quant1_2, mrd_D_quant3_2) =
                turb.mrdpp(mrd_x_2, mrd_data_2)
            (mrd_x_3, mrd_data_3, time_middle_3) = turb.completemrd(blockdata3, "w", "T", M, MRDshift)
            (mrd_time_3, mrd_D_median_3, mrd_D_min_3, mrd_D_max_3, mrd_D_quant1_3, mrd_D_quant3_3) =
                turb.mrdpp(mrd_x_3, mrd_data_3)
            (mrd_x_4, mrd_data_4, time_middle_4) = turb.completemrd(blockdata4, "w", "T", M, MRDshift)
            (mrd_time_4, mrd_D_median_4, mrd_D_min_4, mrd_D_max_4, mrd_D_quant1_4, mrd_D_quant3_4) =
                turb.mrdpp(mrd_x_4, mrd_data_4)
            #=if size(blockdata5, 1)>0
                (mrd_x_5, mrd_data_5, time_middle_5) = turb.completemrd(blockdata5, "w", "T", M, MRDshift)
                try
                    (mrd_time_5,mrd_D_median_5,mrd_D_min_5,mrd_D_max_5,mrd_D_quant1_5, mrd_D_quant3_5) = 
                    turb.mrdpp(mrd_x_5, mrd_data_5)
                catch e
                end
            end=#

            #Plot result of MRD
            startstring = Dates.format(blockstarttime, "yyyymmddHHMMSS")
            endstring = Dates.format(blockendtime, "yyyymmddHHMMSS")
            filename = string(startstring, "_to_", endstring, "_MRD.png")

            fig = PyPlot.figure()#figsize=(12,7))
            ax = fig.add_subplot(111)
            ax.set_title(string("MRD ", blockstarttime, " - ", blockendtime))
            ax.set_xlabel(L"t_{avg}~\mathrm{[s]}")
            ax.set_ylabel(L"$C_{w\theta}~[\cdot 10^{-3}~\mathrm{K~m~s^{-1}}]$")
            ax.grid(true)
            PyPlot.xscale("log")
            ax.plot(Dates.value.(mrd_time_1) ./ 1000, mrd_D_median_1 .* 1000, label="T1IRG")
            ax.fill_between(Dates.value.(mrd_time_1) ./ 1000, mrd_D_quant1_1 .* 1000, mrd_D_quant3_1 .* 1000, alpha=0.3)
            ax.plot(Dates.value.(mrd_time_2) ./ 1000, mrd_D_median_2 .* 1000, label="T2IRG")
            ax.fill_between(Dates.value.(mrd_time_2) ./ 1000, mrd_D_quant1_2 .* 1000, mrd_D_quant3_2 .* 1000, alpha=0.3)
            ax.plot(Dates.value.(mrd_time_3) ./ 1000, mrd_D_median_3 .* 1000, label="T2LCSAT")
            ax.fill_between(Dates.value.(mrd_time_3) ./ 1000, mrd_D_quant1_3 .* 1000, mrd_D_quant3_3 .* 1000, alpha=0.3)
            ax.plot(Dates.value.(mrd_time_4) ./ 1000, mrd_D_median_4 .* 1000, label="T2UCSAT")
            ax.fill_between(Dates.value.(mrd_time_4) ./ 1000, mrd_D_quant1_4 .* 1000, mrd_D_quant3_4 .* 1000, alpha=0.3)
            #=if @isdefined mrd_D_median_5
                ax.plot(Dates.value.(mrd_time_5)./1000, mrd_D_median_5.*1000, label="Kaijo")
                ax.fill_between(Dates.value.(mrd_time_5)./1000, mrd_D_quant1_5.*1000, mrd_D_quant3_5.*1000, alpha=0.3)
            end=#
            ax.legend()
            PyPlot.figtext(0.002, 0.005, statsstring)
            #PyPlot.show()
            PyPlot.savefig(joinpath(MRDoutfolder, filename))
            PyPlot.close()
        end #if
    end #suppress
    println(blockstarttime, " - ", blockendtime, " done")
end #while