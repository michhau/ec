######################################################
###        MODULE FOR HANDELING TURBULENCE DATA    ###
###             author: Michi Haugeneder           ###
######################################################
#=
Provides functions for the handeling of the measured turbulence data,
i.e. CSAT and IRGASON. Containing I/O and processing.
Also containing rawdata processing for TC data.
=#
module turb

using Dates, DelimitedFiles, CSV, DataFrames, LsqFit,
    Statistics, ProgressMeter, StatsBase, NCDatasets,
    Interpolations
import FastRunningMedian

include( "general.jl")
import .gen
export readperiodfile, nan2missing!, missing2nan!, loadrawgeneric, loadt1raw, loadt2raw, loadkaijoraw,
    createtimestamped3Dwind, csvtodataframe, saveturbasnetcdf, readturbasnetcdf,
    makecontinuous, findnearest, extendtofinertimeseries!, splitdaynight, simplewinddir, qualcontrolflags!,
    sonicqualcontrol, despiking, printmissstats, repositionnanmask!, drdf!, doublerotation,
    detrend, parametersblocksplitting, blockevaluation, interpolatemissing, winddir,
    detectgaps, blockapply, OSHD_SHF, contflux, turbflux, turbfluxdrperperiod, avgflux, advect

######################################################
###                LOADING AND I/O                 ###
######################################################
"""
    readperiodfile(file::String)::DataFrame

Read in periods file to get start/end of measurement periods.
"""
function readperiodfile(file::String)::DataFrame
    @info("Loading measurement period information from ", file)
    intercsv = CSV.File(file; header=0, skipto=2) |> Tables.matrix
    dateformat = DateFormat("yyyy-mm-dd HH:MM:SS")
    df = DataFrame(periodstart=DateTime.(intercsv[:, 1],
            dateformat), periodend=DateTime.(intercsv[:, 2], dateformat),
        comment=intercsv[:, 3])
    println("Done")
    return df
end

"""
    nan2missing!(data::DataFrame)

OLD Replace NaN with missing in DataFrame in place.
Except first column (time)
"""
function nan2missing!(data::DataFrame)
    for col in eachcol(data)
        replace!(col, NaN => missing)
    end
end

"""
    missing2nan!(data::DataFrame)

Replace missing with NaN in DataFrame in-place.
"""
function missing2nan!(data::DataFrame)
    for col in eachcol(data)
        if !(Int <: eltype(col))
            replace!(col, missing => NaN)
        else
            replace!(col, missing => -9999)
        end
    end
end

#########     generic I/O for EC data        ########
"""
    loadrawgeneric(source::String)

Read-in raw turbulence data
"""
function loadrawgeneric(source::String)
    @info("Loading generic eddy covariance raw data")
    labelsfromfile = readlines(source)[2]
    println("Column labels: ", labelsfromfile)
    @info("You need to assign those labels properly in the next step. Time as first column is already extracted.")
    df = CSV.File(source; header=0, skipto=5, tasks=Threads.nthreads()) |> Tables.matrix
    dateformat = DateFormat("yyyy-mm-dd HH:MM:SS.ss")
    timeofmeasure = DateTime.(df[:, 1], dateformat)
    println("Done")
    return timeofmeasure, df[:, 2:end]
    # return df[:,1], df[:,8:end]
end

#########       I/O for tower 1/2 data        ########
"""
    loadt1raw(source::String)

Read-in raw turbulence data from tower 1
"""
function loadt1raw(source::String)
    @info("Loading tower 1 turbulence raw data")
    df = CSV.File(source; header=0, skipto=5, tasks=Threads.nthreads()) |> Tables.matrix
    dateformat = DateFormat("yyyy-mm-dd HH:MM:SS.ss")
    timeofmeasure = DateTime.(df[:, 1], dateformat)
    println("Done")
    return timeofmeasure, df[:, 8:end]
    # return df[:,1], df[:,8:end]
end

"""
    loadt2raw(source::String)

Read-in raw turbulence data from tower 2
"""
function loadt2raw(source::String)
    @info("Loading tower 2 turbulence raw data")
    #df = CSV.File(source; header=0, skipto=5, tasks = Threads.nthreads())|> Tables.matrix
    df = CSV.File(source; header=0, tasks=Threads.nthreads()) |> Tables.matrix
    #dateformat = DateFormat("yyyy-mm-dd HH:MM:SS.ss")
    #timeofmeasure = DateTime.(df[:,1], dateformat)
    println("Done")
    return df[:, 1], df[:, 13:end]
end

#########   I/O for Kaijo turbulence data     ########
"""
    loadkaijoraw(source::String)

Load 3D-ultrasonic rawdata from Kaijo and return timestamps and data.
"""
function loadkaijoraw(source::String)
    @info("Loading Kaijo raw data from ", source)
    @warn("Be careful with axis rotation!! Check image for possible reorientation.")
    df = CSV.File(source; header=0, skipto=5, tasks=Threads.nthreads()) |> Tables.matrix
    dateformat = DateFormat("yyyy-mm-dd HH:MM:SS.ss")
    timeofmeasure = DateTime.(df[:, 1], dateformat)
    println("Done")
    return timeofmeasure, df[:, 3:end]
end

######### I/O for TJK-station turbulence data ########
"""
    load3Dwindraw(source::String)::Array

Not exposed. Load 3D-ultrasonic rawdata into a Float-Array.
"""
function load3Dwindraw(source::String)::Array
    @info("Loading ultrasonic raw data from ", source)
    df = CSV.File(source; header=0, skipto=5, tasks=Threads.nthreads()) |> Tables.matrix
    #df = readdlm(source, ',', Float64, '\n'; skipstart = 4, use_mmap = true)
    #row of first 0 (referrs to corresponding entry in the timestamp file)
    rowfirstzero = findfirst(df[:, 1] .== 0)
    #withdraw rows prior to "rowfirstzero" and return resulting array
    println("Done")
    return df[rowfirstzero:end, :]
end

"""
    loadtimestampsraw(source::String)

Not exposed. Load timestamps for idx=0 from 3D-ultrasonic timestamp-file.
Return [rawtimestamps, rawtimestampdata (for checking)].
"""
function loadtimestampsraw(source::String)
    @info("Loading raw timestamp data for idx=0")
    df = readdlm(source, ',', String, '\n'; skipstart=4, use_mmap=true)
    #extract first column to be DateTime
    dateformat = DateFormat("yyyy-mm-dd HH:MM:SS.ss")
    timeofmeasure = DateTime.(df[:, 1], dateformat)
    #only u-component for quality control
    df = parse.(Float64, df[:, [3, 4]])
    println("Done")
    return timeofmeasure[df[:, 1].==0], df[df[:, 1].==0, 2]
end

"""
    createtimestamps(
    rawwinddata::Array, rawtimestampdata::Vector, rawtimestamps::Vector, timestep)::Vector

Not exposed. Creating the timestamps for the single 3D-wind-records from timestamp-file.
"""
function createtimestamps(
    rawwinddata::Array, rawtimestampdata::Vector, rawtimestamps::Vector, timestep)::Vector
    @info("Creating timestamps...")
    #index for accessing the timestamp array to get the according time stamps
    #if "timestampidx"=2 then skip first row
    timestampidx = 1
    #create array with final timestamps of same size as winddata-rows
    timestamps = Array{DateTime}(undef, size(rawwinddata, 1))

    for i in 1:size(rawwinddata, 1)
        #global timestampidx
        if rawwinddata[i, 1] != 0
            timestamps[i] = timestamps[i-1] + timestep
        else
            #Quality control: u-wind-component
            if rawtimestampdata[timestampidx] == rawwinddata[i, 2]
                timestamps[i] = rawtimestamps[timestampidx]
                timestampidx += 1
            else
                @warn("Quality control at ", rawtimestamps[timestampidx], " (row ", i, ") in windfile failed. Timestamp calculation stopped")
                println("u in timestamp-file: ", rawtimestampdata[timestampidx])
                println("u,v,w,T in windfile: ", rawwinddata[i, 2:end])
                break
                #timestamps[i] = rawtimestamps[timestampidx]
                #timestampidx += 1
            end
        end
        if i == size(rawwinddata, 1)
            println("Timestamps created! All quality checks passed!")
        end
    end
    return timestamps
end

"""
    createtimestamped3Dwind(windrawsource::String, timestampfile::String, timestep)

Unifies all the loading and creating functions for the timestamps.
Timestep corresponds to measurement frequency.
"""
function createtimestamped3Dwind(windrawsource::String, timestampfile::String, timestep)
    println("Creating data with timestamps and 3D-wind.")
    println("HINT: Usually only needs to be executed once. Then read from created file.")
    rawwinddata = load3Dwindraw(windrawsource)
    (rawtimestamps, rawtimestampdata) = loadtimestampsraw(timestampfile)
    timestamps = createtimestamps(
        rawwinddata, rawtimestampdata, rawtimestamps, timestep)
    return timestamps, rawwinddata[:, 2:end]
end

"""
    csvtodataframe(source::String)

Read data from .csv-file and convert it to a DataFrame that is returned.
Usual way of importing data after .csv-file is being created once.
"""
function csvtodataframe(source::String)
    @info("Reading data into DataFrame", source)
    df = DataFrame(CSV.File(source))
    allowmissing!(df)
    for col in eachcol(df)
        replace!(col, NaN => missing)
    end
    return df
end

"""
    saveturbasnetcdf(data::DataFrame, target::String, deflatelvl::Int64=5)

Saving preprocessed turbulence data as NetCDF4
"""
function saveturbasnetcdf(data::DataFrame, target::String, deflatelvl::Int64=5)
    @info("Saving turbulence data to NetCDF4-file")
    ds = NCDataset(target, "c")
    cols = names(data)
    for col in cols
        data_col = collect(data[:, col])
        defDim(ds, col, size(data, 1))
        defVar(ds, col, data_col, (col,); shuffle=true, deflatelevel=deflatelvl)
        #defVar(ds, name, data, ("row", "col", "frame", "u, w"); shuffle=true, deflatelevel=deflatelvl)

    end
    close(ds)
end

"""
    readturbasnetcdf(source::String, perstart::DateTime, perend::DateTime)::DataFrame

Read turbulence data given start and endtime from NetCDF4-file.
"""
function readturbasnetcdf(source::String, perstart::DateTime, perend::DateTime)::DataFrame
    @info("Reading turbulence data from NetCDF-file")
    ds = Dataset(source, "r")

    #find limiting indices
    tmptime = ds["time"]
    time = tmptime[:]
    closeststart = minimum(abs.(skipmissing(time .- perstart)))
    closestend = minimum(abs.(skipmissing(time .- perend)))
    startidx = findfirst(x -> x == closeststart, abs.(time .- perstart))
    endidx = findfirst(x -> x == closestend, abs.(time .- perend))

    #load
    df = DataFrame()
    cols = map(name, varbyattrib(ds))
    for col in cols
        df[:, col] = ds[col][startidx:endidx]
    end
    close(ds)
    return df
end

######################################################
"""
    makecontinuous(df::DataFrame)

Make a continuous time series and fill with missings
"""
function makecontinuous(df::DataFrame)
    @info("Making data series continuous...")
    gaps = detectwrongtimestep(df, Millisecond(50))
    if size(gaps, 1) == 0
        return df
    end
    size_dfnew = Int((df.time[end] - df.time[1]) / Millisecond(50)) + 1
    dfnew = similar(df, size_dfnew)
    dfnew[1:gaps.idx_before_gap[1], :] = df[1:gaps.idx_before_gap[1], :]
    newidx = 0
    timeoffset = Millisecond(0) #offset between original data and new one (new-orig)
    for igap in 1:size(gaps, 1)
        #if igap == 1
        #    newidx = 0 #for debugging
        #end
        beforegap = gaps.idx_before_gap[igap]
        if igap == 1
            newidx += beforegap
        else
            newidx += beforegap - gaps.idx_before_gap[igap-1]
        end
        timebefore = dfnew.time[newidx] #gaps.time_before_gap[igap]
        gaplength_raw = (gaps.gaplength[igap] - timeoffset) / Millisecond(50) - 1
        if mod(gaplength_raw, 1) !== 0.0#!isa(gaplength_raw, Int)
            gaplengthtmp = round(Int, (gaps.gaplength[igap] - timeoffset) / Millisecond(50)) - 1
            timeoffset = (gaplengthtmp - gaplength_raw) * Millisecond(50)
            gaplength = gaplengthtmp
        else
            gaplength = round(Int, gaplength_raw)
            timeoffset = Millisecond(0)
        end
        for jmis in 1:gaplength
            dfnew[newidx+jmis, 1] = timebefore .+ Millisecond(50) .* jmis
            dfnew[newidx+jmis, 2:end] .= missing
        end
        newidx += gaplength
        if igap < size(gaps, 1)
            sizenogap = gaps.idx_before_gap[igap+1] - beforegap
            dfnew[(1:sizenogap).+newidx, 1] .= collect(1:sizenogap) .* Millisecond(50) .+ dfnew[newidx, 1]
            dfnew[(1:sizenogap).+newidx, 2:end] .= df[beforegap+1:gaps.idx_before_gap[igap+1], 2:end]
        else
            dfnew[newidx+1:end, 1] = collect(1:(size(dfnew, 1)-newidx)) .* Millisecond(50) .+ dfnew[newidx, 1]
            dfnew[newidx+1:end, 2:end] = df[beforegap+1:end, 2:end]
        end
    end
    return dfnew
end

"""
    findnearest(targettime::DateTime, searchvec::Vector)::Vector{Int}

Find nearest indices for given targettime in searchvector.
"""
function findnearest(targettime::DateTime, searchvec::Vector)::Vector{Int}
    return findall(x -> x == minimum(abs.(searchvec .- targettime)), abs.(searchvec .- targettime))
end

"""
    extendtofinertimeseries!(finaldata::AbstractVector, finaltime::Vector{DateTime}, initialdata::AbstractVector, initialtime::Vector{DateTime})

Given initialtime-vector and finaltime-vector extend initialdata to finaldata.
No interpolation, just look for nearest neighbor in initialtime.
"""
function extendtofinertimeseries!(finaldata::AbstractVector, finaltime::Vector{DateTime}, initialdata::AbstractVector, initialtime::Vector{DateTime})
    lookidx = findfirst(x->x==minimum(abs.(initialtime .- finaltime[1])), abs.(initialtime .- finaltime[1]))
    for idx in 1:size(finaldata, 1)
        if lookidx == size(initialtime, 1)
            finaldata[idx] = initialdata[lookidx]
        elseif abs(finaltime[idx] - initialtime[lookidx]) <= abs(finaltime[idx] - initialtime[lookidx+1])
            finaldata[idx] = initialdata[lookidx]
        else
            lookidx += 1
            finaldata[idx] = initialdata[lookidx]
        end
    end
end

"""
    splitdaynight(data::DataFrame, daystart::Time, dayend::Time)

Split data into day and night parts
"""
function splitdaynight(data::DataFrame, daystart::Time, dayend::Time)
    dataday = data[daystart.<=Dates.Time.(data.time).<=dayend, :]
    datanight = data[.!(daystart .<= Dates.Time.(data.time) .<= dayend), :]
    return dataday, datanight
end

"""
    simplewinddir(u::Number, v::Number)::Number

Winddirection in deg for u and v. u>0&&v>0 => dir∈]0,90[
"""
function simplewinddir(u::Number, v::Number)::Number
    return mod(rad2deg(atan(-v, u)) + 360, 360)
end

"""
    winddir(datain::DataFrame)::DataFrame

Calculate wind direction of DataFrame
"""
function winddir(datain::DataFrame)::DataFrame
    alpha = zeros(Float64, size(datain, 1))
    for idx in 1:size(datain, 1)
        if !(ismissing(datain.u[idx]) || ismissing(datain.v[idx]))
            abshorwind = sqrt(datain.u[idx]^2 + datain.v[idx]^2)
            beta1 = rad2deg(acos(datain.u[idx] / abshorwind))
            beta2 = rad2deg(asin(datain.v[idx] / abshorwind))
            if beta2 < 0
                alpha[idx] = 360 - beta1
            else
                alpha[idx] = beta1
            end
        else
            alpha[idx] = NaN
        end
    end
    return DataFrame(time=datain.time, α=alpha)
end

"""
    qualcontrolflags!(evaldf::DataFrame)

Replace bad data with missing according to sonic and gas analyser
diagnostic flags and unphysical data
"""
function qualcontrolflags!(evaldf::DataFrame)
    if count(x -> x == "diagsonic", names(evaldf)) > 0
        evaldf[findall(x -> x > 42, collect(skipmissing(evaldf.diagsonic))), 2:end] .= missing
        @info("Sonic diagnostic flag found and used for quality control")
    end
    if count(x -> x == "diagirg", names(evaldf)) > 0
        evaldf[findall(x -> x > 42, collect(skipmissing(evaldf.diagirg))), 2:end] .= missing
        @info("Gas Analyser diagnostic flag found and used for quality control")
    end
end

"""
    sonicqualcontrol(data::DataFrame, umin=-15.00001, umax=15.00001,
    vmin=-15.00001, vmax=15.00001, wmin=-2.500001, wmax=2.500001, Tmin=-7.00001,
    Tmax=21.00001, qh2omin=0, qh2omax=25)

Quality control of the ultrasonic data based on threshold value.
Set to missing
"""
function sonicqualcontrol(data::DataFrame, umin=-15.00001, umax=15.00001,
    vmin=-15.00001, vmax=15.00001, wmin=-2.500001, wmax=2.500001, Tmin=-7.00001,
    Tmax=21.00001, qh2omin=0, qh2omax=25)
    @info("Quality control for u,v,w,T...")
    k = 0
    for i in 1:size(data, 1)
        if isless(umax, data.u[i]) || isless(vmax, data.v[i]) || isless(wmax, data.w[i]) || isless(Tmax, data.T[i]) ||
           isless(data.u[i], umin) || isless(data.v[i], vmin) || isless(data.w[i], wmin) || isless(data.T[i], Tmin)
            data[i, 2:end] .= missing
            k = k + 1
        end
    end
    println(k, " rows set to missing (", round(k / size(data, 1) * 1000, digits=2), "‰)")
    if any(names(data) .== "h2o") #IRGASON
        @info("IRGASON detected")
        nrirgmis = 0
        for i in 1:size(data, 1)
            cond = isless(data.h2o[i], qh2omin) || isless(qh2omax, data.h2o[i])
            if cond
                data[i, [:h2o, :co2]] .= missing
            end
            nrirgmis = count(cond)
        end
        println("IRG: ", nrirgmis, " rows (only :h2o and :co2) set to missing (", round(nrirgmis / size(data, 1) * 1000, digits=2), "‰)")
    end
    return data
end

"""
    despiking(datain::DataFrame, windowwidth=6000, maxsteps=10, breakcrit=1.05)::DataFrame

Despiking algorithm after Sigmund et al. (2022)
"""
function despiking(datain::DataFrame, windowwidth=6000, maxsteps=10, breakcrit=1.05)::DataFrame
    @info("Despiking...")
    #number of elements for the running median
    if windowwidth % 2 == 0 #make odd, so running median returns same length
        windowwidth += 1
    end

    #criterion to compare to
    criterion = 6 / 0.6745

    turb.missing2nan!(datain)
    disallowmissing!(datain)

    #BitVector containing 1 (spike) or 0 (no spike)
    spike = BitVector(undef, size(datain, 1))
    fill!(spike, 0)
    spikeirg = similar(spike)
    fill!(spikeirg, 0)

    #create vector to look only next to detected spikes in next iterations
    tolook = BitVector(undef, size(datain, 1))
    fill!(tolook, 1)

    if any(names(datain) .== "h2o")
        colstodespike = [:u, :v, :w, :T, :h2o]
    else
        colstodespike = [:u, :v, :w, :T]
    end

    devtomedian = datain[:, colstodespike]

    nrsteps = 0
    nrspikes = zeros(Int64, maxsteps)
    nrspikesirg = zeros(Int64, maxsteps)
    nrspikestot = zeros(Int64, maxsteps)
    mad_data = zeros(Float64, size(datain, 1), length(colstodespike))

    while nrsteps < maxsteps
        nrsteps += 1
        for (idx, icol) in enumerate(colstodespike)
            tmp = FastRunningMedian.running_median(datain[:, icol], windowwidth)
            devtomedian[:, icol] = abs.(datain[:, icol] .- tmp)
            mad_data[:, idx] = FastRunningMedian.running_median(devtomedian[:, icol], windowwidth)#mad(filter(!isnan, datain[:,icol]))
        end

        #set previously detected spikes to equal the median
        devtomedian[spike, [:u, :v, :w, :T]] .= 0
        if any(colstodespike .== "h2o")
            devtomedian[spikeirg, :h2o] .= 0
        end

        progtext = string("iteration ", nrsteps, " ")
        p = Progress(size(devtomedian, 1) - 2, desc=progtext)
        Threads.@threads for jdx in 2:size(devtomedian, 1)-1
            if tolook[jdx]
                leftside = 0
                leftsideirg = 0
                for jcol in 1:size(devtomedian, 2)
                    a = devtomedian[jdx, jcol] - mean(filter(!isnan, [devtomedian[jdx-1, jcol], devtomedian[jdx+1, jcol]]))
                    if colstodespike[jcol] != :h2o
                        leftside += a / mad_data[jdx, jcol]
                    else
                        leftsideirg = a / mad_data[jdx, jcol]
                    end
                end
                if leftside > criterion
                    spike[jdx] = 1
                end
                if leftsideirg > criterion
                    spikeirg[jdx] = 1
                end
                next!(p)
            end
        end
        nrspikes[nrsteps] = count(spike)
        nrspikesirg[nrsteps] = count(spikeirg)
        nrspikestot[nrsteps] = nrspikes[nrsteps] + nrspikesirg[nrsteps]
        @show nrspikes[nrsteps]
        @show nrspikesirg[nrsteps]
        @show nrspikestot[nrsteps]

        #only check neigbouring to spikes in next step
        spiketot = spike .| spikeirg
        spikeloc = findall(x -> x, spiketot)
        fill!(tolook, 0)
        tolook[spikeloc.-1] .= 1
        tolook[spikeloc.+1] .= 1

        devtomedian[spike, [:u, :v, :w, :T]] .= 0
        if nrspikesirg[nrsteps] != 0
            devtomedian[spikeirg, :h2o] .= 0
        end
        datain[spike, [:u, :v, :w, :T]] .= NaN
        if nrspikesirg[nrsteps] != 0
            datain[spikeirg, :h2o] .= NaN
        end
        if nrsteps > 1 && nrspikestot[nrsteps] / nrspikestot[nrsteps-1] < breakcrit
            nrsteps = maxsteps + 1 # break
        end
    end
    return datain
end

"""
    printmissstats(df::DataFrame)

Show percentages for u,v,w,T for missing data
"""
function printmissstats(df::DataFrame)
    missingu = count(x -> ismissing(x), df.u)
    missingv = count(x -> ismissing(x), df.v)
    missingw = count(x -> ismissing(x), df.w)
    missingT = count(x -> ismissing(x), df.T)
    println("Missing data (u,v,w,T [‰])")
    println(round(1000 * missingu / size(df, 1), digits=2), ", ",
        round(1000 * missingv / size(df, 1), digits=2), ", ",
        round(1000 * missingw / size(df, 1), digits=2), ", ",
        round(1000 * missingT / size(df, 1), digits=2))
end

"""
    repositionnanmask!(data::DataFrame; timetoreposition=Minute(45))

Apply a NaN-mask to tower1 and tower2 data when sensors were repositioned.
Only tower1 and tower2!
"""
function repositionnanmask!(data::DataFrame; timetoreposition=Minute(45))
    newpositions = [DateTime(2021, 05, 25, 14, 00, 00), DateTime(2021, 05, 28, 09, 00, 00),
        DateTime(2021, 05, 31, 10, 23, 00), DateTime(2021, 06, 08, 09, 45, 00)]
    for i in 1:length(newpositions)
        idcstonan = (newpositions[i]-timetoreposition) .<= data.time .< newpositions[i]
        for j in 1:size(data, 2)
            if Float64 <: eltype(data[:,j])
                data[idcstonan, j] .= NaN
            end
        end
    end
end

"""
    drdf!(data::DataFrame; blockdur=Minute(30), periodwise=true, gapthresh=Minute(10))

Perform in-place double rotation on DataFrame
(using col. names 'u','v','w') to set mean(v)=mean(w)=0., blockdur=duration of blocks, periodwise refers to gaps.
"""
function drdf!(data::DataFrame; blockdur=Minute(30), periodwise=true, gapthresh=Minute(10))
    endidcs = zeros(Int64, 0)
    startidcs = zeros(Int64, 0)
    currnanidx = 1
    println("Double rotation for blocks of ", blockdur)
    blockduridx = round(Int, blockdur/Millisecond(50))
    startidx = 1
    if periodwise
        gaps = detectnanandgap(data, gapthresh)
        nanendidcs = gaps.idx_before_gap
        while startidx < size(data, 1)
            push!(startidcs, startidx)
            if currnanidx<=size(nanendidcs, 1) && nanendidcs[currnanidx] <= startidx+blockduridx-1
                push!(endidcs, nanendidcs[currnanidx])
                startidx = gaps.idx_after_gap[currnanidx]
                currnanidx += 1
            else
                push!(endidcs, startidx+blockduridx-1)
                startidx+=blockduridx-1
            end
        end
        endidcs[end] = size(data, 1)
        @info("Double rotation period-wise. Considering data gaps (e.g. due to reposition).")
        println(size(gaps, 1) + 1, " periods")
    else
        @info("Performing in-place Double Rotation without considering data gaps (DR over gaps as well).")
        while startidx < size(data, 1)
            push!(startidcs, startidx)
            push!(endidcs, startidx+blockduridx-1)
            startidx+=blockduridx-1
        end
        endidcs[end] = size(data, 1)

    end

    for per in axes(endidcs, 1)
        startidx = startidcs[per]
        endidx = endidcs[per]

        datatouse = data[startidx:endidx, :]

        #creating necessary temporary arrays
        data1 = Array{Union{Missing,Float64}}(missing, size(datatouse, 1), 2)
        data2 = Array{Union{Missing,Float64}}(missing, size(datatouse, 1), 3)

        #calculating averages
        meanu = mean(filter(!isnan, skipmissing(datatouse[:, :u])))
        meanv = mean(filter(!isnan, skipmissing(datatouse[:, :v])))
        meanw = mean(filter(!isnan, skipmissing(datatouse[:, :w])))

        #=
        #calculate mean windspeed and -direction
        mean_wndspd = sqrt(meanu^2 + meanv^2 + meanw^2)
        mean_dir = mod(rad2deg(atan(-meanv, meanu)), 360)
        @info("Info:", mean_wndspd, mean_dir)
        =#

        #calculating first rotation angle alpha[rad]
        alpha = atan(meanv, meanu)
        #@show rad2deg(alpha)
        #rotating the windfield to obtain meanv=0 (w-component stays the same)
        data1[:, 1] = datatouse[:, :u] .* cos(alpha) + datatouse[:, :v] .* sin(alpha)
        data1[:, 2] = -datatouse[:, :u] .* sin(alpha) + datatouse[:, :v] .* cos(alpha)
        #calculating second rotation angle beta[rad]
        beta = atan(meanw, mean(filter(!isnan, skipmissing(data1[:, 1]))))
        #@show rad2deg(beta)
        #rotating the windfield to obtain meanw=0 (v-component stays the same)
        data2[:, 1] = data1[:, 1] .* cos(beta) + datatouse[:, :w] .* sin(beta)
        data2[:, 2] = data1[:, 2]
        data2[:, 3] = -data1[:, 1] .* sin(beta) + datatouse[:, :w] .* cos(beta)

        #overwrite the input-data
        data[startidx:endidx, :u] = data2[:, 1]
        data[startidx:endidx, :v] = data2[:, 2]
        data[startidx:endidx, :w] = data2[:, 3]
    end
end

"""
    doublerotation(data)

Perform the double rotation for given 'data' (time,u,v,w) to set mean(v)=mean(w)=0.
Also return mean windspeed and -direction.
"""
function doublerotation(data)
    #creating necessary temporary arrays
    data1 = Array{Union{Missing,Float64}}(missing, size(data, 1), 2)
    data2 = Array{Union{Missing,Float64}}(missing, size(data, 1), 3)

    #calculating averages
    meanu = mean(skipmissing(data[:, 2]))
    meanv = mean(skipmissing(data[:, 3]))
    meanw = mean(skipmissing(data[:, 4]))

    #calculate mean windspeed and -direction
    mean_wndspd = sqrt(meanu^2 + meanv^2 + meanw^2)
    mean_dir = mod(atan(-meanv, meanu), 360)

    #calculating first rotation angle alpha[rad]
    alpha = atan(meanv, meanu)
    #rotating the windfield to obtain meanv=0 (w-component stays the same)
    data1[:, 1] = data[:, 2] .* cos(alpha) + data[:, 3] .* sin(alpha)
    data1[:, 2] = -data[:, 2] .* sin(alpha) + data[:, 3] .* cos(alpha)
    #calculating second rotation angle beta[rad]
    beta = atan(meanw, mean(skipmissing(data1[:, 1])))
    #rotating the windfield to obtain meanw=0 (v-component stays the same)
    data2[:, 1] = data1[:, 1] .* cos(beta) + data[:, 4] .* sin(beta)
    data2[:, 2] = data1[:, 2]
    data2[:, 3] = -data1[:, 1] .* sin(beta) + data[:, 4] .* cos(beta)
    #@info("DR check", mean(skipmissing(data2[:,2])))
    return data2, mean_wndspd, mean_dir
end

"""
    interpolatemissing(data::DataFrame, threshgapsize::Int=20)

Linearly interpolate missing values. If gapsize>threshgapsize => do nothing.
"""
function interpolatemissing(data::DataFrame, threshgapsize::Int=20)
    rowlength = size(data, 1)
    collength = size(data, 2)

    nrreplaced = 0

    #determine rows containing missing values
    missingarr = BitArray(undef, rowlength, collength - 1)
    missingvec = BitArray(undef, rowlength)
    for i in 1:collength-1
        missingarr[:, i] = ismissing.(data[:, i+1]) .| isnan.(data[:, i+1])
    end
    for i in 1:rowlength
        missingvec[i] = any(missingarr[i, :])
    end

    #determine if cols are Int values
    colisint = Int .<: eltype.(data[:, i] for i in 1:collength)

    #handle missing values in first or last entry
    if missingvec[1]
        #check for gap size if < threshgapsize
        nextnonmissing = findnext(x -> x == 0, missingvec, 1)
        if !isnothing(nextnonmissing)
            gapsize = nextnonmissing - 1
            if gapsize <= threshgapsize
                nrreplaced += 1
                for jcol in 2:collength
                    if colisint[jcol]
                        data[1, jcol] = round(Int, mean(skipmissing(data[:, jcol])))
                    else
                        data[1, jcol] = mean(skipmissing(data[:, jcol]))
                    end
                end
            end
        end
    end
    if missingvec[end]
        #check for gap size if < threshgapsize
        prevnonmissing = findprev(x -> x == 0, missingvec, size(missingvec, 1))
        if !isnothing(prevnonmissing)
            gapsize = abs(prevnonmissing - rowlength) - 1
            if gapsize <= threshgapsize
                nrreplaced += 1
                for jcol in 2:collength
                    if colisint[jcol]
                        data[end, jcol] = round(Int, mean(skipmissing(data[:, jcol])))
                    else
                        data[end, jcol] = mean(skipmissing(data[:, jcol]))
                    end
                end
            end
        end
    end
    #handle missing values in the bulk
    j = 1
    while j <= rowlength - 1
        j += 1
        if missingvec[j]
            #check for gap size if < threshgapsize
            nextnonmissing = findnext(x -> x == 0, missingvec, j)
            if !isnothing(nextnonmissing)
                gapsize = nextnonmissing - j
                if gapsize <= threshgapsize
                    nrreplaced += gapsize
                    for icol in 2:collength
                        fct = (data[nextnonmissing, icol] - data[j-1, icol]) / (gapsize + 1)
                        offs = data[j-1, icol]
                        if colisint[icol]
                            for kstep in 1:gapsize
                                data[j+kstep-1, icol] = round(Int, offs + fct * kstep)
                            end
                        else
                            for kstep in 1:gapsize
                                data[j+kstep-1, icol] = offs + fct * kstep
                            end
                        end
                    end
                end
                j = j + gapsize - 1
            else
                j = rowlength #all missing to end, so skip rest    
            end
        end
    end
    replaced = nrreplaced / count(missingvec)
    println("Replaced ", round(replaced * 100, digits=2),
        "% of missing values.")
    return data
end

"""
    detrend(invec::Vector)::Vector

Following Stull p.307ff. Detrend data with linear fit
"""
function detrend(invec::Vector)::Vector
    model(x, p) = p[1] * x .+ p[2]
    p0 = [0.0001, 3.0]
    fit = curve_fit(model, collect(1:size(invec, 1)), invec, p0)
    param = fit.param
    return invec - (param[1] * collect(1:size(invec, 1)) .+ param[2])
end

"""
    parametersblocksplitting(blocklength, timestep, sizedf::Integer)

Determine parameters for block splitting
"""
function parametersblocksplitting(blocklength, timestep, sizedf::Integer)
    rowsperblock = ceil(Int, Millisecond(blocklength) / timestep)
    M = floor(Int, log2(rowsperblock))
    Mdur = canonicalize(Dates.CompoundPeriod(2^M * timestep))
    mrd_discarded = canonicalize(blocklength - Mdur)
    nrblcks = ceil(Int, sizedf / rowsperblock)
    lengthlastblck = canonicalize(Dates.CompoundPeriod((sizedf -
                                                        rowsperblock * (nrblcks - 1)) * timestep))
    return rowsperblock, M, Mdur, mrd_discarded, nrblcks, lengthlastblck
end

"""
    blockevaluation(evaldf::DataFrame, blocklength, timestep)

Evaluate data blocks and get information on blockdur, wT, mean wind, meanT, winddir.
"""
function blockevaluation(evaldf::DataFrame, blocklength, timestep)
    #parameters for splitting the data into blocks
    (rowsperblock, M, Mdur, mrd_discarded, nrblocks, lengthlastblock) =
        turb.parametersblocksplitting(blocklength, timestep, size(evaldf, 1))

    println("Taking ", nrblocks - 1, " blocks of ", blocklength, " + ")
    println("last block length is ", lengthlastblock, ". (currently skipped)")

    #create DataFrame to be filled with data from the blocks
    eval_block = DataFrame(time_middle=DateTime[], time_start=DateTime[],
        duration=Dates.CompoundPeriod[], nrrecs=Int[], rate_missing=Float64[], mean_T=Float64[],
        mean_wndspd=Float64[], mean_dir=Float64[], wT=Float64[])

    rotdata = DataFrame(time=DateTime[], u=Float64[], v=Float64[], w=Float64[], T=Float64[])
    allowmissing!(rotdata)
    @info("Calculation of properties for blocks...")
    #iterator over the blocks
    @showprogress for i in 1:nrblocks-1
        #println("Block ", i, "/", nrblocks)
        startidx = 1 + (i - 1) * rowsperblock
        if i < nrblocks
            endidx = startidx + rowsperblock - 1
        else
            endidx = size(evaldf, 1)
        end

        #nr. of (non-)missing records
        nrmis = count(ismissing, evaldf.u[startidx:endidx])
        nrnonmis = length(evaldf.u[startidx:endidx]) - nrmis
        #@info("Percent missing", round(Int, 100*nrmis/(nrmis+nrnonmis)))

        #mean virtual temperature
        meanT = mean(skipmissing(evaldf.T[startidx:endidx]))

        #creating temporary array for the following calculations
        evaldf_dr = Array{Union{Missing,Float64}}(missing, endidx - startidx + 1, 4)
        #evaldf_dr = zeros(Float64, endidx-startidx + 1, 4)
        evaldf_dr[:, 4] = evaldf.T[startidx:endidx]

        #double rotation
        (evaldf_dr[:, 1:3], mean_wndspd, mean_dir) =
            turb.doublerotation(evaldf[startidx:endidx, 1:4])  #components w/o T-col

        #interpolate the missing values
        if nrmis / length(evaldf_dr[:, 2]) < 0.40
            evaldf_dr = turb.interpolatemissing(evaldf_dr)
        else
            evaldf_dr[:, 2:end] .= 0
        end

        #detrend
        for j in 2:size(evaldf_dr, 2)
            evaldf_dr[:, j] = turb.detrend(collect(evaldf_dr[:, j]))
        end

        #calculate flux for the block
        wT = mean(skipmissing(evaldf_dr[:, 3] .* evaldf_dr[:, 4]))

        #create timestamps
        time_start = evaldf.time[startidx]
        time_middle = evaldf[round(Int, (startidx + endidx) / 2), 1]
        block_dur = canonicalize(Dates.CompoundPeriod((endidx - startidx) * timestep))

        #write all the variables to the DataFrame
        push!(eval_block, [time_middle, time_start, block_dur, nrnonmis + nrmis,
            nrmis / (nrnonmis + nrmis), meanT, mean_wndspd, mean_dir, wT])

        #write the rotated data to the dataframe
        for idx in 1:rowsperblock
            push!(rotdata, [evaldf[startidx+idx-1, 1], evaldf_dr[idx, 1],
                evaldf_dr[idx, 2], evaldf_dr[idx, 3], evaldf[startidx+idx-1, 5]])
        end
    end
    return eval_block, rotdata
end

"""
    detectgaps(data::DataFrame, gapthresh::Period)

detect gaps bigger than 'gapthresh' in data
"""
function detectgaps(data::DataFrame, gapthresh::Period)
    gaps = DataFrame(idx_before_gap=Int64[], time_before_gap=DateTime[], gaplength=Period[])
    for i in 1:size(data, 1)-1
        if data.time[i+1] - data.time[i] > gapthresh
            push!(gaps, [i, data.time[i], data.time[i+1] - data.time[i]])
        end
    end
    return gaps
end

"""
    detectnanandgap(data::DataFrame, gapthresh::Period)

Detect also nan-filled gaps addiitonally to above function 'detectgaps'
"""
function detectnanandgap(data::DataFrame, gapthresh::Period)
    gaps = DataFrame(idx_before_gap=Int64[], idx_after_gap=Int64[])
    for i in 1:size(data, 1)-1
        if data.time[i+1] - data.time[i] > gapthresh
            push!(gaps, [i, i+1])
        end
    end

    idxthresh = gapthresh/Millisecond(50)
    nanidcs = findall(x->(ismissing(x) | isnan(x)), data.u)

    lengthcurrser = 1
    serstartidx = 1

    for i in 2:length(nanidcs)
        if nanidcs[i] == nanidcs[i-1] +1
            lengthcurrser += 1
        else
            if lengthcurrser >= idxthresh
                push!(gaps, [serstartidx-1, nanidcs[i]])
            end
            lengthcurrser = 1
            serstartidx = nanidcs[i]
        end
    end

    sort(gaps, :idx_before_gap)
end

"""
    detectwrongtimestep(data::DataFrame, timestep::Period)

detect timesteps not equal to 'timestep' in data
"""
function detectwrongtimestep(data::DataFrame, timestep::Period)
    gaps = DataFrame(idx_before_gap=Int64[], time_before_gap=DateTime[], gaplength=Period[])
    for i in 1:size(data, 1)-1
        if data.time[i+1] - data.time[i] != timestep
            push!(gaps, [i, data.time[i], data.time[i+1] - data.time[i]])
        end
    end
    return gaps
end

"""
    blockapply(f::Function, datasource::DataFrame, startidx::Integer, endidx::Integer, dur::Period)

Apply a function func (i.e. mrd) to blocks (duration -> dur) of data from
datasource. Last block shorter. If endidx==-1 use all data.
"""
function blockapply(f::Function, datasource::DataFrame, startidx::Integer, endidx::Integer, dur::Period)
    ###################################################################
    #splitting the data into blocks
    if endidx == -1
        endidx = size(datasource, 1)
    end
    data = datasource[startidx:endidx, :]
    measuretimestep = data[3, 1] - data[2, 1]
    rowsperblock = ceil(Int, Millisecond(dur) / measuretimestep)
    nrblocks = ceil(Int, size(data, 1) / rowsperblock)
    lengthlastblock = canonicalize(Dates.CompoundPeriod((size(data, 1) - rowsperblock * (nrblocks - 1)) * measuretimestep))
    println("Applying Block stuff to ", nrblocks - 1, " blocks of ", dur, " + ")
    println("last block length is ", lengthlastblock, ".")
    ##################################################################
    #apply the function func on the data on the block
    for i in 1:nrblocks
        startidx = 1 + (i - 1) * rowsperblock
        if i < nrblocks
            endidx = startidx + rowsperblock
        else
            endidx = size(data, 1)
        end
        data[startidx:endidx] = f(data[startidx:endidx])
    end
    ##################################################################
    #return the transformed data
    return data
end

######################################################
###                FLUX CALCULATIONS               ###
######################################################
"""
    OSHD_SHF(ta::Vector, tsrf::Vector, ua::Vector, z0::Float64,
    zt1::Number, zu1::Number, dh::Number=0.0)

Calculate 'modelled' flux according to OSHD-FSM with measured data
"""
function OSHD_SHF(ta::Vector, tsrf::Vector, ua::Vector, z0::Float64,
    zt1::Number, zu1::Number, dh::Number=0.0)
    #parameters
    ρ = 1.2 #kg m^{-3}
    cp = 1005 #J kg^{-1} K^{-1}
    bstb = 5 #atmospheric stability
    g = 9.81 #m s^{-2}

    #step1: friction velocity
    cd = 0.4 / (log((zu1 - dh) / z0)^2)
    ustar = sqrt(cd) .* ua

    #step2: bulk Richardson number
    RiB = 9.81 .* (ta .- tsrf) .* (zu1 - dh)^2 ./ ((zt1 - dh) .* ta .* ua.^2)
    RiB[RiB .> 0.2] .= 0.2 #limit stability

    fh = zeros(Float64, length(RiB))
    #step3: stability factor
    for idx in 1:length(RiB)
        if RiB[idx] > 0
            fh[idx] = (1 + 3 * bstb * RiB[idx] * sqrt(1 + bstb * RiB[idx]))^(-1)
        else #RiB<0
            fh[idx] = 1 - 3 * bstb * RiB[idx] * (1 + 3 * bstb^2 * cd * sqrt(-RiB[idx] * zu1 / z0))^(-1)
        end
    end

    #step4: eddy diffusivity
    kh = fh .* 0.4 .* ustar ./ (log(zt1 / (0.1 * z0)))

    #step5: sensible heat flux H
    h = cp * ρ .* kh .* (tsrf .- ta)

    return h, RiB, tsrf.-ta
end

"""
    contflux(X::Vector, Y::Vector, numofele::Integer)::Vector

Calculate continuous flux (X'Y')
"""
function contflux(X::Vector, Y::Vector, numofele::Integer)::Vector
    Xavg = gen.movingaverage(X, numofele)
    Yavg = gen.movingaverage(Y, numofele)
    Xdash = X - Xavg
    Ydash = Y - Yavg
    return Xdash .* Ydash
end

"""
    turbflux(data::DataFrame, reyavgtime::Period, p_over_p0::Float64=1013 / 798)::DataFrame

Calculate turbulent flux and Obukhov length
"""
function turbflux(data::DataFrame, reyavgtime::Period, p_over_p0::Float64=1013 / 798)::DataFrame
    leng = size(data, 1)
    fluxout = DataFrame("time" => data.time, "wT" => fill(NaN, leng), "wq" => fill(NaN, leng),
        "uw" => fill(NaN, leng), "vw" => fill(NaN, leng), "uT" => fill(NaN, leng),
        "u_star" => fill(NaN, leng), "uu" => fill(NaN, leng), "vv" => fill(NaN, leng),
        "ww" => fill(NaN, leng), "TT" => fill(NaN, leng), "tke" => fill(NaN, leng), "T_pot_20Hz" => fill(NaN, leng),
        "L_highfreq" => fill(NaN, leng))

    timestep = data.time[2] - data.time[1]
    #check
    checkidx = round(Int, (size(data, 1) - 1) * rand()) + 1
    timestepcheck = data.time[checkidx] - data.time[checkidx-1]
    if timestep != timestepcheck
        @warn("timestep check (sonic data) was not successfull!")
    end
    avgidcs = maximum([1 round(Int, Millisecond(reyavgtime) / Millisecond(timestep))])
    avg_u = gen.movingaverage(data.u, avgidcs)
    avg_v = gen.movingaverage(data.v, avgidcs)
    avg_w = gen.movingaverage(data.w, avgidcs)
    avg_T = gen.movingaverage(data.T, avgidcs)
    avg_T_pot = avg_T .* (p_over_p0^(2 / 7))
    dev_u = data.u .- avg_u
    dev_v = data.v .- avg_v
    dev_w = data.w .- avg_w
    dev_T = data.T .- avg_T
    fluxout.wT = dev_w .* dev_T
    fluxout.uw = dev_u .* dev_w
    fluxout.vw = dev_v .* dev_w
    fluxout.uT = dev_u .* dev_T
    fluxout.u_star = (fluxout.uw .^ 2 + fluxout.vw .^ 2) .^ (1 / 4)
    fluxout.uu = dev_u .^ 2
    fluxout.vv = dev_v .^ 2
    fluxout.ww = dev_w .^ 2
    fluxout.TT = dev_T .* dev_T
    fluxout.tke = (fluxout.uu .+ fluxout.vv .+ fluxout.ww) .* 0.5
    fluxout.T_pot_20Hz = data.T .* (p_over_p0^(2 / 7))
    fluxout.L_highfreq = -((fluxout.u_star .^ 3) .* avg_T_pot) ./ fluxout.wT ./ (0.4 * 9.81)#Obukhov-length
    if "h2o" in names(data)
        avg_q = gen.movingaverage(data.h2o, avgidcs)
        dev_q = data.h2o .- avg_q
        fluxout.wq = dev_w .* dev_q
    end
    return fluxout
end

"""
    turbfluxdrperperiod(data::DataFrame, reyavgtime::Period, dr_period::Period, p_over_p0::Float64=1013 / 798)::DataFrame

Calculate turbulent flux and double rotate every given interval
"""
function turbfluxdrperperiod(data::DataFrame, reyavgtime::Period, dr_period::Period, p_over_p0::Float64=1013 / 798)::DataFrame
    leng = size(data, 1)
    fluxout = DataFrame("time" => data.time, "wT" => fill(NaN, leng), "wq" => fill(NaN, leng),
        "uw" => fill(NaN, leng), "vw" => fill(NaN, leng), "uT" => fill(NaN, leng),
        "u_star" => fill(NaN, leng), "uu" => fill(NaN, leng), "vv" => fill(NaN, leng),
        "ww" => fill(NaN, leng), "tke" => fill(NaN, leng), "T_pot_20Hz" => fill(NaN, leng),
        "L_highfreq" => fill(NaN, leng))

    timestep = data.time[2] - data.time[1]
    #check
    checkidx = round(Int, (size(data, 1) - 1) * rand()) + 1
    timestepcheck = data.time[checkidx] - data.time[checkidx-1]
    if timestep != timestepcheck
        @warn("timestep check (sonic data) was not successfull!")
    end
    avgidcs = maximum([1 round(Int, Millisecond(reyavgtime) / Millisecond(timestep))])
    recsperdrperiod = round(Int, dr_period / timestep)
    startidx = 1
    endidx = 1
    nrperiods = ceil(Int, size(data, 1) / recsperdrperiod)
    for i in 1:nrperiods
        startidx = 1 + (i - 1)*recsperdrperiod
        if i !== nrperiods
            endidx = startidx + recsperdrperiod -1
        else
            endidx = size(data, 1)
        end

        datatmp = copy(data[startidx:endidx, :])
        drdf!(datatmp);
        avg_u = gen.movingaverage(datatmp.u, avgidcs)
        avg_v = gen.movingaverage(datatmp.v, avgidcs)
        avg_w = gen.movingaverage(datatmp.w, avgidcs)
        avg_T = gen.movingaverage(datatmp.T, avgidcs)
        avg_T_pot = avg_T .* (p_over_p0^(2 / 7))
        dev_u = datatmp.u .- avg_u
        dev_v = datatmp.v .- avg_v
        dev_w = datatmp.w .- avg_w
        dev_T = datatmp.T .- avg_T
        fluxout.wT[startidx:endidx] = dev_w .* dev_T
        fluxout.uw[startidx:endidx] = dev_u .* dev_w
        fluxout.vw[startidx:endidx] = dev_v .* dev_w
        fluxout.uT[startidx:endidx] = dev_u .* dev_T
        fluxout.u_star[startidx:endidx] = (fluxout.uw[startidx:endidx] .^ 2 + fluxout.vw[startidx:endidx] .^ 2) .^ (1 / 4)
        fluxout.uu[startidx:endidx] = dev_u .^ 2
        fluxout.vv[startidx:endidx] = dev_v .^ 2
        fluxout.ww[startidx:endidx] = dev_w .^ 2
        fluxout.tke[startidx:endidx] = (fluxout.uu[startidx:endidx] .+ fluxout.vv[startidx:endidx] .+ fluxout.ww[startidx:endidx]) .* 0.5
        fluxout.T_pot_20Hz[startidx:endidx] = datatmp.T .* (p_over_p0^(2 / 7))
        fluxout.L_highfreq[startidx:endidx] = -((fluxout.u_star[startidx:endidx] .^ 3) .* avg_T_pot) ./ fluxout.wT[startidx:endidx] ./ (0.4 * 9.81)#Obukhov-length
        if "h2o" in names(data)
            avg_q = gen.movingaverage(datatmp.h2o, avgidcs)
            dev_q = datatmp.h2o .- avg_q
            fluxout.wq[startidx:endidx] = dev_w .* dev_q
        end
    end
    return fluxout
end

"""
    avgflux(data::DataFrame, peri::Period)

Average the turbulent fluxes in the DataFrame
"""
function avgflux(data::DataFrame, peri::Period)
    ele = round(Int, Millisecond(peri) / Millisecond(50))
    fluxavg = similar(data)
    fluxavg.time = data.time
    for iname in names(data)[2:end]
        fluxavg[:, iname] = gen.movingaverage(data[:, iname], ele)
    end
    return fluxavg
end

"""
    obukhov(fluxin::DataFrame, avgtime::Period)::DataFrame

Calculate Obukhov length for moving-averaged inputs
"""
function obukhov(fluxin::DataFrame, avgtime::Period)::DataFrame
    leng = size(fluxin, 1)
    Lout = DataFrame("time" => fluxin.time, "L" => fill(NaN, leng))
    fluxavg = avgflux(fluxin, avgtime)
    Lout.L = -((fluxavg.u_star .^ 3) .* fluxavg.T_pot_20Hz) ./ fluxavg.wT ./ (0.4 * 9.81)#Obukhov-length
    return Lout
end

"""
    advect(ec1::DataFrame, ec2::DataFrame, quantity::String, windcomp::String, avgtime::Period, limangle::Number, Δy::Number)::DataFrame

TBW
"""
function advect(ec1::DataFrame, ec2::DataFrame, quantity::String, windcomp::String, avgtime::Period, limangle::Number, Δy::Number)::DataFrame
    #from Mott20(HEFEX): HA = -\frac{ΔT}{Δy}̅v with ̅v=(v₁+v₂)/2; ΔT = T₁-T₂
    @info("Advection. Direction: (1.quantity-2.quantity)*winddir")

    #generate output DataFrame
    advdf = DataFrame(time=ec1.time, adv=fill(NaN, size(ec1.time)), meanwind=fill(NaN, size(ec1.time)), Δq=fill(NaN, size(ec1.time)))

    if !(windcomp == "u" || windcomp == "v")
        @error("windcomponent must be u or v!")
        return advdf
    end

    #check arguments
    if !in(quantity, names(ec1)) || !in(quantity, names(ec2))
        @error("Advected 'quantity' needs to be a column of 'ec1' and 'vec2'!")
        return advdf
    end
    if !in(windcomp, names(ec1)) || !in(windcomp, names(ec2))
        @error("Wind direction 'windcomp' needs to be a column of 'vec1' and 'vec2'!")
        return advdf
    end
    #check timeshift between ec1.time and ec2.time
    timeshiftidx1 = ec1.time[1] - ec2.time[1]
    timeshiftidxend = ec1.time[end] - ec2.time[end]
    @show timeshiftidx1, timeshiftidxend

    avgidx = round(Int, Millisecond(avgtime) / Millisecond(50))

    diffquant = gen.movingaverage(ec1[:, quantity], avgidx) - gen.movingaverage(ec2[:, quantity], avgidx)
    meanwind = (gen.movingaverage(ec1[:, windcomp], avgidx) + gen.movingaverage(ec2[:, windcomp], avgidx)) ./ 2
    if windcomp == "u"
        meancrosswind = (gen.movingaverage(ec1.v, avgidx) + gen.movingaverage(ec2.v, avgidx)) ./ 2
    elseif windcomp == "v"
        meancrosswind = (gen.movingaverage(ec1.u, avgidx) + gen.movingaverage(ec2.u, avgidx)) ./ 2
    end

    #calculate limit for angle corridor
    @info("Defining angle corridor on ec2")
    corrlim = sqrt.(meanwind .^ 2 + meancrosswind .^ 2) .* cos(deg2rad(limangle))
    meanwind[-corrlim.<meanwind.<corrlim] .= NaN

    #actual advection calculation
    advdf.adv = diffquant .* meanwind ./ Δy
    advdf.meanwind = meanwind
    advdf.Δq = diffquant


    return advdf
end
end #module
