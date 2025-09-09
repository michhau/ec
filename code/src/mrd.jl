######################################################
###       MODULE CONTAINING FUNCTIONS FOR MRD      ###
###             author: Michi Haugeneder           ###
######################################################
#=
Functions that are needed to compute the
Multi-Resolution Flux Decomposition.
=#

module MRD

using Dates, Statistics, NCDatasets, Distributed
using DataFrames, ProgressMeter

include("turb_data.jl")
include("general.jl")
import .turb
import .gen
export mrd, completemrd, mrdpp, completenomrd, mrd_std, mrd_std_block,
    generateseries, fastnomrd, nomrd, contnomrd,
    catdatanomrd, nomrd2d, write2dmrddatatonetcdf,
    read2dmrdfromnetcdf, create2dmrdtitle

"""
    detectgapsdatatimearray(data_time::Vector, gapthresh::Period)::DataFrame

detect gaps bigger than 'gapthresh' in data (Array)
"""
function detectgapsdatatimearray(data_time::Vector, gapthresh::Period)::DataFrame
    gaps = DataFrame(idx_before_gap=Int64[], time_before_gap=DateTime[], gaplength=Period[])
    for i in axes(data_time, 1)[1:end-1]
        if abs(data_time[i+1] - data_time[i]) > gapthresh
            push!(gaps, [i, data_time[i], data_time[i+1] - data_time[i]])
        end
    end
    return gaps
end

"""
    detectgaps(data::DataFrame, gapthresh::Period)

detect gaps bigger than 'gapthresh' in data
"""
function detectgaps(data::DataFrame, gapthresh::Period)
    gaps = DataFrame(idx_before_gap=Int64[], time_before_gap=DateTime[], gaplength=Period[])
    for i in axes(data, 1)[1:end-1]
        if data.time[i+1] - data.time[i] > gapthresh
            push!(gaps, [i, data.time[i], data.time[i+1] - data.time[i]])
        end
    end
    return gaps
end

"""
    mrd(data_a::Vector, data_b::Vector, M::Integer, Mx::Integer)

NEW Orthogonal Multiresolution Flux Decomposition. Adapted from Ivana Stiperski's code.
See Vickers&Mahrt 2003 and Howell&Mahrt 1997. With uncertainty estimation.
"""
function mrd(data_a::Vector, data_b::Vector, M::Integer, Mx::Integer)
    D = zeros(Float64, M - Mx)
    Dstd = similar(D)
    data_a2 = copy(data_a)
    data_b2 = copy(data_b)
    for ims in 0:(M-Mx)
        ms = M - ims
        l = 2^ms
        nw = round(Int, (2^M) / l)
        wmeans_a = zeros(Float64, nw)
        wmeans_b = similar(wmeans_a)
        for i in 1:nw
            k = round(Int, (i - 1) * l + 1)
            wmeans_a[i] = mean(data_a2[k:i*l])
            wmeans_b[i] = mean(data_b2[k:i*l])
            data_a2[k:i*l] .-= wmeans_a[i]
            data_b2[k:i*l] .-= wmeans_b[i]
        end
        if nw > 1
            D[ms+1] = mean(wmeans_a .* wmeans_b)
            Dstd[ms+1] = std(wmeans_a .* wmeans_b, mean=D[ms+1])
        end
    end
    return D, Dstd
end

"""
    completemrd(data::DataFrame, col1::String, col2::String, M::Integer, shift::Integer; normed::Bool=false)

performing a continuous MRD of wθ covariance for a DataFrame,
requiring fields 'time', 'col1'(e.g.='w'), and 'col2'(e.g.= 'T'). 'shift'(=M) allows for continuous MRD
"""
function completemrd(data::DataFrame, col1::String, col2::String, M::Integer, shift::Integer; normed::Bool=false)
    @info("MRD for DataFrame")
    #calculation of timestep as diff between second and first entry
    timestep = Millisecond(data.time[2] - data.time[1])
    #check with middle and middle+1 entry
    checktime = Millisecond(data.time[round(Int, size(data, 1) / 2)+1] - data.time[round(Int, size(data, 1) / 2)])
    if timestep != checktime
        @warn("Timestep and check-timestep do not agree! Careful!")
    end
    blocklength = 2^M * timestep
    println("M=", M, ", according to blocklength of ", canonicalize(Dates.CompoundPeriod(blocklength)))

    #shift between subsequent MRDs
    timeshift = shift * timestep

    #gaps in data
    gaps = detectgaps(data, Second(10))
    push!(gaps, [size(data, 1) + 1, data.time[end] + Second(1), Second(99)])
    @show size(gaps, 1)

    fx = ones(Float64, size(data, 1))
    if normed
        fx_str = string(col1, col2)
        fx_str_ret = string(col2, col1)
        tmp_fx = turb.turbflux(data, timestep * 2^11)
        if in(fx_str, names(tmp_fx))
            fx = gen.movingaverage(tmp_fx[:, fx_str], 2^M)
        elseif in(fx_str_ret, names(tmp_fx))
            fx = gen.movingaverage(tmp_fx[:, fx_str_ret], 2^M)
        else
            @warn("Norming: Fluxes not found in fluxfile! Applying no norm!")
        end
    end

    mrd_x = Millisecond.(Dates.value.(2 .^ collect(1:M) .* timestep))
    data_cont_mrd = Array{Float64}(undef, length(mrd_x), 0)
    time_middle = Vector{DateTime}(undef, 0)

    startidx = 1
    endidx = 2^M
    gapidx = 1
    nrblocks = 0
    normfct = 1.0
    p = ProgressUnknown("MRD blocks:")
    while endidx + shift <= size(data, 1)
        if nrblocks != 0
            startidx += shift
            endidx += shift
        end
        if endidx <= gaps.idx_before_gap[gapidx]
            datatouse1 = data[startidx:endidx, col1]
            datatouse2 = data[startidx:endidx, col2]
            time_middle_idx = data.time[startidx+round(Int, (endidx - startidx) / 2)]#(2^(M-1)).+nrblocks*shift]
            push!(time_middle, time_middle_idx)
            #calculate MRD of w and T
            (mrd_data_tmp, ) = mrd(datatouse1, datatouse2, M, 0)
            if normed
                normfct = sum(mrd_data_tmp[1:11])/fx[startidx+round(Int, (endidx - startidx) / 2)]
                mrd_data_tmp ./= normfct
            end
            data_cont_mrd = cat(data_cont_mrd, mrd_data_tmp, dims=2)
            nrblocks += 1
        else
            startidx = gaps.idx_before_gap[gapidx] + 1 - shift
            endidx = startidx + (2^M) - 1
            theonrstepsyet = round(Int, 1 + ((data.time[startidx+shift] - data.time[1] - blocklength).value / timeshift.value))
            nrstepsyet = length(time_middle)
            @show theonrstepsyet
            @show nrstepsyet
            nrstepstodo = theonrstepsyet - nrstepsyet
            for istep in 1:nrstepstodo
                if istep == 1 && size(time_middle, 1) == 0
                    push!(time_middle, data.time[1] + timeshift)
                else
                    push!(time_middle, time_middle[end] + timeshift)
                end
                data_cont_mrd = cat(data_cont_mrd, fill(NaN, M), dims=2)
                nrblocks += 1
            end
            gapidx += 1
        end
        ProgressMeter.next!(p)
    end
    ProgressMeter.finish!(p)
    #@show nrblocks
    return mrd_x, data_cont_mrd, time_middle
end

"""
    mrdpp(mrd_x::Vector, mrd::Array)

Postprocessing of 'completemrd' to obtain quantiles etc.
"""
function mrdpp(mrd_x::Vector, mrd::Array)
    mrd_D_median = fill(NaN, size(mrd, 1))
    mrd_D_max = similar(mrd_D_median)
    mrd_D_min = similar(mrd_D_median)
    mrd_D_quantile1 = similar(mrd_D_median)
    mrd_D_quantile3 = similar(mrd_D_median)

    for irow in axes(mrd, 1)
        mrd_D_median[irow] = median(filter(!isnan, mrd[irow, :]))
        mrd_D_min[irow] = minimum(filter(!isnan, mrd[irow, :]))
        mrd_D_max[irow] = maximum(filter(!isnan, mrd[irow, :]))
        quant = quantile(filter(!isnan, mrd[irow, :]), [0.25, 0.75])
        mrd_D_quantile1[irow] = quant[1]
        mrd_D_quantile3[irow] = quant[2]
    end
    return mrd_x, mrd_D_median, mrd_D_min, mrd_D_max, mrd_D_quantile1, mrd_D_quantile3
end

"""
    completenomrd(data::DataFrame, col1::String, col2::String, M::Integer, shift::Integer, nrpoints::Int64; normed::Bool=false)

performing a continuous NO-MRD of wθ covariance for a DataFrame,
requiring fields 'time', 'col1'(e.g.='w'), and 'col2'(e.g.= 'T'). 'shift'(=M) allows for continuous MRD
"""
function completenomrd(data::DataFrame, col1::String, col2::String, M::Integer, shift::Integer, nrpoints::Int64; normed::Bool=false)
    @info("NO-MRD for DataFrame")
    #calculation of timestep as diff between second and first entry
    timestep = Millisecond(data.time[2] - data.time[1])
    #check with middle and middle+1 entry
    checktime = Millisecond(data.time[round(Int, size(data, 1) / 2)+1] - data.time[round(Int, size(data, 1) / 2)])
    if timestep != checktime
        @warn("Timestep and check-timestep do not agree! Careful!")
    end
    blocklength = 2^M * timestep
    println("M=", M, ", according to blocklength of ", canonicalize(Dates.CompoundPeriod(blocklength)))

    #shift between subsequent MRDs
    timeshift = shift * timestep

    #gaps in data
    gaps = detectgaps(data, Second(10))
    push!(gaps, [size(data, 1) + 1, data.time[end] + Second(1), Second(99)])
    @show size(gaps, 1)

    rowsize = 2*nrpoints
    mrd_x = Array{Millisecond}(undef, rowsize, 0)
    data_cont_mrd = Array{Float64}(undef, rowsize, 0)
    time_middle = Vector{DateTime}(undef, 0)
    maxsizerow = 0

    fx = ones(Float64, size(data, 1))
    if normed
        fx_str = string(col1, col2)
        fx_str_ret = string(col2, col1)
        tmp_fx = turb.turbflux(data, timestep * 2^11)
        if in(fx_str, names(tmp_fx))
            fx = gen.movingaverage(tmp_fx[:, fx_str], 2^M)
        elseif in(fx_str_ret, names(tmp_fx))
            fx = gen.movingaverage(tmp_fx[:, fx_str_ret], 2^M)
        else
            @warn("Norming: Fluxes not found in fluxfile! Applying no norm!")
        end
    end

    startidx = 1
    endidx = 2^M
    gapidx = 1
    nrblocks = 0
    normpart = 1
    p = ProgressUnknown("MRD blocks:")
    while (endidx + shift <= size(data, 1) || nrblocks == 0)
        if nrblocks != 0
            startidx += shift
            endidx += shift
        end
        if endidx <= gaps.idx_before_gap[gapidx]
            datatouse1 = data[startidx:endidx, col1]
            datatouse2 = data[startidx:endidx, col2]
            time_middle_idx = data.time[startidx+round(Int, (endidx - startidx) / 2)]#(2^(M-1)).+nrblocks*shift]
            push!(time_middle, time_middle_idx)
            #calculate MRD of w and T
            (x_tmp, mrd_data_tmp, normpart) = contnomrd(datatouse1, datatouse2, nrpoints, normed)
            if normed
                normfct = normpart/fx[startidx+round(Int, (endidx - startidx) / 2)]
                mrd_data_tmp ./= normfct
            end            
            maxsizerow = maximum([maxsizerow, size(x_tmp, 1)])
            data_cont_mrd = cat(data_cont_mrd, vcat(mrd_data_tmp, zeros(Float64, rowsize-length(mrd_data_tmp))), dims=2)
            mrd_x = cat(mrd_x, vcat(x_tmp,fill(Millisecond(-99), rowsize-length(x_tmp))), dims=2)
            nrblocks += 1
        else
            startidx = gaps.idx_before_gap[gapidx] + 1 - shift
            endidx = startidx + (2^M) - 1
            theonrstepsyet = round(Int, 1 + ((data.time[startidx+shift] - data.time[1] - blocklength).value / timeshift.value))
            nrstepsyet = length(time_middle)
            @show theonrstepsyet
            @show nrstepsyet
            nrstepstodo = theonrstepsyet - nrstepsyet
            for istep in 1:nrstepstodo
                if istep == 1 && size(time_middle, 1) == 0
                    push!(time_middle, data.time[1] + timeshift)
                else
                    push!(time_middle, time_middle[end] + timeshift)
                end
                data_cont_mrd = cat(data_cont_mrd, fill(NaN, M), dims=2)
                mrd_x = cat(mrd_x, mrd_x[:,end])
                nrblocks += 1
            end
            gapidx += 1
        end
        ProgressMeter.next!(p)
    end
    ProgressMeter.finish!(p)
    #@show nrblocks
    return mrd_x[1:maxsizerow, :], data_cont_mrd[1:maxsizerow, :], time_middle
end

"""
    mrd_std(a::Vector, b::Vector, M::Integer)

MRD including std on different pairs of scales (fig.4 in Nilsson2014).
Code from Iva Stiperski
"""
function mrd_std(a::Vector, b::Vector, M::Integer)
    D = fill(NaN, M, M)
    if length(a) !== 2^M
        @error("Data length (a) must be 2^M!")
        return D
    elseif length(b) !== 2^M
        @error("Data length (b) must be 2^M!")
        return D
    end

    a_save = copy(a)
    b_save = copy(b)
    norm_std = std((a .- mean(a)).*(b .- mean(b)))

    for order in 1:2
        if order == 1 # compute non-diagonal terms decomposing b in smaller pieces than a
            a = copy(a_save)
            b = copy(b_save)
        elseif order == 2 # compute non-diagonal terms decomposing a in smaller pieces than b
            a = copy(b_save)
            b = copy(a_save)
        end

        for imsa in 0:M
            msa = M - imsa;
            la = 2^msa
            nwa=(2^M)/la

            b = copy(b_save)
            sumab = zeros(M)

            #case where scale of b is larger than scale of a
            #we substract the mean of b and don't store the mean
            for imsb in 0:M-msa-1
                msb = M - imsb
                lb = 2^msb
                nwb = (2^M)/lb

                for ib in 1:nwb
                    kb = Int((ib-1)*lb+1)
                    zb = b[kb]
                    for jb::Int in kb+1:kb+lb-1
                        zb = zb + b[jb]
                    end
                    zb = zb/lb
                    for jb::Int in kb:ib*lb
                        b[jb] = b[jb] - zb
                    end
                end
            end

            for ia in 1:nwa
                ka = Int((ia-1)*la+1)
                za = a[ka]
                for ja::Int in ka+1:ka+la-1
                    za = za + a[ja]
                end
                za = za/la
                for ja::Int in ka:ia*la
                    a[ja] = a[ja] - za
                end

                #case where scale of b is smaller than scale of a. we compute
                #the sume of the products of the mean
                for imsb in M-msa:M
                    msb = M-imsb
                    lb = 2^msb
                    nwb = (2^M)/lb

                    for ib in 1:2^(msa-msb)
                        kb = Int((ia - 1)*la + (ib - 1)*lb + 1)
                        zb = b[kb]
                        for jb::Int in kb+1:kb+lb-1
                            zb = zb + b[jb]
                        end
                        zb = zb/lb
                        if msb<M
                            sumab[msb+1] = sumab[msb+1]+(za*zb)^2
                        end
                        for jb::Int in kb:kb+lb-1
                            b[jb] = b[jb] - zb
                        end
                    end
                end
            end

            for imsb in M-msa:M
                msb=M-imsb
                lb = 2^msb
                nwb = (2^M)/lb
                if nwb>1 && nwa>1
                    if order == 1
                        D[msa+1, msb+1] = sqrt(sumab[msb+1]/nwb)/norm_std
                    elseif order == 2
                        D[msb+1, msa+1] = sqrt(sumab[msb+1]/nwb)/norm_std
                    end
                end
            end
        end
    end
    return D
end

"""
    mrd_std_block(data_tmp::DataFrame, starttime::DateTime, endtime::DateTime, M::Int, col1::String, col2::String, exp_col2::Int)::Array{Float64, 2}

Apply the different-scale-MRD (Nilsson14) to a block of data
and return an averaged array D.
"""
function mrd_std_block(data_tmp::DataFrame, starttime::DateTime, endtime::DateTime, M::Int, col1::String, col2::String, exp_col2::Int)::Array{Float64, 2}

    startidx = findfirst(data_tmp.time .>= starttime)
    endidx = findlast(data_tmp.time .<= endtime)

    tmpstart = startidx
    tmpend = startidx + 2^M - 1
    nrmrds = 0
    D = zeros(Float64, M, M)

    while tmpend <= endidx
        data_touse = copy(data_tmp[tmpstart:tmpend, :])
        turb.drdf!(data_touse, periodwise = false)
        Dtmp = mrd_std(data_touse[:,col1], data_touse[:, col2] .^exp_col2, M)
        D += Dtmp
        nrmrds += 1
        tmpstart += 2^M
        tmpend += 2^M
    end

    println("Different frequencies MRD done for ", nrmrds, " blocks.")
    notdone = Dates.canonicalize((endidx - tmpstart)*Millisecond(50))
    println("The last ", notdone, " of the block were not used.")

    D ./= nrmrds

    return D
end

"""
    nrpoints2modesformrd(datalen::Integer, nrpoints::Integer)

Helper function for nomrd/fastnomrd to get modes. Not exported.
"""
function nrpoints2modesformrd(datalen::Integer, nrpoints::Integer)

    if nrpoints < 3
        nrpoints = 3
        @warn("Nr. of Points must be >=3! Setting it to 3.")
    end

    if mod(datalen, 2) != 0 #make lengths even
        datalen -= 1
    end

    nroffset = -1
    ctmp = collect(range(0, log2(datalen / 2), nrpoints + nroffset))
    ctmp = vcat(round.(Int, 2 .^ ctmp), [datalen / 2])
    unique!(ctmp)

    len = length(ctmp) - 1
    while len < nrpoints
        nroffset += 1
        ctmp = collect(range(0, log2(datalen / 2), nrpoints + nroffset))
        ctmp = vcat(round.(Int, 2 .^ ctmp), [datalen / 2])
        unique!(ctmp)
        len = length(ctmp) - 1
    end
    return ctmp[end:-1:1]
end

"""
    generateseries(start::Int64, datlength::Int64)

Generate series for multiple (so continuous) no-mrds
"""
function generateseries(start::Int64, datlength::Int64)
    series = zeros(Int64, ceil(Int, log2(datlength))+1)
    a = start
    idx = 1
    while a <= datlength
        series[idx] = a
        a *= 2
        idx += 1
    end
    return series[idx-1:-1:1]
end

"""
    fastnomrd(data_a_in::Vector, data_b_in::Vector, series::Vector{Int64}, nanthresh=0.8)

NO-MRD without calculating quantiles.
"""
function fastnomrd(data_a_in::Vector, data_b_in::Vector, series::Vector{Int64}, nanthresh=0.8)

    data_a = copy(data_a_in)
    data_b = copy(data_b_in)
    timestep = Millisecond(50) #only needed for output

    if length(data_a) == 0 || length(data_b) == 0
        return [0], [0], [0], [0], [0], [0], [0]
    end

    #println("Non-orthogonal MRD in ", length(kdecr), " resolution modes...")
    D = zeros(Float64, length(series) - 1)
    wmeans_a = fill(NaN, length(data_a))
    wmeans_b = fill(NaN, length(data_b))
    last_sum_a = 0.0
    last_sum_b = 0.0
    nrnans = 0
    nrnanwdws = 0
    #tmp = fill(NaN, size(data_a, 1))
    for (kidx, k) in enumerate(series)
        nw = length(data_a) - k + 1
        fill!(wmeans_a, NaN)
        fill!(wmeans_b, NaN)
        nrnanwdws = 0
        last_sum_a = 0.0
        last_sum_b = 0.0
        for i in 1:nw
            if i == 1
                idcs = 1:k
                nrnans = count(x -> isnan(x), data_a[idcs])
                if nrnans < k
                    last_sum_a = sum(filter(!isnan, data_a[idcs]))
                    last_sum_b = sum(filter(!isnan, data_b[idcs]))
                    wmeans_a[1] = last_sum_a / (k - nrnans)
                    wmeans_b[1] = last_sum_b / (k - nrnans)
                else
                    wmeans_a[1] = NaN
                    wmeans_b[1] = NaN
                    nrnanwdws += 1
                    last_sum_a = 0.0
                    last_sum_b = 0.0
                end
            else
                idx = k + i - 1
                if isnan(data_a[idx]) || isnan(data_b[idx])
                    nrnans += 1
                    corrforward_a = 0
                    corrforward_b = 0
                else
                    corrforward_a = data_a[idx]
                    corrforward_b = data_b[idx]
                end
                if isnan(data_a[i-1]) || isnan(data_b[i-1])
                    nrnans -= 1
                    corrbackward_a = 0
                    corrbackward_b = 0
                else
                    corrbackward_a = data_a[i-1]
                    corrbackward_b = data_b[i-1]
                end

                if nrnans / k < nanthresh
                    last_sum_a += (corrforward_a - corrbackward_a)
                    last_sum_b += (corrforward_b - corrbackward_b)
                    wmeans_a[i] = last_sum_a / (k - nrnans)
                    wmeans_b[i] = last_sum_b / (k - nrnans)
                else
                    last_sum_a += (corrforward_a - corrbackward_a)
                    last_sum_b += (corrforward_b - corrbackward_b)
                    wmeans_a[i] = NaN
                    wmeans_b[i] = NaN
                    nrnanwdws += 1
                end
            end
        end

        #write to result vectors
        #=if kidx == 8
            tmp[1:nw] = wmeans_a[1:nw] .* wmeans_b[1:nw]
        end=#
        if kidx > 1
            if nrnanwdws/nw < nanthresh
                D[length(series)-kidx+1] = mean(filter(!isnan, wmeans_a[1:nw] .* wmeans_b[1:nw]))
            else
                D[length(series)-kidx+1] = NaN
            end
        end


        #subtract the calculated means preparing the next step (except last step)
        if kidx < length(series)
            data_a[1:nw] .-= wmeans_a[1:nw]
            data_a[nw+1:end] .-= wmeans_a[nw]
            data_b[1:nw] .-= wmeans_b[1:nw]
            data_b[nw+1:end] .-= wmeans_b[nw]
        end
    end

    x_return = series[end-1:-1:1] .* timestep

    return x_return, D
end

"""
    fastnomrd(data_a_in, series::Vector{Int64}, nanthresh=0.8)

NO-MRD without calculating quantilesfor one variable.
"""
function fastnomrd(data_a_in, series::Vector{Int64}, nanthresh=0.8)

    data_a = copy(data_a_in)
    timestep = Millisecond(50) #only needed for output

    if length(data_a) == 0
        return [0], [0], [0], [0], [0], [0], [0]
    end

    #println("Non-orthogonal MRD in ", length(kdecr), " resolution modes...")
    D = zeros(Float64, length(series) - 1)
    wmeans_a = fill(NaN, length(data_a))
    last_sum_a = 0.0
    nrnans = 0
    nrnanwdws = 0
    #tmp = fill(NaN, size(data_a, 1))
    for (kidx, k) in enumerate(series)
        nw = length(data_a) - k + 1
        fill!(wmeans_a, NaN)
        nrnanwdws = 0
        last_sum_a = 0.0
        for i in 1:nw
            if i == 1
                idcs = 1:k
                nrnans = count(x -> isnan(x), data_a[idcs])
                if nrnans < k
                    last_sum_a = sum(filter(!isnan, data_a[idcs]))
                    wmeans_a[1] = last_sum_a / (k - nrnans)
                else
                    wmeans_a[1] = NaN
                    nrnanwdws += 1
                    last_sum_a = 0.0
                end
            else
                idx = k + i - 1
                if isnan(data_a[idx])
                    nrnans += 1
                    corrforward_a = 0
                else
                    corrforward_a = data_a[idx]
                end
                if isnan(data_a[i-1])
                    nrnans -= 1
                    corrbackward_a = 0
                else
                    corrbackward_a = data_a[i-1]
                end

                if nrnans / k < nanthresh
                    last_sum_a += (corrforward_a - corrbackward_a)
                    wmeans_a[i] = last_sum_a / (k - nrnans)
                else
                    last_sum_a += (corrforward_a - corrbackward_a)
                    wmeans_a[i] = NaN
                    nrnanwdws += 1
                end
            end
        end

        #write to result vectors
        #=if kidx == 8
            tmp[1:nw] = wmeans_a[1:nw] .* wmeans_b[1:nw]
        end=#
        if kidx > 1
            if nrnanwdws/nw < nanthresh
                D[length(series)-kidx+1] = mean(filter(!isnan, wmeans_a[1:nw]))
            else
                D[length(series)-kidx+1] = NaN
            end
        end


        #subtract the calculated means preparing the next step (except last step)
        if kidx < length(series)
            data_a[1:nw] .-= wmeans_a[1:nw]
            data_a[nw+1:end] .-= wmeans_a[nw]
        end
    end

    x_return = series[end-1:-1:1] .* timestep

    return x_return, D
end

"""
    nomrd(data_a_in::Vector, data_b_in::Vector, series::Vector{Int64}, nanthresh=0.5)

Non-orthogonal Multiresolution Flux Decomposition.
See Howell&Mahrt 1997. GLG (9) and comment few lines further.
"""
function nomrd(data_a_in::Vector, data_b_in::Vector, series::Vector{Int64}, nanthresh=0.5)

    data_a = copy(data_a_in)
    data_b = copy(data_b_in)
    timestep = Millisecond(50) #only needed for output

    if length(data_a) == 0 || length(data_b) == 0
        return [0], [0], [0], [0], [0], [0], [0]
    end

    #println("Non-orthogonal MRD in ", length(kdecr), " resolution modes...")
    D = zeros(Float64, length(series) - 1)
    Dq1 = similar(D)    #first quartile
    Dq3 = similar(D)   #third quartile
    wmeans_a = fill(NaN, length(data_a))
    wmeans_b = fill(NaN, length(data_b))
    last_sum_a = 0.0
    last_sum_b = 0.0
    nrnans = 0
    for (kidx, k) in enumerate(series)
        nw = length(data_a) - k + 1
        fill!(wmeans_a, NaN)
        fill!(wmeans_b, NaN)
        last_sum_a = 0.0
        last_sum_b = 0.0
        for i in 1:nw
            if i == 1
                idcs = 1:k
                nrnans = count(x -> isnan(x), data_a[idcs])
                if nrnans < k
                    last_sum_a = sum(filter(!isnan, data_a[idcs]))
                    last_sum_b = sum(filter(!isnan, data_b[idcs]))
                    wmeans_a[1] = last_sum_a / (k - nrnans)
                    wmeans_b[1] = last_sum_b / (k - nrnans)
                else
                    wmeans_a[1] = NaN
                    wmeans_b[1] = NaN
                    last_sum_a = 0.0
                    last_sum_b = 0.0
                end
            else
                idx = k + i - 1
                if isnan(data_a[idx]) || isnan(data_b[idx])
                    nrnans += 1
                    corrforward_a = 0
                    corrforward_b = 0
                else
                    corrforward_a = data_a[idx]
                    corrforward_b = data_b[idx]
                end
                if isnan(data_a[i-1]) || isnan(data_b[i-1])
                    nrnans -= 1
                    corrbackward_a = 0
                    corrbackward_b = 0
                else
                    corrbackward_a = data_a[i-1]
                    corrbackward_b = data_b[i-1]
                end

                if nrnans / k < nanthresh
                    last_sum_a += (corrforward_a - corrbackward_a)
                    last_sum_b += (corrforward_b - corrbackward_b)
                    wmeans_a[i] = last_sum_a / (k - nrnans)
                    wmeans_b[i] = last_sum_b / (k - nrnans)
                else
                    last_sum_a += (corrforward_a - corrbackward_a)
                    last_sum_b += (corrforward_b - corrbackward_b)
                    wmeans_a[i] = NaN
                    wmeans_b[i] = NaN
                end
            end
        end

        #write to result vectors
        if kidx > 1
            tmp = 1
            try
                tmp = quantile(filter(!isnan, wmeans_a[1:nw] .* wmeans_b[1:nw]), [0.25, 0.5, 0.75])
            catch e
                tmp = [NaN, NaN, NaN]
            end
            D[length(series)-kidx+1] = tmp[2]
            Dq1[length(series)-kidx+1] = tmp[1]
            Dq3[length(series)-kidx+1] = tmp[3]
        end

        #subtract the calculated means preparing the next step (except last step)
        if kidx < length(series)
            data_a[1:nw] .-= wmeans_a[1:nw]
            data_b[1:nw] .-= wmeans_b[1:nw]
            data_a[nw+1:end] .-= wmeans_a[nw]
            data_b[nw+1:end] .-= wmeans_b[nw]
        end
    end

    x_return = series[end-1:-1:1] .* timestep

    #normalization
    #binwidth = ones(Float64, length(x_return))
    #=for i in 1:length(binwidth)-1 
        binwidth[i] = Dates.value(kdecr_return[i+1]-kdecr_return[i])*1000
    end
    binwidth[end] = Dates.value(size(data_a2, 1)*Millisecond(50)-kdecr_return[end])*1000
    =#
    weightedsum = 1#abs.(sum(binwidth .* D))
    D_norm = D ./ weightedsum
    Dq1_norm = Dq1 ./ weightedsum
    Dq3_norm = Dq3 ./ weightedsum

    return x_return, D_norm, Dq1_norm, Dq3_norm, D, Dq1, Dq3
end

"""
    contnomrd(data_a::Vector, data_b::Vector, nrsteps::Int64, normed::Bool=false)

Multiple no-MRDs to obtain a quasi-continuous no-MRD
"""
function contnomrd(data_a::Vector, data_b::Vector, nrsteps::Int64, normed::Bool=false)
    x_out = Vector{Millisecond}(undef, 0)
    D = Vector{Float64}(undef, 0)

    startidx = 1
    normpart = 1
    while length(x_out) < nrsteps
        ser = generateseries(startidx, size(data_a, 1))
        if length(ser) < 2
            break
        end
        (tmp_x, tmp_D) = fastnomrd(data_a, data_b, ser)
        if normed &&startidx == 1
            normpart = sum(tmp_D[1:11])
        end
        x_out = vcat(x_out, tmp_x)
        D = vcat(D, tmp_D)

        startidx += 2
    end

    #sort the output
    perm = sortperm(x_out)
    x_out = x_out[perm]
    D = D[perm]

    if count(x->isnan(x), D) >= 2
        D .= NaN
    end

    return x_out, D, normpart
end

"""
    contnomrd(data_a, nrsteps::Int64, normed::Bool=false)

Multiple no-MRDs to obtain a quasi-continuous no-MRD for one variable
"""
function contnomrd(data_a, nrsteps::Int64, normed::Bool=false)
    x_out = Vector{Millisecond}(undef, 0)
    D = Vector{Float64}(undef, 0)

    startidx = 1
    normpart = 1
    while length(x_out) < nrsteps
        ser = generateseries(startidx, size(data_a, 1))
        if length(ser) < 2
            break
        end
        (tmp_x, tmp_D) = fastnomrd(data_a, ser)
        if normed &&startidx == 1
            normpart = sum(tmp_D[1:11])
        end
        x_out = vcat(x_out, tmp_x)
        D = vcat(D, tmp_D)

        startidx += 2
    end

    #sort the output
    perm = sortperm(x_out)
    x_out = x_out[perm]
    D = D[perm]

    if count(x->isnan(x), D) >= 2
        D .= NaN
    end

    return x_out, D, normpart
end

"""
    catdatanomrd(data::DataFrame, col1::String, col2::String, nrpoints::Integer)

Apply Non-Orthogonal MRD to categorized data. Run per block.
"""
function catdatanomrd(data::DataFrame, col1::String, col2::String, nrpoints::Integer)

    if size(data, 1) == 0
        return [Millisecond(0)], [0], [0], [0]
    end

    timestep = Millisecond(50)

    #gaps in data
    gaps = detectgaps(data, Second(1))
    #artificially append last row
    push!(gaps, [size(data, 1), data.time[end] + Second(1), Second(99)])

    #find out longest & shortest block without gap
    maxblocklen = gaps.idx_before_gap[1]
    minblocklen = gaps.idx_before_gap[1]
    for i in axes(gaps, 1)[2:end]
        maxblocklen = maximum([maxblocklen, gaps.idx_before_gap[i] - gaps.idx_before_gap[i-1] - 1])
        minblocklen = minimum([minblocklen, gaps.idx_before_gap[i] - gaps.idx_before_gap[i-1] - 1])
    end

    if isodd(maxblocklen)
        maxblocklen -= 1
    end
    if isodd(minblocklen)
        minblocklen -= 1
    end

    nomrd_x = collect(timestep:timestep:maxblocklen*timestep)
    mrd_Dn = Array{Float64}(undef, length(nomrd_x), 0)
    q1 = similar(mrd_Dn)
    q3 = similar(mrd_Dn)

    if maxblocklen - minblocklen == 0
        blocklen = maxblocklen
        @show blocklen

        sizetoassign = round(Int, nrpoints * 1.2)
        nomrd_x = zeros(DateTime, sizetoassign)
        mrd_Dn = Array{Float64}(undef, sizetoassign, 0)
        q1 = similar(mrd_Dn)
        q3 = similar(mrd_Dn)

        @showprogress "Calculating NO-MRDs for blocks..." for gi in axes(gaps, 1)

            startidx = 1
            endidx = gaps.idx_before_gap[gi]
            if gi > 1
                startidx = gaps.idx_before_gap[gi-1] + 1
            end

            (t, dn_tmp) = contnomrd(
                data[startidx:endidx, col1],
                data[startidx:endidx, col2], nrpoints)

            dn = fill(NaN, size(mrd_Dn, 1))
            #q1_t = similar(dn)
            #q3_t = similar(dn)
            dn[1:length(dn_tmp)] = dn_tmp
            #q1_t[1:length(dq1_tmp)] = dq1_tmp
            #q3_t[1:length(dn_tmp)] = dq3_tmp


            if gi == 1
                nomrd_x[1:length(t)] = t
            end
            mrd_Dn = cat(mrd_Dn, dn, dims=2)
            #q1 = cat(q1, q1_t, dims=2)
            #q3 = cat(q3, q3_t, dims=2)
        end

    else #linear interpolation to come to common grid
        println("using interpolation to come to common x-points")
        @show maxblocklen
        @show minblocklen

        @showprogress "Calculating NO-MRDs for blocks..." for gi in axes(gaps, 1)
            startidx = 1
            endidx = gaps.idx_before_gap[gi]
            if gi > 1
                startidx = gaps.idx_before_gap[gi-1] + 1
            end

            if endidx - startidx <= 0.2 * maxblocklen
                continue
            end

            (t, dn) = contnomrd(
                data[startidx:endidx, col1],
                data[startidx:endidx, col2], nrpoints)

            if length(t) == 0
                continue
            end

            itpn = interpolate((Dates.value.(t),), dn, Gridded(Linear()))
            #itpq1 = interpolate((Dates.value.(t),), dq1n, Gridded(Linear()))
            #itpq3 = interpolate((Dates.value.(t),), dq3n, Gridded(Linear()))


            tmp_dn = fill(NaN, length(nomrd_x))
            #tmp_q1 = similar(tmp_dn)
            #tmp_q3 = similar(tmp_dn)
            t_idx1 = findfirst(x -> x == t[1], nomrd_x)
            t_idxend = findfirst(x -> x == t[end], nomrd_x)
            tmp_dn[t_idx1:t_idxend] = itpn.(Dates.value.(collect(t[1]:timestep:t[end])))
            #tmp_q1[t_idx1:t_idxend] = itpq1.(Dates.value.(collect(t[1]:timestep:t[end])))
            #tmp_q3[t_idx1:t_idxend] = itpq3.(Dates.value.(collect(t[1]:timestep:t[end])))

            mrd_Dn = cat(mrd_Dn, tmp_dn, dims=2)
            #q1 = cat(q1, tmp_q1, dims=2)
            #q3 = cat(q3, tmp_q3, dims=2)

        end
    end

    #statistics
    dn_median = zeros(Float64, size(mrd_Dn, 1))
    dn_q1 = similar(dn_median)
    dn_q3 = similar(dn_median)

    for jdx in axes(dn_median, 1)
        ana_data = filter(!isnan, mrd_Dn[jdx, :])
        if size(ana_data, 1) > 1
            (dn_q1[jdx], dn_median[jdx], dn_q3[jdx]) = quantile(ana_data, [0.25, 0.5, 0.75])
        elseif size(ana_data, 1) == 1 && size(q1, 2) > 0
            #dn_q1[jdx] = q1[jdx]
            #dn_q3[jdx] = q3[jdx]
            dn_median[jdx] = ana_data[1]
        else
            (dn_q1[jdx], dn_median[jdx], dn_q3[jdx]) = (NaN, NaN, NaN)
        end
    end
    return nomrd_x, dn_median, dn_q1, dn_q3
end

"""
    calc_mrds(res_rm::RemoteChannel, actcol::Int64, data::AbstractArray, data_time::AbstractVector, nrpoints::Int)

Called on single workers.
"""
function calc_mrds(res_rm::RemoteChannel, actcol::Int64, data::AbstractArray, data_time::AbstractVector, nrpoints::Int)
    time_mid = data_time[round(Int, (size(data_time, 1)) / 2)]
    if size(data, 2) == 2
        (x_tmp, mrd_data_tmp, a) = contnomrd(data[:, 1], data[:, 2], nrpoints)
    elseif size(data, 2) == 1
        (x_tmp, mrd_data_tmp, a) = contnomrd(data, nrpoints)
    end
    if myid() == 2
        println(actcol)
    end
    put!(res_rm, (actcol, size(x_tmp, 1), time_mid, x_tmp, mrd_data_tmp))
end

"performing a quasi-continuous 2D non-orthogonal MRD for a DataFrame.
'len' data length of single MRDs as fraction of size(data, 1), 'shift'(as fraction of len) 0<shift<=1"
function nomrd2d(data::Array, data_time::Vector, meanwind::Vector{Float64}, nrpoints::Integer, len::Float64, shift::Float64; drperblock::Bool=false)
    @info("2D NO-MRD")

    #check for suitable values
    if !(0 < shift <= 1)
        @warn("0 < shift <= 1 not fulfilled. Setting shift to 0.1")
        shift = 0.1
    end

    if !(0 < len <= 1)
        @warn("0 < len <= len not fulfilled. Setting len = 0.05")
        len = 0.05
    end

    #calculation of timestep as diff between second and first entry
    timestep = Millisecond(data_time[2] - data_time[1])
    #check with middle and middle+1 entry
    checktime = Millisecond(data_time[round(Int, size(data, 1) / 2)+1] - data_time[round(Int, size(data, 1) / 2)])
    if timestep != checktime
        @warn("Timestep and check-timestep do not agree! Careful!")
    end
    blocklength = round(Int, len * size(data, 1)) * timestep
    blocklen_idx = round(Int, len * size(data, 1))
    println("length=", round(len * 100, digits=2), "% of data length, according to blocklength of ", canonicalize(Dates.CompoundPeriod(blocklength)))

    #shift between subsequent MRDs
    shift_idx = round(Int, shift * blocklen_idx)

    #gaps in data
    gaps = detectgapsdatatimearray(data_time, Second(10))
    @show size(gaps, 1)
    push!(gaps, [size(data, 1) + 1, data_time[end] + Second(1), Second(99)])

    #calculate nr. of no-mrds
    nrnomrds = ceil(Int, 1 + (1 - len) / (shift * len))

    #estimate maximum dimensions for arrays to allocate space
    sizerow = round(Int, 1.8 * nrpoints)
    sizecol = round(Int, 1.2 * (nrnomrds + 2 * (size(gaps, 1) + 2)))

    actualsizecol = 1 #Threads.Atomic{Int}(0)

    gapidx = 1
    gapencountered = false
    lastmrd = false

    startidcs = zeros(Int64, sizecol) #Vector{Int64}(undef, sizecol)
    endidcs = zeros(Int64, sizecol) #Vector{Int64}(undef, sizecol)

    startidcs[1] = 1
    endidcs[1] = blocklen_idx

    while !lastmrd

        actualsizecol += 1

        if !gapencountered
            startidcs[actualsizecol] = startidcs[actualsizecol-1] + shift_idx
            endidcs[actualsizecol] = endidcs[actualsizecol-1] + shift_idx
        else
            startidcs[actualsizecol] = gaps.idx_before_gap[gapidx-1] + 1
            endidcs[actualsizecol] = startidcs[actualsizecol] + shift_idx - 1
            gapencountered = false
        end

        if endidcs[actualsizecol] > size(data, 1)
            endidcs[actualsizecol] = size(data, 1)
            lastmrd = true
        end

        if endidcs[actualsizecol] > gaps.idx_before_gap[gapidx]
            endidcs[actualsizecol] = gaps.idx_before_gap[gapidx]
            gapencountered = true
            gapidx += 1
        end

        if endidcs[actualsizecol] - startidcs[actualsizecol] < 0.1 * blocklen_idx
            @warn("using super-short (<10% of normal) block. Maybe needs adjustment in code.")
        end

    end

    startidcs = startidcs[1:actualsizecol]
    endidcs = endidcs[1:actualsizecol]

    meanwindout = zeros(Float64, actualsizecol)
    for jx in 1:lastindex(meanwindout)
        meanwindout[jx] = mean(filter(!isnan, meanwind[startidcs[jx]:endidcs[jx]]))
    end

    @show actualsizecol
    println("Tasks per worker: ", round(Int, actualsizecol/length(workers())))
    println()

    res = Channel{Tuple{Int64,Int64,DateTime,Vector{DateTime},Vector{Float64}}}(actualsizecol)
    res_rm = RemoteChannel(() -> res)

    @sync @distributed for i in 1:actualsizecol
        if drperblock
            tmpdat = copy(data[startidcs[i]:endidcs[i], :])
            turb.drdf!(tmpdat, periodwise=false)
            data2 = Array{Float64,2}(undef, size(tmpdat, 1), length(cols_for_mrd))
            for (colidx, colstring) in enumerate(cols_for_mrd)
                data2[:, colidx] = tmpdat[:, colstring]
            end
            calc_mrds(res_rm, i, data2, data_time[startidcs[i]:endidcs[i]], nrpoints)
        else
            calc_mrds(res_rm, i, data[startidcs[i]:endidcs[i], :], data_time[startidcs[i]:endidcs[i]], nrpoints)
        end
    end

    #take from result-Channel and write to arrays
    #allocate Arrays
    mrd_x = Array{Millisecond,2}(undef, sizerow, actualsizecol)
    data_mrd = Array{Float64,2}(undef, sizerow, actualsizecol)
    time_middle = Vector{DateTime}(undef, actualsizecol)
    maxsizerow = Vector{Int64}(undef, actualsizecol)

    fill!(mrd_x, Millisecond(-9999))
    fill!(data_mrd, NaN)
    fill!(time_middle, DateTime(1990, 03, 28, 03, 14, 15))
    fill!(maxsizerow, 0)

    while !isempty(res)
        (idx, maxrow, time_mid, x_tmp, mrd_data_tmp) = take!(res)
        mrd_x[1:maxrow, idx] .= x_tmp
        data_mrd[1:maxrow, idx] .= mrd_data_tmp
        time_middle[idx] = time_mid
        maxsizerow[idx] = maxrow
    end

    maxrow = maximum(maxsizerow)
    #cut arrays to actual/maximum dimensions
    mrd_x = mrd_x[1:maxrow[], :] #1:actualsizecol]
    data_mrd = data_mrd[1:maxrow[], :] # 1:actualsizecol]

    return mrd_x, data_mrd, time_middle, meanwindout
end

"""
    write2dmrddatatonetcdf(x::AbstractArray, mrd_data::AbstractArray, time_middle::AbstractVector,
    meanwind::Vector{Float64}, evalstart::DateTime, evalend::DateTime, cols_for_mrd::Vector,
    nrpoints::Int64, blocklen::Float64, shift::Float64, sourcefile::String, target::String, deflatelvl::Int64=5)

Save the results from the distributed 2D-MRD-computation to netcdf-file
"""
function write2dmrddatatonetcdf(x::AbstractArray, mrd_data::AbstractArray, time_middle::AbstractVector,
    meanwind::Vector{Float64}, evalstart::DateTime, evalend::DateTime, cols_for_mrd::Vector,
    nrpoints::Int64, blocklen::Float64, shift::Float64, sourcefile::String, target::String, deflatelvl::Int64=5)
    ds = NCDataset(target, "c")
    defDim(ds, "single_mrds", size(time_middle, 1))
    defDim(ds, "mrd_timescale", size(x, 1))
    defDim(ds, "tmp", 1)
    defDim(ds, "nr_cols_for_mrd", size(cols_for_mrd, 1))
    defVar(ds, "time_middle", time_middle, ("single_mrds",); shuffle=true, deflatelevel=deflatelvl)
    xnc = defVar(ds, "x", Dates.value.(x), ("mrd_timescale", "single_mrds",); shuffle=true, deflatelevel=deflatelvl)
    xnc.attrib["units"] = "Millisecond"
    defVar(ds, "mrd_data", mrd_data, ("mrd_timescale", "single_mrds",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "meanwind", meanwind, ("single_mrds", ); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "evalstart", [evalstart], ("tmp",))
    defVar(ds, "evalend", [evalend], ("tmp",))
    defVar(ds, "cols_for_mrd", cols_for_mrd, ("nr_cols_for_mrd",))
    defVar(ds, "nrpoints", nrpoints, ())
    defVar(ds, "blocklen", blocklen, ())
    defVar(ds, "shift", shift, ())
    defVar(ds, "sourcefile", sourcefile, ())
    close(ds)
end

"""
    read2dmrdfromnetcdf(source::String)

Read 2dmrd-data from NetCDF4-file.
"""
function read2dmrdfromnetcdf(source::String)
    @info("Reading 2D-MRD data from NetCDF-file")
    ds = Dataset(source, "r")
    #load
    tmptime_mid = ds["time_middle"]
    time_middle = tmptime_mid[:]
    xtmp = ds["x"]
    x = Millisecond.(xtmp[:])
    mrd_data_tmp = ds["mrd_data"]
    mrd_data = mrd_data_tmp[:]
    meanwind_tmp = ds["meanwind"]
    meanwind = meanwind_tmp[:]
    evalsttmp = ds["evalstart"]
    evalstart = evalsttmp[:][1]
    evalendtmp = ds["evalend"]
    evalend = evalendtmp[:][1]
    colstmp = ds["cols_for_mrd"]
    cols = colstmp[:]
    nrpointstmp = ds["nrpoints"]
    nrpoints = nrpointstmp[1]
    blocklentmp = ds["blocklen"]
    blocklen = blocklentmp[1]
    shifttmp = ds["shift"]
    shift = shifttmp[1]
    sourcetmp = ds["sourcefile"]
    sourcefile = sourcetmp[1]
    close(ds)
    return time_middle, x, mrd_data, meanwind,
    evalstart, evalend, cols, nrpoints, blocklen,
    shift, sourcefile
end

"""
    instrumentnamefromsourcefile(sourcefile::String, regexp=r"[a-zA-Z0-9]*.nc")::String

Extract instrument name from source file
"""
function instrumentnamefromsourcefile(sourcefile::String, regexp=r"\/[a-zA-Z0-9]*\.nc")::String
    inst = match(regexp, sourcefile).match[2:end-3]
    if inst[end-1:end] == "df"
        inst = inst[1:end-2]
    end
    if inst[end-1:end] == "db"
        inst = inst[1:end-2]
    end
    return uppercase(inst)
end

"""
    create2dmrdtitle(sourcefile::String, evalstart::DateTime,
    evalend::DateTime)::String

Create title string for 2D-MRD plot
"""
function create2dmrdtitle(sourcefile::String, evalstart::DateTime,
    evalend::DateTime)::String
    instname = instrumentnamefromsourcefile(sourcefile)
    evalstartstring = Dates.format(evalstart, "dd.mm. HH:MM")
    if Date(evalstart) == Date(evalend)
        evalendstring = Dates.format(evalend, "HH:MM")
    else
        evalendstring = Dates.format(evalend, "dd.mm. HH:MM")
    end
    return string(instname, " ", evalstartstring, " \u2013 ", evalendstring)
end

end # module MRD