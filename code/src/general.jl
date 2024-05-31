######################################################
###                GENERAL MODULE                  ###
###            author: Michi Haugeneder            ###
######################################################
#=
Module containing general functions applicable to all measurements
=#

module gen

using CSV, DataFrames, Dates, StatsBase

export movingaverage, backlookingmovavg, getinstrumentheights

"Create moving average (michi 16.09.2021); omitting NaNs!!"
function movingaverage(X::Vector, numofele::Integer)
    BackDelta = div(numofele, 2)
    ForwardDelta = isodd(numofele) ? div(numofele, 2) : div(numofele, 2) - 1
    len = size(X, 1)
    #create vector with vec_isnan[i]=1 if X[i]=NaN, 0 otherwise
    vec_isnan = isnan.(X)
    if numofele >= len
        #println("#avg elements >= size(vector). Returning mean.")
        return ones(Float64, len) .* mean(filter(!isnan, X))
    else
        Y = ones(Float64, size(X, 1))
        firstnonnan = findfirst(x -> x == false, vec_isnan)
        len = findlast(x -> x == false, vec_isnan)
        if isnothing(firstnonnan)
            Y .= NaN
            return Y
        elseif firstnonnan > 1
            Y[1:firstnonnan-1] .= NaN
        end
        if len < length(X)
            Y[len:end] .= NaN
        end
        if (len - firstnonnan) +1 <= numofele
            return Y.* mean(filter(!isnan, X))
        end
        n = firstnonnan
        summed = sum(filter(!isnan, X[(0:ForwardDelta-1).+n]))
        curr_nans = sum(vec_isnan[(0:ForwardDelta-1).+n])
        curr_nonnans = ForwardDelta - curr_nans #count(x -> x == 0, vec_isnan[(0:ForwardDelta-1).+n])
        if n <= (BackDelta + firstnonnan)
            for n in (0:BackDelta) .+ firstnonnan
                #@info("Loop1")
                #if not NaN
                if !vec_isnan[n+ForwardDelta]
                    summed += X[n+ForwardDelta]
                    curr_nonnans += 1
                    Y[n] = summed / curr_nonnans
                    #if NaN
                else
                    curr_nans += 1
                    if n > 1
                        Y[n] = Y[n-1]
                    else
                        Y[n] = NaN
                    end
                end
            end
        end
        n = BackDelta - firstnonnan + 1
        for n in (BackDelta+firstnonnan+1):(len-ForwardDelta)
            #@info("Loop2")
            curr_nans += vec_isnan[n+ForwardDelta]
            curr_nans -= vec_isnan[n-BackDelta-1]
            curr_nonnans += !vec_isnan[n+ForwardDelta]
            curr_nonnans -= !vec_isnan[n-BackDelta-1]
            if !vec_isnan[n+ForwardDelta]
                summed += X[n+ForwardDelta]
            end
            if !vec_isnan[n-BackDelta-1]
                summed -= X[n-BackDelta-1]
            end
            if curr_nonnans > 0
                Y[n] = summed / curr_nonnans
            else
                Y[n] = NaN
            end
        end
        n = len - ForwardDelta
        for n in len-ForwardDelta+1:len
            #@info("Loop3")
            curr_nans -= vec_isnan[n-BackDelta-1]
            curr_nonnans -= !vec_isnan[n-BackDelta-1]
            if !vec_isnan[n-BackDelta-1]
                summed -= X[n-BackDelta-1]
            end
            if curr_nonnans > 0
                Y[n] = summed / curr_nonnans
            else
                Y[n] = NaN
            end
        end
    end
    return Y
end

"Only back-looking moving average."
function backlookingmovavg(X::Vector, numofele::Integer)::Vector
    len = size(X, 1)
    #create vector with vec_isnan[i]=1 if X[i]=NaN, 0 otherwise
    vec_isnan = isnan.(X)
    if numofele >= len
        println("#avg elements >= size vector! Returning original vector!")
        return X
    else
        Y = fill(NaN, size(X, 1))
        summed = 0
        curr_nans = 0
        n = 1
        while n <= numofele    #first part until numofele elements available
            #if not NaN
            if !vec_isnan[n]
                summed += X[n]
                Y[n] = summed / (n - curr_nans)
                #if NaN
            else
                curr_nans += 1
                if n > 1
                    Y[n] = Y[n-1]
                else
                    Y[n] = NaN
                end
            end
            n += 1
        end
        while n <= len    #second part until end (length of input vector)
            curr_nans += vec_isnan[n]
            curr_nans -= vec_isnan[n-numofele]
            if !vec_isnan[n]
                summed += X[n]
            end
            if !vec_isnan[n-numofele]
                summed -= X[n-numofele]
            end
            if numofele > curr_nans
                Y[n] = summed / (numofele - curr_nans)
            else
                Y[n] = NaN
            end
            #println(curr_nans)
            n += 1
        end
    end
    return Y
end
#=
"Get the height above snow surface for the instruments 
with a given source file. Interpolate between noted heights."
function getinstrumentheights(source::String, timevec::Vector)
    @info("Loading instrument heights. Using defined column names!
    Linearly interpolating between noted values")
    df = CSV.File(source; header=0, skipto=12, tasks = Threads.nthreads())|> Tables.matrix
    dateformat = DateFormat("yyyy-mm-dd HH:MM:SS")
    timenoted = DateTime.(df[:,1], dateformat)
    inivec = zeros(size(timevec, 1))
    outdf = DataFrame(time = timevec, t1tc1=inivec, t1tc2=inivec,
    t1tc3=inivec, t1tc4=inivec, t1tc5=inivec, t1irg=inivec, 
    t2tc1=inivec, t2tc2=inivec, t2tc3=inivec, t2tc4=inivec, 
    t2tc5=inivec, t2tc6=inivec, t2tc7=inivec, t2tc8=inivec,
    t2tc9=inivec, t2tc10=inivec, t2irg=inivec, t2csat1=inivec,
    t2csat2=inivec, t2vent=inivec, t3tc1=inivec, t3tc2=inivec,
    t3tc3=inivec, t3tc4=inivec, t3tc5=inivec, t3tc6=inivec,
    t3tc7=inivec, t3tc8=inivec, t3tc9=inivec, t3tc10=inivec,
    t3vent=inivec)
    for i in 1:size(timenoted, 1)
        pos = findall(x->x>=timenoted[i], outdf[:,1])
        println(size(pos, 1))
        outdf[pos[1],2:end] = df[i, 8:end]
    end
    println("Done")
    return outdf
end
=#
end #module
