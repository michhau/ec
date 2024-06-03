######################################################
###        MODULE FOR FOURIER-TRAFO OF IRDATA      ###
###            author: Michi Haugeneder            ###
######################################################
module ft

using LsqFit, FFTW, Statistics, Dates
using PyCall, ProgressMeter
#using PyPlot

export detrend, belltaper, ftforwindow, fourierplus, spectralenergydensity,
findcutofffreq, lowpass, ibelltaper, lowpassallinone,
fouriercutoffmatrix, logavg

######################################################
###                  PREPROCESSING                 ###
######################################################
"After Stull p.307ff. Detrend data with linear fit"
function detrend(invec::Vector)::Vector
    model(x, p) = p[1] * x .+ p[2]
    p0 = [0.0001, 3.0]
    fit = curve_fit(model, collect(1:size(invec, 1)), invec, p0)
    param = fit.param
    return invec - (param[1]*collect(1:size(invec, 1)) .+ param[2])
end

"After Stull p.310 Apply bell tapering to smoothen edges."
function belltaper(invec::Vector)::Vector
    #check equation 8.4.3 on p.310 in Stull
    outvec = similar(invec)
    N = size(invec, 1)
    lim1 = round(Int, 0.1*N)
    lim2 = round(Int, 0.9*N)
    for i in eachindex(invec)
        if (i <= lim1) || (i >= lim2)
            outvec[i] = invec[i] * (sin(5*pi*i/N))^2
        else
            outvec[i] = invec[i]
        end
    end
    return outvec
end

######################################################
###               FOURIER TRANSFORM                ###
######################################################
"Apply pixel-wise FT to window and return average"
function ftforwindow(data::Array, duration::Number, avginterval::AbstractFloat)
    Soff_avg = zeros(Float64, ceil(Int, size(data, 3)/2.0))
    freq_avg = zeros(Float64, ceil(Int, size(data, 3)/2.0))
    freqtimesSoff = zeros(Float64, ceil(Int, size(data, 3)/2.0))
    count = 0
    for idx_col in 1:size(data,2)
        for idx_row in 1:size(data, 1)
            count += 1
            rawvec = filter(!isnan, data[idx_row, idx_col,:])
            vec1 = ft.detrend(rawvec)
            vect = ft.belltaper(vec1)
            (bla, ft1) = ft.fourierplus(vect)
            (Soff1, freq1) = ft.spectralenergydensity(ft1, duration)
            #freqtimesSoff += freq1.*Soff1
            Soff_avg += Soff1
            freq_avg += freq1
        end
    end
    Soff_avg = Soff_avg ./count
    freq_avg = freq_avg ./count
    #freqtimesSoff = freqtimesSoff ./(count)
    freqtimesSoff1 = freq_avg.*Soff_avg
    return ft.logavg(freq_avg, freqtimesSoff1,avginterval)
end

"Fourier-trafo with the necessary folding according to stull"
function fourierplus(invec::Vector)
    #compare JuliaFFT-doc (https://juliamath.github.io/AbstractFFTs.jl/stable/api/#Public-Interface-1)
    #with Stull page 303 (8.4.1b) to see the division by N
    FTvec = fft(invec)/length(invec) #(8.4.1b with n->n+1 (Stull n=0 -> here n=1))

    #folding back the second half according to Stull p. 313
    #resulting in dexcrete spectral intensity (or energy) E_A(n)
    nyfreq = ceil(Int, length(FTvec)/2.0)  #Nyquist frequency
    EAvec = zeros(Float64, nyfreq)
    EAvec[1] = abs.(FTvec[1])^2

    if mod(length(FTvec),2) == 0 #even
        EAvec[2:end-1] = 2 * abs.(FTvec[2:nyfreq-1]).^2
        EAvec[end] = abs.(FTvec[nyfreq])^2
    else
        EAvec[2:end] = 2 * abs.(FTvec[2:nyfreq]).^2
    end
    return FTvec*length(invec), EAvec
end

"Converting to spectral energy density (c.f. (8.6.2b) in Stull)"
function spectralenergydensity(FTin::Vector, duration)
    freq = zeros(Float64, length(FTin))
    for i in 2:length(freq)
        freq[i] = (i-1)/duration
    end
    Sout = similar(FTin)
    Sout[1] = FTin[1]
    Sout[2:end] = FTin[2:end]*duration
    return Sout, freq
end

"For lowpass-filter, below a certain freq 'flim' return idx of next bigger entry in freq-vector"
function findcutofffreq(freqin::Vector, flim)::Integer
    #find index idxflim of limit freq >= flim in freqin
    idxflim = size(freqin, 1)
    for i in 1:size(freqin, 1)
        if freqin[i] >= flim
            idxflim = i
            break
        end
    end
    #println("Desired LPF cutoff-frequency: ", flim, " Hz,")
    #println("actual LPF cutoff-frequency: ", freqin[idxflim], " Hz")
    return idxflim
end

"Apply a low-pass filter and return the Fourier-backtransformed data"
function lowpass(ftrawin::Vector, idxlim::Integer)::Vector
    lpft = ftrawin
    lpft[idxlim:end-idxlim+2] .= 0 #setting desired freqs to 0

    #backtrafo
    lpift = ifft(lpft)
    lpiftreal = zeros(Float64, length(lpift))
    for i in eachindex(lpift)
        if imag(lpift[i]) > exp10(-13)
            println("ERROR: Return of inverse FT is not REAL. Please check.")
            println("imaginary part: ", imag(lpift[i]))
            println("EXIT. Returning 0!")
            return [0.0]
        end
    end
    lpiftreal = real.(lpift)
    return lpiftreal
end

"Inverse bell-tapering after Stull."
function ibelltaper(invec::Vector)::Vector
    #check equation 8.4.3 on p.310 in Stull
    outvec = similar(invec)
    N = size(invec, 1)
    lim1 = round(Int, 0.1*N)
    lim2 = round(Int, 0.9*N)
    for i in eachindex(invec)
        if (i <= lim1) || (i >= lim2)
            outvec[i] = invec[i] / (sin(5*pi*i/N))^2
        else
            outvec[i] = invec[i]
        end
    end
    return outvec
end

######################################################
###                   ALL-IN-ONE                   ###
######################################################
"Apply low-pass filter to input vector. All-in-one-function."
function lowpassallinone(invect::Vector, dur, fco)::Vector
    #apply bell taper and detrending
    vec1 = belltaper(detrend(invect))

    #apply Fourier-Trafo
    (FTraw, FTvec) = fourierplus(vec1)
    (temp,freq) = spectralenergydensity(FTvec, dur)

    #apply inverse Fourier-Trafo
    idxfco = findcutofffreq(freq, fco)

    #back-trafo
    lpift = lowpass(FTraw, idxfco)
    #apply inverse bell tapering
    lpf = ibelltaper(lpift)
    return lpf
end

"Apply 'lowpassallinone' for a matrix"
function fouriercutoffmatrix(inmatrix::Array, dur, fco)
    out = similar(inmatrix)
    IRdetrend = similar(inmatrix)
    @info("Lowpass filtering for array...")
    p = Progress(size(inmatrix, 2))
    update!(p, 0)
    jj=Threads.Atomic{Int}(0)
    l=Threads.SpinLock()
    Threads.@threads for i in 1:size(inmatrix, 2)
      #  println(i, "/", size(inmatrix, 2))
        for j in 1:size(inmatrix, 1)
            if !isnan(maximum(inmatrix[j,i,:]))
                IRdetrend[j,i,:] = detrend(inmatrix[j,i,:])
                out[j,i,:] = lowpassallinone(inmatrix[j,i,:], dur, fco)
            else
                out[j,i,:] .= NaN
                IRdetrend[j,i,:] .= NaN
            end
        end
        Threads.atomic_add!(jj,1)
        Threads.lock(l)
        update!(p,jj[])
        Threads.unlock(l)
    end
    return out, IRdetrend
end

######################################################
###                 PLOTTING STUFF                 ###
######################################################
"Average input vector to one value per 'interval' to plot on a log x-axis without a mess.
'inputvecx' und 'inputvecy' need to be the same length. "
function logavg(inputvecx::Vector, inputvecy::Vector,
    interval)
    #create empty output vectors
    outvecx = Float64[]
    outvecy = Float64[]
    #find limits of vector
    logvecmin = log10(minimum(filter(!isnan, skipmissing(inputvecx[2:end]))))
    logvecmax = log10(maximum(filter(!isnan, skipmissing(inputvecx[2:end]))))
    intervals = collect(logvecmin:interval:logvecmax)
    for iint in 1:(length(intervals)-1)
        idxs = findall(x->intervals[iint] .<= log10(x) .< intervals[iint+1], inputvecx)
        dataxininterval = inputvecx[idxs]
        datayininterval = inputvecy[idxs]
        if length(dataxininterval) > 1
            dataxtmp = maximum(filter(!isnan, skipmissing(dataxininterval)))
            dataytmp = mean(filter(!isnan, skipmissing(datayininterval)))
            push!(outvecx, dataxtmp)
            push!(outvecy, dataytmp)
        elseif length(dataxininterval) == 1
            push!(outvecx, dataxininterval[1])
            push!(outvecy, datayininterval[1])
        end
    end
    return outvecx, outvecy
end


end #module
