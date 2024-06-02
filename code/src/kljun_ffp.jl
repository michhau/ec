######################################################
###           FLUX FOOTPRINT CALCULATION           ###
###            author: Michi Haugeneder            ###
######################################################
#=
from: Kljun et al. (2015) `A simple two-dimensional 
parameterisation for Flux Footprint Prediction (FFP)`

adapted the Python script to JULIA
=#
module kljun

using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
LogNorm = pyimport("matplotlib.colors")
import PyPlot

export ffp, plot_ffp

"""
    footprintkljun(zm, z0, umean, h, ol, sigmav, ustar)

Compute the flux footprint according to Kljun et al. (2015).

# Arguments
- `zm::Number`: Measurement height above displacement height (i.e. z-d) [m]
- `z0`: Roughness length [m]; enter `nothing` if not known 
- `umean`: Mean wind speed at zm [m/s]; enter `nothing` if not known 
         Either z0 or umean is required. If both are given,
         z0 is selected to calculate the footprint
- `h::Number`: Boundary layer height [m]
- `ol::Number`: Obukhov length [m]
- `sigmav::Number`: standard deviation of lateral velocity fluctuations [ms-1]
- `ustar::Number`: friction velocity [ms-1]

# optional inputs
- `wind_dir=nothing`: wind direction in degrees (of 360) for rotation of the footprint    
- `rs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]`: Percentage of source area for which to provide contours, must be between 10% and 90%. 
           Can be either a single value (e.g., "80") or a list of values (e.g., "[10, 20, 30]")
           Expressed either in percentages ("80") or as fractions of 1 ("0.8"). 
           Set to "NaN" for no output of percentages
- `nx::Int=1000`: Integer scalar defining the number of grid elements of the scaled footprint.
           Large nx results in higher spatial resolution and higher computing time.
           nx must be >=600.
- `rslayer::Bool=false`: Calculate footprint even if zm within roughness sublayer: set rslayer = true
           Note that this only gives a rough estimate of the footprint as the model is not 
           valid within the roughness sublayer. z0 is needed for estimation of the RS.
- `crop::Bool=false`: Crop output area to size of the 80% footprint or the largest r given if crop=true

# Output
- `x_ci_max`: x location of footprint peak (distance from measurement) [m]
- `x_ci`: x array of crosswind integrated footprint [m]
- `f_ci`: array with footprint function values of crosswind integrated footprint [m-1] 
- `x_2d`: x-grid of 2-dimensional footprint [m], rotated if wind_dir is provided
- `y_2d`: y-grid of 2-dimensional footprint [m], rotated if wind_dir is provided
- `f_2d`: footprint function values of 2-dimensional footprint [m-2]
- `rs`: percentage of footprint as in input, if provided
- `fr`: footprint value at r, if r is provided
- `xr`: x-array for contour line of r, if r is provided
- `yr`: y-array for contour line of r, if r is provided
- `flag_err`: false if no error, true in case of error
"""
function ffp(zm::Number, z0, umean,
    h::Number, ol::Number, sigmav::Number, ustar::Number,
    wind_dir=nothing, rs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    rslayer::Bool=false, nx::Int=1000, crop::Bool=false)

    #########################################
    # Input check
    flag_err = false

    # Check existence of required input pars
    if any(isnothing.([zm, h, ol, sigmav, ustar])) || (isnothing(z0) && isnothing(umean))
        throw(ArgumentError("At least one argument is missing."))
    end

    # Check passed values
    if zm <= 0
        throw(ArgumentError("zm (measurement height) must be larger than zero."))
    end
    #if z0 is not None and umean is None and z0 <= 0.: raise_ffp_exception(3)
    if h <= 10
        throw(ArgumentError("h (BL height) must be larger than 10 m."))
    end
    if zm > h
        throw(ArgumentError("zm (measurement height) must be smaller than h (BL height)"))
    end
    if !isnothing(z0) && isnothing(umean) && zm <= 12.5 * z0
        if rslayer
            @warn("zm (measurement height) should be above roughness sub-layer (12.5*z0).")
        else
            throw(ArgumentError("No valid zm (measurement height above displacement height) passed."))
        end
    end
    if zm / ol <= -15.5
        throw(ArgumentError("zm/ol (measurement height to Obukhov length ratio) must be equal or larger than -15.5."))
    end
    if sigmav <= 0
        throw(ArgumentError("sigmav (standard deviation of crosswind) must be larger than zero."))
    end
    if ustar <= 0.1
        throw(ArgumentError("ustar (friction velocity) must be >=0.1."))
    end
    if !isnothing(wind_dir) && !(0 <= wind_dir <= 360)
        throw(ArgumentError("wind_dir (wind direction) must be >=0 and <=360."))
    end
    if nx < 600
        throw(ArgumentError("nx must be >= 600."))
    end

    # Resolve ambiguity if both z0 and umean are passed (defaults to using z0)
    if !any(isnothing.([z0, umean]))
        @info("Using z0, ignoring umean.")
    end


    #########################################
    # Handle rs
    if !isnothing(rs)
        # Check that rs is a vector, otherwise make it a vector
        if isa(rs, Number)
            if 0.9 < rs <= 1 || 90 < rs <= 100
                rs = [0.9]
            end
        end
        if !isa(rs, Vector)
            throw(ArgumentError("rs should be a vector of numbers."))
        end

        # If rs is passed as percentages, normalize to fractions of one
        if maximum(rs) >= 1
            rs = [x / 100 for x in rs]
        end

        # Eliminate any values beyond 0.9 (90%) and inform user
        if maximum(rs) > 0.9
            @warn("Values > 0.9 for rs are neglected")
            rs = [rs .<= 0.9]
        end

        # Sort levels in ascending order
        rs = sort(rs)
    end


    #########################################
    # Model parameters
    a = 1.4524
    b = -1.9914
    c = 1.4622
    d = 0.1359
    ac = 2.17
    bc = 1.66
    cc = 20.0

    xstar_end = 30
    oln = 5000 #limit to L for neutral scaling
    k = 0.4 #von Karman


    #########################################
    # Scaled X* for crosswind integrated footprint
    xstar_ci_param = collect(d:(xstar_end-d)/(nx+1):xstar_end)[2:end]

    # Crosswind integrated scaled F* 
    fstar_ci_param = a .* (xstar_ci_param .- d) .^ b .* exp.(-c ./ (xstar_ci_param .- d))
    ind_notnan = .!isnan.(fstar_ci_param)
    fstar_ci_param = fstar_ci_param[ind_notnan]
    xstar_ci_param = xstar_ci_param[ind_notnan]

    # Scaled sig_y*
    sigystar_param = ac .* sqrt.(bc .* xstar_ci_param .^ 2 ./ (1 .+ cc .* xstar_ci_param))

    #########################################
    # Real scale x and f_ci
    if !isnothing(z0)
        # Use z0
        if ol <= 0 || ol >= oln
            xx = (1 - 19.0 * zm / ol)^0.25
            psi_f = log((1 + xx^2) / 2.0) + 2.0 * log((1 + xx) / 2.0) - 2.0 * atan(xx) + pi / 2
        elseif ol > 0 && ol < oln
            psi_f = -5.3 * zm / ol
        end
        x = xstar_ci_param * zm / (1.0 - (zm / h)) * (log(zm / z0) - psi_f)
        if log(zm / z0) - psi_f > 0
            x_ci = x
            f_ci = fstar_ci_param ./ zm .* (1 - (zm / h)) / (log(zm / z0) - psi_f)
        else
            throw(ErrorException("log(zm/z0)-psi_f <= 0!"))
            flag_err = true
        end
    else
        # Use umean if z0 not available
        x = xstar_ci_param * zm / (1.0 - zm / h) * (umean / ustar * k)
        if umean / ustar > 0
            x_ci = x
            f_ci = fstar_ci_param / zm * (1.0 - zm / h) / (umean / ustar * k)
        else
            throw(ErrorException("umean/ustar <= 0!"))
            flag_err = true
        end
    end
    #Maximum location of influence (peak location)
    xstarmax = -c / b + d
    if !isnothing(z0)
        x_ci_max = xstarmax * zm / (1.0 - (zm / h)) * (log(zm / z0) - psi_f)
    else
        x_ci_max = xstarmax * zm / (1.0 - (zm / h)) * (umean / ustar * k)
    end
    #Real scale sig_y
    if abs(ol) > oln
        ol = -1E6
    end
    if ol <= 0   #convective
        scale_const = 1E-5 * abs(zm / ol)^(-1) + 0.80
    elseif ol > 0  #stable
        scale_const = 1E-5 * abs(zm / ol)^(-1) + 0.55
    end
    if scale_const > 1
        scale_const = 1.0
    end
    sigy = sigystar_param ./ scale_const .* zm .* sigmav / ustar
    sigy[sigy.<0] .= NaN

    #Real scale f(x,y)
    dx = x_ci[2] - x_ci[1]
    y_pos = collect(0:dx:(length(x_ci)/2)*dx*1.5)
    #f_pos = Array{Float64}(undef, len(f_ci), len(y_pos))
    #fill!(f_pos, NaN)
    f_pos = Array{Float64}(undef, length(f_ci), length(y_pos))
    fill!(f_pos, NaN)
    for ix in 1:length(f_ci)
        f_pos[ix, :] = f_ci[ix] * 1 / (sqrt(2 * pi) * sigy[ix]) * exp.(-y_pos .^ 2 / (2 * sigy[ix]^2))
    end

    #Complete footprint for negative y (symmetrical)
    y_neg = -y_pos[end:-1:1]
    f_neg = f_pos[:, end:-1:1]
    y = vcat(y_neg[1:end-1], y_pos)
    f = transpose(vcat(transpose(f_neg[:, 1:end-1]), transpose(f_pos)))

    #Matrices for output
    x_2d = Array{Float64}(undef, length(x), length(y))
    y_2d = similar(x_2d)
    for colidx in 1:length(y)
        x_2d[:, colidx] = x
    end
    for rowidx in 1:length(x)
        y_2d[rowidx, :] = y
    end
    #np.tile(x[:,None], (1,len(y)))
    #y_2d = np.tile(y.t,(len(x),1))
    f_2d = f

    #########################################
    # Derive footprint ellipsoid incorporating R% of the flux, if requested,
    # starting at peak value.
    dy = dx
    xrs = 0
    yrs = 0
    clevs = 0
    rs_dummy = 0

    if isnothing(rs) && crop
        rs_dummy = [0.8]
    else
        rs_dummy = rs
        clevs = get_contour_levels(f_2d, dx, dy, rs_dummy)
        frs = [item[3] for item in clevs]

        #test-run to get dimensions
        (tmpa, ) = get_contour_vertices(x_2d, y_2d, f_2d, minimum(frs))
        length_tmpa = length(tmpa)
        xrs = Array{Float64}(undef, length_tmpa, length(frs))
        fill!(xrs, NaN)
        yrs = similar(xrs)
        fill!(yrs, NaN)

        for (ix, fr) in enumerate(frs)
            (xr, yr) = get_contour_vertices(x_2d, y_2d, f_2d, fr)
            if size(xr, 1) == 0 #isnan(xr)
                frs[ix] = NaN
            end
            xrs[1:length(xr),ix] = xr
            yrs[1:length(yr),ix] = yr
        end
    end


    #########################################
    # Crop domain and footprint to the largest rs value
    if crop
        if !isnothing(rs)
            dminx = floor(minimum(xrs[:,end]))
            dmaxx = ceil(maximum(xrs[:,end]))
            dminy = floor(minimum(yrs[:,end]))
            dmaxy = ceil(maximum(yrs[:,end]))
        else
            dminx = floor(minimum(xrs))
            dmaxx = ceil(maximum(xrs))
            dminy = floor(minimum(yrs))
            dmaxy = ceil(maximum(yrs))
        end
        #jrange = np.where((y_2d[0] >= dminy) & (y_2d[0] <= dmaxy))[0]
        jrange = findall(x -> x == 1, dminy .<= y_2d[1, :] .<= dmaxy)
        #jrange = np.concatenate(([jrange[0]-1], jrange, [jrange[-1]+1]))
        jrange = vcat([jrange[1] - 1], jrange, [jrange[end] + 1])
        #jrange = jrange[np.where((jrange>=0) & (jrange<=y_2d.shape[0]-1))[0]]
        jrange = jrange[1 .<= jrange .<= size(y_2d, 1)]
        #irange = np.where((x_2d[:,0] >= dminx) & (x_2d[:,0] <= dmaxx))[0]  
        irange = findall(x -> x == 1, dminx .<= x_2d[:, 1] .<= dmaxx)
        #irange = np.concatenate(([irange[0]-1], irange, [irange[-1]+1]))
        irange = vcat([irange[1] - 1], irange, [irange[end] + 1])
        #irange = irange[np.where((irange>=0) & (irange<=x_2d.shape[1]-1))[0]]
        irange = irange[1 .<= irange .<= size(y_2d, 1)]
        #jrange = [[it] for it in jrange]
        x_2d = x_2d[irange, jrange]
        y_2d = y_2d[irange, jrange]
        f_2d = f_2d[irange, jrange]
    end


    #########################################
    #Rotate 3d footprint if requested
    if !isnothing(wind_dir)
        wind_dir = wind_dir * pi / 180
        dist = sqrt.(x_2d .^ 2 + y_2d .^ 2)
        angle = atan.(y_2d, x_2d)
        x_2d = dist .* sin.(wind_dir .- angle)
        y_2d = dist .* cos.(wind_dir .- angle)

        if !isnothing(rs)
            for (ix, r) in enumerate(rs)
                xr_lev = xrs[:,ix] #[x for x in xrs[ix] if !isnan(x)]
                yr_lev = yrs[:,ix] #[x for x in yrs[ix] if !isnan(x)]
                dist = sqrt.(xr_lev .^ 2 + yr_lev .^ 2)
                angle = atan.(yr_lev, xr_lev)
                xr = dist .* sin.(wind_dir .- angle)
                yr = dist .* cos.(wind_dir .- angle)
                xrs[1:length(xr),ix] = xr
                yrs[1:length(yr),ix] = yr    
            end
        end
    end


    return x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d, rs, frs, xrs, yrs, flag_err
end

#=
"Check ffp_climatology input for consistency"
function check_ffp_inputs(ustar::Number, sigmav::Number,
    h::Number, ol::Number, wind_dir::Number,
    zm::Number)::Bool
    if zm <= 0
        throw(ArgumentError("zm (measurement height) must be larger than zero."))
        return false
    end
    #if z0 is not None and umean is None and z0 <= 0.: raise_ffp_exception(3)
    if h <= 10
        throw(ArgumentError("h (BL height) must be larger than 10 m."))
        return false
    end
    if zm > h
        throw(ArgumentError("zm (measurement height) must be smaller than h (BL height)"))
        return false
    end
    if zm / ol <= -15.5
        throw(ArgumentError("zm/ol (measurement height to Obukhov length ratio) must be equal or larger than -15.5."))
        return false
    end
    if sigmav <= 0
        throw(ArgumentError("sigmav (standard deviation of crosswind) must be larger than zero."))
        return false
    end
    if ustar <= 0.1
        throw(ArgumentError("ustar (friction velocity) must be >=0.1."))
        return false
    end
    if any(!(0 <= wind_dir <= 360) || isnan(winddir))
        throw(ArgumentError("wind_dir (wind direction) must be >=0 and <=360."))
        return false
    end
    return true
end

"""
    ffp_climatology(zm, z0, umean, h, ol, sigmav, ustar, wind_dir, domain, dx, dy, nx, ny, rs, rslayer, smooth_data, crop, pulse, verbosity, fig)

Compute the flux footprint climatology according to Kljun et al. (2015).

# Arguments
All vectors need to be of equal length (one value for each time step)
- `zm::Vector`      : Measurement height above displacement height (i.e. z-d) [m]
- `umean::Vector`   : Mean wind speed at zm [m/s]; enter `nothing` if not known 
                      Either z0 or umean is required. If both are given,
                      z0 is selected to calculate the footprint
- `h::Vector`       : Boundary layer height [m]
- `ol::Vector`      : Obukhov length [m]
- `sigmav::Vector`  : standard deviation of lateral velocity fluctuations [ms-1]
- `ustar::Vector`   : friction velocity [ms-1]
- `wind_dir::Vector`: wind direction in degrees (of 360) for rotation of the footprint

# optional inputs
- `domain::Vector=[-1000.0 1000.0 -1000.0 1000.0]`:
            Domain size as an array of [xmin xmax ymin ymax] [m].
            Footprint will be calculated for a measurement at [0 0 zm] m
            Default is smallest area including the r% footprint or [-1000 1000 -1000 1000]m,
            whichever smallest (80% footprint if r not given).
-  `dx::Float64=2.0`, `dy::Float64=2.0`:
            Cell size of domain [m]
            Small dx, dy results in higher spatial resolution and higher computing time
            Default is dx = dy = 2 m. If only dx is given, dx=dy.
- `nx::Int=1000`, `ny::Int=1000`:
            Two integer scalars defining the number of grid elements in x and y
            Large nx/ny result in higher spatial resolution and higher computing time
            Default is nx = ny = 1000. If only nx is given, nx=ny.
            If both dx/dy and nx/ny are given, dx/dy is given priority if the domain is also specified.
- `rs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]`:
            Percentage of source area for which to provide contours, must be between 10% and 90%. 
            Can be either a single value (e.g., "80") or a list of values (e.g., "[10, 20, 30]")
            Expressed either in percentages ("80") or as fractions of 1 ("0.8"). 
            Set to "NaN" for no output of percentages
- `rslayer::Bool=false`: Calculate footprint even if zm within roughness sublayer: set rslayer = true
            Note that this only gives a rough estimate of the footprint as the model is not 
            valid within the roughness sublayer. z0 is needed for estimation of the RS.

- `smooth_data::Bool=true`:
            Apply convolution filter to smooth footprint climatology if smooth_data=true (default)
- `crop::Bool=false`:
            Crop output area to size of the 80% footprint or the largest r given if crop=true
- `pulse::Int=0`:
            Display progress of footprint calculations every pulse-th footprint (e.g., "100")
- `fig::Bool=false`:
            Plot an example figure of the resulting footprint (on the screen): set fig = 1. 
            Default is false (i.e. no figure). 

# Output
#- `x_ci_max`: x location of footprint peak (distance from measurement) [m]
#- `x_ci`: x array of crosswind integrated footprint [m]
#- `f_ci`: array with footprint function values of crosswind integrated footprint [m-1] 
- `x_2d`: x-grid of 2-dimensional footprint [m], rotated if wind_dir is provided
- `y_2d`: y-grid of 2-dimensional footprint [m], rotated if wind_dir is provided
- `fclim_2d`: normalized footprint function values of 2-dimensional footprint [m-2]
- `rs`: percentage of footprint as in input, if provided
- `fr`: footprint value at r, if r is provided
- `xr`: x-array for contour line of r, if r is provided
- `yr`: y-array for contour line of r, if r is provided
- `n`: number of footprints calculated and included in footprint climatology
- `flag_err`: false if no error, true in case of error
"""
function ffp_climatology(zm::Vector{Number}, umean::Vector{Number},
    h::Vector{Number}, ol::Vector{Number}, sigmav::Vector{Number}, ustar::Vector{Number},
    wind_dir::Vector{Number}, domain=[nothing], dx=nothing, dy=nothing,
    nx=nothing, ny=nothing, rs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    rslayer::Bool=false, smooth_data::Bool=true, crop::Bool=false, pulse::Int=-100, fig::Bool=false)

    #########################################
    # Input check
    flag_err = false

    # Check existence of required input pars
    if any(isnothing.([zm, h, ol, sigmav, ustar])) || (isnothing(umean))
        throw(ArgumentError("At least one argument is missing."))
    end

    # Check that all vectors have same length, if not raise an error and exit
    ts_len = length(ustar)
    if any(length(lst) != ts_len for lst in [sigmav, zm, wind_dir, h, ol])
        throw(ArgumentError("Passed data arrays (ustar, wind_dir, sigmav, zm, h, ol) dont all have the same length."))
        exit()
    end

    #########################################
    # Handle rs
    if !isnothing(rs)
        # Check that rs is a vector, otherwise make it a vector
        if isa(rs, Number)
            if 0.9 < rs <= 1 || 90 < rs <= 100
                rs = [0.9]
            end
        end
        if !isa(rs, Vector)
            throw(ArgumentError("rs should be a vector of numbers."))
        end

        # If rs is passed as percentages, normalize to fractions of one
        if maximum(rs) >= 1
            rs = [x / 100 for x in rs]
        end

        # Eliminate any values beyond 0.9 (90%) and inform user
        if maximum(rs) > 0.9
            @warn("Values > 0.9 for rs are neglected")
            rs = [rs .<= 0.9]
        end

        # Sort levels in ascending order
        rs = sort(rs)
    end

    #########################################
    # Define computational domain
    # Check passed values and make some smart assumptions
    if isa(dx, Number) && isnothing(dy)
        dy = dx
    end
    if isa(dy, Number) && isnothing(dx)
        dx = dy
    end
    if !all(isa(item, Number) for item in [dx, dy])
        dx = nothing
        dy = nothing
    end
    if isa(nx, Int) && isnothing(ny)
        ny = nx
    end
    if isa(ny, Int) && isnothing(nx)
        nx = ny
    end
    if !all(isa(item, Int) for item in [nx, ny])
        nx = nothing
        ny = nothing
    end
    if !isa(domain, Vector{Number}) || length(domain) != 4
        domain = [nothing]
    end

    if all(isnothing(item) for item in [dx, nx, domain])
        # If nothing is passed, default domain is a square of 2 Km size centered
        # at the tower with pizel size of 2 meters (hence a 1000x1000 grid)
        domain = [-1000., 1000., -1000., 1000.]
        dx = dy = 2.
        nx = ny = 1000
    elseif !any(isnothing.(domain))
        # If domain is passed, it takes the precendence over anything else
        if !isnothing(dx)
            # If dx/dy is passed, takes precendence over nx/ny
            nx = round(Int,((domain[2]-domain[1]) / dx))
            ny = round(Int,((domain[4]-domain[3]) / dy))
        else
            # If dx/dy is not passed, use nx/ny (set to 1000 if not passed)
            if isnothing(nx)
                nx = 1000
                ny = 1000
            end
            # If dx/dy is not passed, use nx/ny1
            dx = (domain[2]-domain[1]) / float(nx)
            dy = (domain[4]-domain[3]) / float(ny)
        end
    elseif !isnothing(dx) && !isnothing(nx)
        # If domain is not passed but dx/dy and nx/ny are, define domain
        domain = [-nx*dx/2.0, nx*dx/2.0, -ny*dy/2.0, ny*dy/2.0]
    elseif !isnothing(dx)
        # If domain is not passed but dx/dy is, define domain and nx/ny
        domain = [-1000.0, 1000.0, -1000.0, 1000.0]
        nx = round(Int((domain[2]-domain[1]) / dx))
        ny = round(Int((domain[4]-domain[3]) / dy))
    elseif !isnothing(nx)
        # If domain and dx/dy are not passed but nx/ny is, define domain and dx/dy
        domain = [-1000.0, 1000.0, -1000.0, 1000.0]
        dx = (domain[2]-domain[1]) / float(nx)
        dy = (domain[4]-domain[3]) / float(nx)
    end

    # Put domain into more convenient vars
    xmin = domain[1]
    xmax = domain[2]
    ymin = domain[3]
    ymax = domain[4]

    # Define pulse if not passed
    if pulse == -100
        if ts_len <= 20
            pulse = 1
        else
            pulse = round(Int(ts_len / 20))
        end
    end

    #########################################
    # Model parameters
    a = 1.4524
    b = -1.9914
    c = 1.4622
    d = 0.1359
    ac = 2.17
    bc = 1.66
    cc = 20.0

    oln = 5000 #limit to L for neutral scaling
    k = 0.4 #von Karman

    #########################################
    # Define physical domain in cartesian and polar coordinates
    # Cartesian coordinates
    x = LinRange(xmin, xmax, nx + 1)
    y = LinRange(ymin, ymax, ny + 1)
    x_2d = x' .* ones(length(y))
    y_2d = ones(length(x))' .* y 

    # Polar coordinates
    # Set theta such that North is pointing upwards and angles increase clockwise
    rho = sqrt.(x_2d.^2 .+ y_2d.^2)
    theta = atan2.(x_2d, y_2d)

    # initialize raster for footprint climatology
    fclim_2d = zeros(Float64, size(x_2d, 1), size(x_2d, 2))

    #########################################
    # Loop on time series

    # Initialize logic array valids to those 'timestamps' for which all inputs are
    # at least present (but not necessarily phisically plausible)
    valid = zeros(Bool, length(ustar))
    for i in 1:length(valids)
        if !any(isnan.([ustar[i], sigmav[i], h[i], ol[i], wind_dir[i], zm[i]]))
            valid[i] = true
        else
            valid[i] = false
        end
    end

    for (ix, (ustari, sigmavi, hi, oli, wind_diri, zmi, umeani)) in 
        enumerate(zip(ustars, sigmavs, hs, ols, wind_dirs, zms, umeans))

        # Counter
        if mod(ix,pulse) == 0:
            println("Calculating footprint ", ix+1, " of ", ts_len)
        end

        valid[ix] = valid[ix] && check_ffp_inputs(ustari, sigmavi, hi, oli, wind_diri, zmi)

        # If inputs are not valid, skip current footprint
        if valid[ix]
         #   raise_ffp_exception(16, verbosity)
         #else
            #########################################
            # Rotate coordinates into wind direction
            if !isnan[wind_diri]
                rotated_theta = theta - wind_dir * np.pi / 180.
            end

            #########################################
            # Create real scale crosswind integrated footprint and dummy for
            # rotated scaled footprint
            fstar_ci_dummy = np.zeros(x_2d.shape)
            f_ci_dummy = np.zeros(x_2d.shape)
            xstar_ci_dummy = np.zeros(x_2d.shape)
            px = np.ones(x_2d.shape)

            xstar_ci_dummy = (rho * np.cos(rotated_theta) / zm * (1. - (zm / h)) / (umean / ustar * k))
            px = np.where(xstar_ci_dummy > d)
            fstar_ci_dummy[px] = a * (xstar_ci_dummy[px] - d)**b * np.exp(-c / (xstar_ci_dummy[px] - d))
            f_ci_dummy[px] = (fstar_ci_dummy[px] / zm * (1. - (zm / h)) / (umean / ustar * k))
            end
            #########################################
            # Calculate dummy for scaled sig_y* and real scale sig_y
            sigystar_dummy = np.zeros(x_2d.shape)
            sigystar_dummy[px] = (ac * np.sqrt(bc * np.abs(xstar_ci_dummy[px])**2 / (1 +
                                  cc * np.abs(xstar_ci_dummy[px]))))

            if abs(ol) > oln
                ol = -1E6
            end
            if ol <= 0   #convective
                scale_const = 1E-5 * abs(zm / ol)**(-1) + 0.80
            elseif ol > 0  #stable
                scale_const = 1E-5 * abs(zm / ol)**(-1) + 0.55
            end
            if scale_const > 1
                scale_const = 1.0
            end

            sigy_dummy = np.zeros(x_2d.shape)
            sigy_dummy[px] = (sigystar_dummy[px] / scale_const * zm * sigmav / ustar)
            sigy_dummy[sigy_dummy < 0] = np.nan

            #########################################
            # Calculate real scale f(x,y)
            f_2d = np.zeros(x_2d.shape)
            f_2d[px] = (f_ci_dummy[px] / (np.sqrt(2 * np.pi) * sigy_dummy[px]) *
                        np.exp(-(rho[px] * np.sin(rotated_theta[px]))**2 / ( 2. * sigy_dummy[px]**2)))

            #########################################
            # Add to footprint climatology raster
            fclim_2d = fclim_2d + f_2d;
        end
    end

    #########################################
    # Scaled X* for crosswind integrated footprint
    xstar_ci_param = collect(d:(xstar_end-d)/(nx+1):xstar_end)[2:end]

    # Crosswind integrated scaled F* 
    fstar_ci_param = a .* (xstar_ci_param .- d) .^ b .* exp.(-c ./ (xstar_ci_param .- d))
    ind_notnan = .!isnan.(fstar_ci_param)
    fstar_ci_param = fstar_ci_param[ind_notnan]
    xstar_ci_param = xstar_ci_param[ind_notnan]

    # Scaled sig_y*
    sigystar_param = ac .* sqrt.(bc .* xstar_ci_param .^ 2 ./ (1 .+ cc .* xstar_ci_param))

    #########################################
    # Real scale x and f_ci
    if !isnothing(z0)
        # Use z0
        if ol <= 0 || ol >= oln
            xx = (1 - 19.0 * zm / ol)^0.25
            psi_f = log((1 + xx^2) / 2.0) + 2.0 * log((1 + xx) / 2.0) - 2.0 * atan(xx) + pi / 2
        elseif ol > 0 && ol < oln
            psi_f = -5.3 * zm / ol
        end
        x = xstar_ci_param * zm / (1.0 - (zm / h)) * (log(zm / z0) - psi_f)
        if log(zm / z0) - psi_f > 0
            x_ci = x
            f_ci = fstar_ci_param ./ zm .* (1 - (zm / h)) / (log(zm / z0) - psi_f)
        else
            throw(ErrorException("log(zm/z0)-psi_f <= 0!"))
            flag_err = true
        end
    else
        # Use umean if z0 not available
        x = xstar_ci_param * zm / (1.0 - zm / h) * (umean / ustar * k)
        if umean / ustar > 0
            x_ci = x
            f_ci = fstar_ci_param / zm * (1.0 - zm / h) / (umean / ustar * k)
        else
            throw(ErrorException("umean/ustar <= 0!"))
            flag_err = true
        end
    end
    #Maximum location of influence (peak location)
    xstarmax = -c / b + d
    if !isnothing(z0)
        x_ci_max = xstarmax * zm / (1.0 - (zm / h)) * (log(zm / z0) - psi_f)
    else
        x_ci_max = xstarmax * zm / (1.0 - (zm / h)) * (umean / ustar * k)
    end
    #Real scale sig_y
    if abs(ol) > oln
        ol = -1E6
    end
    if ol <= 0   #convective
        scale_const = 1E-5 * abs(zm / ol)^(-1) + 0.80
    elseif ol > 0  #stable
        scale_const = 1E-5 * abs(zm / ol)^(-1) + 0.55
    end
    if scale_const > 1
        scale_const = 1.0
    end
    sigy = sigystar_param ./ scale_const .* zm .* sigmav / ustar
    sigy[sigy.<0] .= NaN

    #Real scale f(x,y)
    dx = x_ci[2] - x_ci[1]
    y_pos = collect(0:dx:(length(x_ci)/2)*dx*1.5)
    #f_pos = Array{Float64}(undef, len(f_ci), len(y_pos))
    #fill!(f_pos, NaN)
    f_pos = Array{Float64}(undef, length(f_ci), length(y_pos))
    fill!(f_pos, NaN)
    for ix in 1:length(f_ci)
        f_pos[ix, :] = f_ci[ix] * 1 / (sqrt(2 * pi) * sigy[ix]) * exp.(-y_pos .^ 2 / (2 * sigy[ix]^2))
    end

    #Complete footprint for negative y (symmetrical)
    y_neg = -y_pos[end:-1:1]
    f_neg = f_pos[:, end:-1:1]
    y = vcat(y_neg[1:end-1], y_pos)
    f = transpose(vcat(transpose(f_neg[:, 1:end-1]), transpose(f_pos)))

    #Matrices for output
    x_2d = Array{Float64}(undef, length(x), length(y))
    y_2d = similar(x_2d)
    for colidx in 1:length(y)
        x_2d[:, colidx] = x
    end
    for rowidx in 1:length(x)
        y_2d[rowidx, :] = y
    end
    #np.tile(x[:,None], (1,len(y)))
    #y_2d = np.tile(y.t,(len(x),1))
    f_2d = f

    #########################################
    # Derive footprint ellipsoid incorporating R% of the flux, if requested,
    # starting at peak value.
    dy = dx
    xrs = 0
    yrs = 0
    clevs = 0
    rs_dummy = 0

    if isnothing(rs) && crop
        rs_dummy = [0.8]
    else
        rs_dummy = rs
        clevs = get_contour_levels(f_2d, dx, dy, rs_dummy)
        frs = [item[3] for item in clevs]

        #test-run to get dimensions
        (tmpa, ) = get_contour_vertices(x_2d, y_2d, f_2d, minimum(frs))
        length_tmpa = length(tmpa)
        xrs = Array{Float64}(undef, length_tmpa, length(frs))
        fill!(xrs, NaN)
        yrs = similar(xrs)
        fill!(yrs, NaN)

        for (ix, fr) in enumerate(frs)
            (xr, yr) = get_contour_vertices(x_2d, y_2d, f_2d, fr)
            if size(xr, 1) == 0 #isnan(xr)
                frs[ix] = NaN
            end
            xrs[1:length(xr),ix] = xr
            yrs[1:length(yr),ix] = yr
        end
    end


    #########################################
    # Crop domain and footprint to the largest rs value
    if crop
        if !isnothing(rs)
            dminx = floor(minimum(xrs[:,end]))
            dmaxx = ceil(maximum(xrs[:,end]))
            dminy = floor(minimum(yrs[:,end]))
            dmaxy = ceil(maximum(yrs[:,end]))
        else
            dminx = floor(minimum(xrs))
            dmaxx = ceil(maximum(xrs))
            dminy = floor(minimum(yrs))
            dmaxy = ceil(maximum(yrs))
        end
        #jrange = np.where((y_2d[0] >= dminy) & (y_2d[0] <= dmaxy))[0]
        jrange = findall(x -> x == 1, dminy .<= y_2d[1, :] .<= dmaxy)
        #jrange = np.concatenate(([jrange[0]-1], jrange, [jrange[-1]+1]))
        jrange = vcat([jrange[1] - 1], jrange, [jrange[end] + 1])
        #jrange = jrange[np.where((jrange>=0) & (jrange<=y_2d.shape[0]-1))[0]]
        jrange = jrange[1 .<= jrange .<= size(y_2d, 1)]
        #irange = np.where((x_2d[:,0] >= dminx) & (x_2d[:,0] <= dmaxx))[0]  
        irange = findall(x -> x == 1, dminx .<= x_2d[:, 1] .<= dmaxx)
        #irange = np.concatenate(([irange[0]-1], irange, [irange[-1]+1]))
        irange = vcat([irange[1] - 1], irange, [irange[end] + 1])
        #irange = irange[np.where((irange>=0) & (irange<=x_2d.shape[1]-1))[0]]
        irange = irange[1 .<= irange .<= size(y_2d, 1)]
        #jrange = [[it] for it in jrange]
        x_2d = x_2d[irange, jrange]
        y_2d = y_2d[irange, jrange]
        f_2d = f_2d[irange, jrange]
    end


    #########################################
    #Rotate 3d footprint if requested
    if !isnothing(wind_dir)
        wind_dir = wind_dir * pi / 180
        dist = sqrt.(x_2d .^ 2 + y_2d .^ 2)
        angle = atan.(y_2d, x_2d)
        x_2d = dist .* sin.(wind_dir .- angle)
        y_2d = dist .* cos.(wind_dir .- angle)

        if !isnothing(rs)
            for (ix, r) in enumerate(rs)
                xr_lev = xrs[:,ix] #[x for x in xrs[ix] if !isnan(x)]
                yr_lev = yrs[:,ix] #[x for x in yrs[ix] if !isnan(x)]
                dist = sqrt.(xr_lev .^ 2 + yr_lev .^ 2)
                angle = atan.(yr_lev, xr_lev)
                xr = dist .* sin.(wind_dir .- angle)
                yr = dist .* cos.(wind_dir .- angle)
                xrs[1:length(xr),ix] = xr
                yrs[1:length(yr),ix] = yr    
            end
        end
    end


    return x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d, rs, frs, xrs, yrs, flag_err
end

=#
"Contour levels of f at percentages of f-integral given by rs"
function get_contour_levels(f, dx, dy, rs=nothing)
    #Check input and resolve to default levels in needed
    if !(isa(rs, Vector) || isa(rs, Number))
        rs = collect(0.1:0.1:0.9)
    end
    if isa(rs, Number)
        rs = [rs]
    end

    #Levels
    pclevs = Vector{Float64}(undef, length(rs))
    fill!(pclevs, NaN)
    ars = similar(pclevs)
    fill!(ars, NaN)
    sf = sort(vec(f))[end:-1:1]
    loopidx = findall(x -> (!isnan(x) && !isinf(x)), sf)

    #csf = msf.cumsum().filled(np.nan)*dx*dy
    csf = similar(sf)
    fill!(csf, NaN)
    csf[loopidx[1]] = sf[loopidx[1]]
    for i in loopidx[2:end]
        csf[i] = csf[i-1] + sf[i]
    end
    csf = csf .* (dx * dy)

    for (ix, r) in enumerate(rs)
        dcsf = abs.(csf .- r)
        idxmin_dcsf = last(findmin(filter(!isnan, dcsf)))
        pclevs[ix] = sf[idxmin_dcsf]
        ars[ix] = csf[idxmin_dcsf]
    end

    return [(round(r, digits=3), ar, pclev) for (r, ar, pclev) in zip(rs, ars, pclevs)]
end

function get_contour_vertices(x, y, f, lev)
    cs = PyPlot.contour(x, y, f, [lev])
    PyPlot.close()
    segs = cs.allsegs[1]#[1]
    xr = segs[:, 1] #[vert[1] for vert in segs]
    yr = segs[:, 2] #[vert[2] for vert in segs]
    #Set contour to NaN if it's found to reach the physical domain
    if minimum(x) >= minimum(segs[:, 1]) || maximum(segs[:, 1]) >= maximum(x) || minimum(y) >= minimum(segs[:, 2]) || maximum(segs[:, 2]) >= maximum(y)
        return [NaN], [NaN]
    end
    return xr, yr   # x,y coords of contour points.	
end

"Plot footprint function and contours if request"
function plot_ffp(x_2d, y_2d, fs, clevs=nothing, show_heatmap::Bool=true, normalize=nothing,
    colormap=nothing, line_width=0.5, iso_labels=nothing)

    # If input is a list of footprints, don't show footprint but only contours,
    # with different colors
    if isa(fs, Vector)
        show_heatmap = false
    else
        fs = [fs]
    end

    if isnothing(colormap)
        colormap = cm.turbo
    end
    # Define colors for each contour set
    cs = [colormap(ix) for ix in collect(0:1/(length(fs)-1):1)]

    # Initialize figure
    (fig, ax) = PyPlot.subplots(figsize=(12, 10))
    # fig.patch.set_facecolor('none')
    # ax.patch.set_facecolor('none')

    if !isnothing(clevs)
        # Temporary patch for pyplot.contour requiring contours to be in ascending orders
        clevs = clevs[end:-1:1]

        # Eliminate contour levels that were set to NaN
        # (e.g. because they extend beyond the defined domain)
        clevs = clevs[!isnan.(clevs)]

        # Plot contour levels of all passed footprints
        # Plot isopleth
        levs = [clev for clev in clevs]
        for (f, c) in zip(fs, cs)
            cc = [c] * length(levs)
            if show_heatmap
                cp = ax.contour(x_2d, y_2d, f, levs, colors="w", linewidths=line_width)
            else
                cp = ax.contour(x_2d, y_2d, f, levs, colors=cc, linewidths=line_width)
            end
            # Isopleth Labels
            if !isnothing(iso_labels)
                pers = [String(round(Int, clev[0] * 100), "%") for clev in clevs]
                fmt = Vector{String}(undef, length(cp.levels))
                for (l, s) in zip(cp.levels, pers)
                    fmt[l] = s
                end
                PyPlot.clabel(cp, cp.levels[:], inline=1, fmt=fmt, fontsize=7)
            end
        end
    end

    # plot footprint heatmap if requested and if only one footprint is passed
    if show_heatmap
        if normalize == "log"
            norm = LogNorm()
        else
            norm = nothing
        end

        pcol = 0
        for f in fs
            pcol = PyPlot.pcolormesh(x_2d, y_2d, f, cmap=colormap, norm=norm)
        end
        PyPlot.xlabel("x [m]")
        PyPlot.ylabel("y [m]")
        PyPlot.gca().set_aspect("equal", "box")

        cbar = fig.colorbar(pcol, shrink=1.0, format="%.3e")
        #cbar.set_label('Flux contribution', color = 'k')
    end

    PyPlot.show()
    return fig, ax
end
end #module
