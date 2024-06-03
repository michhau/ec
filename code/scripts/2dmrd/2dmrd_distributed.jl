######################################################
###           DISTRIBUTED 2D-MRD OF EC-DATA        ###
###            author: Michi Haugeneder            ###
######################################################
#=
Create 2D-MRDs of EC-data on multiple cores. Includes
loading of required data.
=#
using Distributed
using Dates, DataFrames, NCDatasets, Statistics

#add workers
if gethostname() == "Michi-T450s" || "x1carbon5"
    if nprocs() == 1
        addprocs(2)
    end
elseif gethostname() == "LINUX24"
    if nprocs() == 1
        addprocs(10)
    end
else #HYPERION
    using SlurmClusterManager
    addprocs(SlurmManager())
end

@everywhere importdir = joinpath(@__DIR__, "..", "..")
@everywhere datapath = "/home/haugened/mrd/data/"
@everywhere outpath = joinpath(datapath, "2dmrd", "db_tot")
include(joinpath(importdir, "src", "turb_data.jl"))
@everywhere include(joinpath(importdir, "src", "mrd.jl"))
import .turb
@everywhere import .MRD

######################################################
#variables

#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
kaijo_outfile_stam = joinpath(datapath, "kaijo")
tower_outfile_stam = joinpath(datapath, "tower", "preproc")

savefilename = joinpath(outpath, "tjk_uu.nc")

path_to_data_to_use = joinpath(tower_outfile_stam, "tjkdf.nc")
cols_for_mrd = ["u", "u"]
drblocks = false
#select data and measurement period to be evaluated
evalstart = DateTime(2021, 05, 20, 23, 00, 00) #20.05.13:00
evalend = DateTime(2021, 06, 10, 23, 00, 00) #11.06.08:00
#evalend = evalstart + Day(10)

#parameters for 2D-MRD
nrpoints = 200 #100
blocklentime = Minute(30)
shifttime = Minute(1)
######################################################

println()
println("-----------S-T-A-R-T-------------")

#calculating blocklen and shift relative to data length
blocklen = Dates.value(Millisecond(blocklentime))/Dates.value(Millisecond(evalend-evalstart))
shift = Dates.value(Millisecond(shifttime))/Dates.value(Millisecond(blocklentime))

######################################################
###            LOADING & PREPROCESSING             ###
######################################################
evaldf = turb.readturbasnetcdf(path_to_data_to_use, evalstart, evalend)

#show statistics about missing data
turb.printmissstats(evaldf)
#interpolate missing
evaldf = turb.interpolatemissing(evaldf)
#double rotation
@info("Double rotation might cause issues when comparing different data sets!")
if path_to_data_to_use[end-7:end-5] == "tjk"
    turb.drdf!(evaldf, periodwise=false)
else
    turb.drdf!(evaldf)
end
turb.missing2nan!(evaldf)

#create seperate vectors for time and data
data_time = replace(evaldf.time, missing => DateTime(1900, 03, 28, 03, 14, 15))

data = Array{Float64,2}(undef, length(data_time), length(cols_for_mrd))
for (colidx, colstring) in enumerate(cols_for_mrd)
    data[:, colidx] = evaldf[:, colstring]
end

#for calculating mean wind per block: MRD time -> spatial domain
meanwind = sqrt.(evaldf.u .^ 2 .+ evaldf.v .^2 .+ evaldf.w .^2)

evaldf = nothing
data = round.(data, digits=4)

######################################################
###       2D QUASI-CONTINUOS NO-MRD                ###
######################################################
@time (x, mrd_data, time_middle, meanwindout) =
MRD.nomrd2d(data, data_time, meanwind, nrpoints, blocklen, shift; drperblock=drblocks)

MRD.write2dmrddatatonetcdf(x, mrd_data, time_middle,
meanwindout, evalstart, evalend, cols_for_mrd,
nrpoints, blocklen, shift, path_to_data_to_use,
savefilename)