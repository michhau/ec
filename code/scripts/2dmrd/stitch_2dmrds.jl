######################################################
###       STITCH SINGLE 2D-MRDS (.NC-FILES)        ###
###          author: Michi Haugeneder              ###
######################################################
#=
Stitch (in time domain) single .nc-2dmrd files to one big one
=#
using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings
import PyPlot, CSV

importdir = joinpath(@__DIR__, "..", "..")
datapath = "/home/haugened/Documents/data/"
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
import .turb
import .MRD

######################################################
#variables

file1 = joinpath(datapath, "2dmrd", "short_periods", "20051300_to_25051300_T2UCSAT.nc")
file2 = joinpath(datapath, "2dmrd", "short_periods", "25051400_to_28050800_T2UCSAT.nc")
file3 = joinpath(datapath, "2dmrd", "short_periods", "28050900_to_31051000_T2UCSAT.nc")
file4 = joinpath(datapath, "2dmrd", "short_periods", "31051100_to_08060900_T2UCSAT.nc")
file5 = joinpath(datapath, "2dmrd", "short_periods",      "08061000_to_end_T2UCSAT.nc")

outfile = joinpath(datapath, "2dmrd", "db_tot", "T2UCSAT_total.nc")

######################################################
println("-----------S-T-A-R-T-------------")
@info("Careful. Script needs adaption for different cases!")

######################################################
###                LOAD AND PLOT                   ###
######################################################
(file1_time_middle, file1_x, file1_mrd_data, file1_meanwind,
    file1_evalstart, file1_evalend, file1_cols_for_mrd,
    file1_nrpoints, file1_blocklen, file1_shift, file1_sourcefile) =
    MRD.read2dmrdfromnetcdf(file1)

(file2_time_middle, file2_x, file2_mrd_data, file2_meanwind,
    file2_evalstart, file2_evalend, file2_cols_for_mrd,
    file2_nrpoints, file2_blocklen, file2_shift, file2_sourcefile) =
    MRD.read2dmrdfromnetcdf(file2)

(file3_time_middle, file3_x, file3_mrd_data, file3_meanwind,
    file3_evalstart, file3_evalend, file3_cols_for_mrd,
    file3_nrpoints, file3_blocklen, file3_shift, file3_sourcefile) =
    MRD.read2dmrdfromnetcdf(file3)

(file4_time_middle, file4_x, file4_mrd_data, file4_meanwind,
    file4_evalstart, file4_evalend, file4_cols_for_mrd,
    file4_nrpoints, file4_blocklen, file4_shift, file4_sourcefile) =
    MRD.read2dmrdfromnetcdf(file4)

(file5_time_middle, file5_x, file5_mrd_data, file5_meanwind,
    file5_evalstart, file5_evalend, file5_cols_for_mrd,
    file5_nrpoints, file5_blocklen, file5_shift, file5_sourcefile) =
    MRD.read2dmrdfromnetcdf(file5)

file1inst = MRD.instrumentnamefromsourcefile(file1_sourcefile)
file2inst = MRD.instrumentnamefromsourcefile(file2_sourcefile)
file3inst = MRD.instrumentnamefromsourcefile(file3_sourcefile)
file4inst = MRD.instrumentnamefromsourcefile(file4_sourcefile)
file5inst = MRD.instrumentnamefromsourcefile(file5_sourcefile)

if !(file1inst == file2inst == file3inst == file4inst == file5inst)
    @warn("Not all instruments the same. Be sure you really want to stitch!")
end

if !(file1_cols_for_mrd == file2_cols_for_mrd == file3_cols_for_mrd == file4_cols_for_mrd == file5_cols_for_mrd)
    @warn("Not all columns the same. Different decompositions! Makes no sense to stitch!")
    exit()
end

overallstart = file1_time_middle[1]
overallend = file5_time_middle[end]

evalstart_out = file1_evalstart
evalend_out = file5_evalend
time_middle_out = collect(overallstart:Minute(1):overallend)
sourcefile_out = file3_sourcefile

maxsizex = maximum([size(file1_x, 1), size(file2_x, 1), size(file3_x, 1), size(file4_x, 1), size(file5_x, 1)])

time_middle_out = Vector{DateTime}(undef, 0)
x_out = Array{Millisecond}(undef, maxsizex, 0)
mrd_data_out = Array{Float64}(undef, maxsizex, 0)
meanwind_out = Vector{Float64}(undef, 0)

files = [file1, file2, file3, file4, file5]
lastentry = DateTime(1995, 03, 28, 03, 14, 15)

for (ix, istr) in enumerate(files)
    (time_middle_tmp, x_tmp, mrd_data_tmp, meanwind_tmp,
        evalstart_tmp, evalend_tmp, cols_for_mrd_tmp,
        nrpoints_tmp, blocklen_tmp, shift_tmp, sourcefile_tmp) =
        MRD.read2dmrdfromnetcdf(istr)
    if ix == 1
        time_middle_out = time_middle_tmp
        x_out = x_tmp
        mrd_data_out = mrd_data_tmp
        meanwind_out = meanwind_tmp
    else
        #fill missing time
        time_fill = collect(lastentry+Second(30):Second(30):time_middle_tmp[1]-Second(30))
        entriestofill = length(time_fill)
        time_middle_out = vcat(time_middle_out, time_fill)
        meanwind_out = vcat(meanwind_out, fill(NaN, entriestofill))
        x_out = cat(x_out, fill(Millisecond(999999), maxsizex, entriestofill); dims=2)
        mrd_data_out = cat(mrd_data_out, fill(NaN, maxsizex, entriestofill); dims=2)

        #fill data
        time_middle_out = vcat(time_middle_out, time_middle_tmp)
        meanwind_out = vcat(meanwind_out, meanwind_tmp)
        x_out = cat(x_out, x_tmp; dims=2)
        mrd_data_out = cat(mrd_data_out, mrd_data_tmp; dims=2)
    end
    lastentry = time_middle_tmp[end]
end

MRD.write2dmrddatatonetcdf(x_out, mrd_data_out, time_middle_out, meanwind_out, evalstart_out,
evalend_out, file1_cols_for_mrd, maxsizex, file1_blocklen, file1_shift, sourcefile_out, outfile)