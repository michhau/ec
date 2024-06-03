######################################################
###     COMBINING SINGLE 2D-MRDS (.NC-FILES)       ###
###          author: Michi Haugeneder              ###
######################################################
#=
Combine multiple single 2D-MRDs for
friction velocity u_star
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

file1 = joinpath(datapath, "2dmrd", "db_tot", "raw", "t2ucsat_uw.nc")
file2 = joinpath(datapath, "2dmrd", "db_tot", "raw", "t2ucsat_vw.nc")
outfile = joinpath(datapath, "2dmrd", "db_tot",      "t2ucsat_ustar.nc")

######################################################
println("-----------S-T-A-R-T-------------")
@info("Careful. Script needs adaption for different cases!")

######################################################
###                    COMBINE                     ###
######################################################
(file1_time_middle, file1_x, file1_mrd_data, file1_meanwind,
    file1_evalstart, file1_evalend, file1_cols_for_mrd,
    file1_nrpoints, file1_blocklen, file1_shift, file1_sourcefile) =
    MRD.read2dmrdfromnetcdf(file1)

(file2_time_middle, file2_x, file2_mrd_data, file2_meanwind,
    file2_evalstart, file2_evalend, file2_cols_for_mrd,
    file2_nrpoints, file2_blocklen, file2_shift, file2_sourcefile) =
    MRD.read2dmrdfromnetcdf(file2)

#combining step
mrd_data_out = (file1_mrd_data .^2 .+ file2_mrd_data .^2) .^(1/4)

MRD.write2dmrddatatonetcdf(file1_x, mrd_data_out, file1_time_middle, file1_meanwind, file1_evalstart,
file1_evalend, ["ustar"], size(file1_x, 1), file1_blocklen, file1_shift, file1_sourcefile, outfile)