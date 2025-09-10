######################################################
###            READ RAW TURBULENCE DATA            ###
###            author: Michi Haugeneder            ###
######################################################
#=
Read the measured rawdata from CSAT/IRGASON/Kaijo and create a file
containing the data for further processing

Workflow to my current (210324) understanding:
- read-in the measured rawdata
- create .csv file for import to DataFrame (only once)

Also containing functions to process the raw TC and vent Tair
data and export it to .csv files
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
#animation = pyimport("matplotlib.animation")
#GridSpec = pyimport("matplotlib.gridspec")
#mpwidgets = pyimport("matplotlib.widgets")

datapath = "/home/michi/Documents/slf/CONTRASTS25/data/turbtower2_slf/PS149_47-1/converted"
importdir = joinpath(@__DIR__, "..", "..")
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .gen

##
# change variables here
#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
#windrawfile = string(datapath, "WindSonic_Duerrboden_2005xx/WindSonic_Dürrboden_combined.dat")
filename = joinpath(datapath, "TOA5_Contra2_sonics_0.dat")
#t1rawfile = joinpath(datapath, "t1_fast_pardenn.dat")
#t2rawfile = joinpath(datapath, "t2tmp", "t2_db_smallac")
#t12_outfile_stam = joinpath(datapath)
#kaijofile = joinpath(datapath, "kaijo", "kaijo.dat")
#kaijo_period_file = joinpath(datapath, "kaijo", "21_kaijo_periods.txt")
#kaijo_outfile_stam = joinpath(datapath, "kaijo", "kaijo")
#windrawfile = joinpath(datapath, "TJK_Dischma_21/WindSonic_Dürrboden2021_combined.dat")
#timestamprawfile = joinpath(datapath, "TJK_Dischma_21/WindSonic_timestamp_Dürrboden2021_combined.DAT")
#outfile_stam = joinpath(datapath, "WindSonic_Duerrboden_2005xx/tjk_sonic_200504-200514")

#blocklength = Minute(30)
##

#########################################################################
###                             RAW DATA I/O                          ###
#########################################################################
println()
println("-----------S-T-A-R-T-------------")
#########################################################################
#uncomment to read in data in filename

#import the raw measurement data and create arrays with data and timestamps
(filetime, filedata) = Main.turb.loadrawgeneric(filename)

#check the column names from the file and enter in the following

#(t2): fill data into DataFrames and output them to file
irg = DataFrame(time = filetime, u = filedata[:,2],
v = filedata[:,3], w = filedata[:,4], T = filedata[:,5], co2 = filedata[:,7],
h2o = filedata[:,8], airpressure = filedata[:,11], diagsonic = filedata[:,6],
diagirg = filedata[:,9], celltemperatur = filedata[:, 10], co2signalstrength = filedata[:,12],h20signalstrength = filedata[:,13])

csat = DataFrame(time = filetime, u = filedata[:,14],
v = filedata[:,15], w = filedata[:,16], T = filedata[:,17], diagsonic = filedata[:,18])

#cut
@info("Adapt the following cutting to your case!")
firstelement = findfirst(x->x>=DateTime(2025,08,26,10,45,00), filetime)
lastelement = findlast(x->x<=DateTime(2025,08,27,15,10,00), filetime)
irg = irg[firstelement:lastelement, :]
csat = csat[firstelement:lastelement, :]

CSV.write(joinpath("/home/michi/Documents/slf/CONTRASTS25/data/processed", "3d_t2_irg.csv"), irg)
CSV.write(joinpath("/home/michi/Documents/slf/CONTRASTS25/data/processed", "3d_t2_csat.csv"), csat)

#########################################################################
#uncomment to read in tower 1 turbulence raw data
#=
#import the raw measurement data and create arrays with data and timestamps
(t1time, t1data) = turb.loadt1raw(t1rawfile)

#(t2): fill data into DataFrames and output them to file
t1irg = DataFrame(time = t1time, u = t1data[:,1],
v = t1data[:,2], w = t1data[:,3], T = t1data[:,4], co2 = t1data[:,5],
h2o = t1data[:,6], airpressure = t1data[:,7], diagsonic = t1data[:,8],
diagirg = t1data[:,9])
CSV.write("/home/haugened/Documents/data/tower/t1irg_new_pardenn.csv", t1irg)
=#
#########################################################################
#=
#uncomment to read in tower 2 turbulence raw data

#import the raw measurement data and create arrays with data and timestamps
(t2time, t2data) = turb.loadt2raw(t2rawfile)

#(t2): fill data into DataFrames and output them to file
t2irg = DataFrame(time = t2time, u = t2data[:,1],
v = t2data[:,2], w = t2data[:,3], T = t2data[:,4], co2 = t2data[:,5],
h2o = t2data[:,6], airpressure = t2data[:,7], diagsonic = t2data[:,8],
diagirg = t2data[:,9])
t2lcsat = DataFrame(time = t2time, u = t2data[:,10], v = t2data[:,11],
w = t2data[:,12], T = t2data[:,13], diagsonic = t2data[:,14])
t2ucsat = DataFrame(time = t2time, u = t2data[:,15], v = t2data[:,16],
w = t2data[:,17], T = t2data[:,18], diagsonic = t2data[:,19])
CSV.write(joinpath(t12_outfile_stam, "t2irg_db_3.csv"), t2irg)
CSV.write(joinpath(t12_outfile_stam, "t2lcsat_db_3.csv"), t2lcsat)
CSV.write(joinpath(t12_outfile_stam, "t2ucsat_db_3.csv"), t2ucsat)
=#
#########################################################################
#uncomment to read in KAIJO raw data
#=
#import the raw measurement data and create arrays with data and timestamps
(kaijotime, kaijodata) = turb.loadkaijoraw(kaijofile)

#fill data into DataFrame and output it to file
kaijodf = DataFrame(time = kaijotime, u = kaijodata[:,1],
v = kaijodata[:,2], w = kaijodata[:,3], T = kaijodata[:,4])
CSV.write(string(kaijo_outfile_stam, ".csv"), kaijodf)
=#
#########################################################################
#uncomment to read in other sonic raw data
#=
#import the raw measurement data and create arrays with data and timestamps
(timestamps, sonicdata) = turb.createtimestamped3Dwind(windrawfile, timestamprawfile, timestep)

#fill data into DataFrame and output it to file
sonicdf = DataFrame(time = timestamps, sonic_u = sonicdata[:,1],
sonic_v = sonicdata[:,2], sonic_w = sonicdata[:,3], sonic_T = sonicdata[:,4])
CSV.write(string(outfile_stam, ".csv"), sonicdf)
=#
#########################################################################
#ventilated air temperature measurement
#=
#tower 2
df2 = CSV.File("/home/haugened/Documents/data/tower/t12_slow_duerrboden.dat"; header=0, skipto=5, ntasks=Threads.nthreads()) |> Tables.matrix
dateformat = DateFormat("yyyy-mm-dd HH:MM:SS")
timeofmeasure2 = DateTime.(df2[:, 1], dateformat)
t2ventdat = replace!(df2[:, 4], Inf => missing)
t2vent = DataFrame(time=timeofmeasure, vent_air_temp=float.(replace!(t2ventdat, NaN => missing)))
turb.saveturbasnetcdf(t2vent, "/home/haugened/Documents/data/tower/vent_air/t2.nc")

#tower 3
df3 = CSV.File("/home/haugened/Documents/data/tower/t3_tairvent_duerrboden.dat"; header=0, skipto=5, ntasks=Threads.nthreads()) |> Tables.matrix
timeofmeasure3 = DateTime.(df3[:, 1], dateformat)
t3ventdat = replace!(df3[:,4], Inf => missing)
t3vent = DataFrame(time=timeofmeasure3, vent_air_temp=float.(replace!(t3ventdat, NaN => missing)))
turb.saveturbasnetcdf(t3vent, "/home/haugened/Documents/data/tower/vent_air/t3.nc")
=#
#########################################################################


#simple plotting stuff
plotkaijodf = replace(datarot[:, 2], missing => NaN)
PyPlot.pygui(true)
fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.set_title("Kaijo data")
ax.set_xlabel("time")
ax.set_ylabel("pot. temperature [deg C]")
ax.grid(true)
ax.plot(datarot[:, 1], plotkaijodf[:])

#########################################################################
println("------------D-O-N-E---------------")
println()
