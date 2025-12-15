######################################################
###     SILVEX II DATA PROCESSING TO NETCDF        ###
###           Reads .dat files and combines        ###
######################################################

using Dates, CSV, DataFrames#, Tables

# Include the turb_data module
importdir = joinpath(@__DIR__, "..", "..")
include(joinpath(importdir, "src", "turb_data.jl"))
import .turb

######################################################
###              CONFIGURATION                     ###
######################################################

# Directory containing the .dat files
data_directory = joinpath("/home/haugened/Documents/data/silvex2/SILVEXII_Silvia2_CSAT3B_2m_preprocessed")

# Output directory for NetCDF files
output_dir = joinpath(data_directory, "netcdf")

# NetCDF filename stem (datetime will be appended)
netcdf_file_stem = "silvex2_silvia2_2m"

# Check for missing timestamps across file boundaries
check_cross_file_timestamps = true

# Expected sampling frequency (Hz)
sampling_frequency = 20  # 20 Hz = 0.05 s interval

# Maximum gap (in number of samples) to fill with NaN
fillnan = Second(10)

######################################################
###              MAIN PROCESSING                   ###
######################################################

@info "Starting SILVEX II data processing..."
@info "Data directory: " * data_directory
@info "Output directory: " * output_dir
@info "Cross-file timestamp check: " * string(check_cross_file_timestamps)
@info "Maximum gap to fill: " * string(fillnan)

# Get all .dat files in the directory (excluding *claudeai.dat files)
all_files = filter(x -> endswith(x, ".dat") && !contains(x, "claudeai"), readdir(data_directory, join=true))

if isempty(all_files)
    @error "No .dat files found in directory: " * data_directory
    error("No data files to process")
end

@info "Found " * string(length(all_files)) * " .dat files"

# Sort files by name (which includes timestamp in filename)
sort!(all_files)

# Expected time step
expected_timestep = Millisecond(1000 / sampling_frequency)  # 50 ms for 20 Hz

# Initialize for period tracking
all_periods = Vector{DataFrame}()
period_start_times = Vector{DateTime}()
current_period = DataFrame(
    time = DateTime[],
    u = Float64[],
    v = Float64[],
    w = Float64[],
    T = Float64[],
    diagsonic = Int64[]
)

# Process each file
last_timestamp = nothing

for (file_idx, filepath) in enumerate(all_files)
    filename = basename(filepath)
    @info "Processing file $file_idx/$(length(all_files)): $filename"
    
    # Read the CSV file
    df = CSV.File(filepath; header=1, delim=',') |> DataFrame
    
    # Parse timestamps
    dateformat = DateFormat("yyyy-mm-ddTHH:MM:SS.s")
    timestamps = df.timestamp
    
    # Parse data columns
    u_data = df.Ux
    v_data = df.Uy
    w_data = df.Uz
    T_data = df.Ts
    diagsonic_data = convert.(Int64, replace(df.diag_sonic, NaN => -9999))
    
    # Check cross-file boundary if enabled
    if check_cross_file_timestamps && !isnothing(last_timestamp)
        time_diff = timestamps[1] - last_timestamp
        gap_duration = Dates.value(time_diff)
        fillnan_duration = Dates.value(fillnan)
        
        if time_diff > expected_timestep
            if gap_duration <= fillnan_duration
                # Fill gap with NaN values
                @info "Filling gap at file boundary ($(time_diff)) with NaN values"
                num_missing = Int(Dates.value(time_diff) / Dates.value(expected_timestep)) - 1
                
                for i in 1:num_missing
                    missing_time = last_timestamp + i * expected_timestep
                    push!(current_period, (missing_time, NaN, NaN, NaN, NaN, -9999))
                end
            else
                # Gap too large, start new period
                @warn "Gap at file boundary ($(time_diff)) exceeds fillnan threshold, starting new period"
                
                # Save current period if not empty
                if nrow(current_period) > 0
                    push!(all_periods, copy(current_period))
                    push!(period_start_times, current_period.time[1])
                end
                
                # Start new period
                current_period = DataFrame(
                    time = DateTime[],
                    u = Float64[],
                    v = Float64[],
                    w = Float64[],
                    T = Float64[],
                    diagsonic = Int64[]
                )
            end
        end
    end
    
    # Process timestamps within file
    for i in 1:length(timestamps)
        current_time = timestamps[i]
        
        # Check for gaps within file (except first row)
        if i > 1
            time_diff = current_time - timestamps[i-1]
            gap_duration = Dates.value(time_diff)
            fillnan_duration = Dates.value(fillnan)
            
            if time_diff > expected_timestep
                if gap_duration <= fillnan_duration
                    # Fill gap with NaN values
                    @info "Filling gap in file $filename at row $i ($(time_diff)) with NaN values"
                    num_missing = Int(Dates.value(time_diff) / Dates.value(expected_timestep)) - 1
                    
                    for j in 1:num_missing
                        missing_time = timestamps[i-1] + j * expected_timestep
                        push!(current_period, (missing_time, NaN, NaN, NaN, NaN, -9999))
                    end
                else
                    # Gap too large, start new period
                    @warn "Gap in file $filename at row $i ($(time_diff)) exceeds fillnan threshold, starting new period"
                    
                    # Save current period if not empty
                    if nrow(current_period) > 0
                        push!(all_periods, copy(current_period))
                        push!(period_start_times, current_period.time[1])
                    end
                    
                    # Start new period
                    current_period = DataFrame(
                        time = DateTime[],
                        u = Float64[],
                        v = Float64[],
                        w = Float64[],
                        T = Float64[],
                        diagsonic = Int64[]
                    )
                end
            end
        end
        
        # Add current row to period
        push!(current_period, (current_time, u_data[i], v_data[i], w_data[i], T_data[i], diagsonic_data[i]))
    end
    
    # Update last timestamp for next iteration
    last_timestamp = timestamps[end]
end

# Add final period if not empty
if nrow(current_period) > 0
    push!(all_periods, current_period)
    push!(period_start_times, current_period.time[1])
end

@info "Total records combined: $(sum(nrow(df) for df in all_periods))"
@info "Number of data periods: $(length(all_periods))"

# Create output directory if it doesn't exist
mkpath(output_dir)

# Save each period as a separate NetCDF file
for (idx, period_data) in enumerate(all_periods)
    if nrow(period_data) == 0
        @warn "Period $idx is empty, skipping..."
        continue
    end
    
    # Format start time for filename
    start_time = period_start_times[idx]
    time_str = Dates.format(start_time, "yyyymmdd_HHMMSS")
    
    # Generate output filename
    output_filename = "$(netcdf_file_stem)_$(time_str).nc"
    output_path = joinpath(output_dir, output_filename)
    
    @info "Saving period $idx ($(nrow(period_data)) records) to: $output_filename"
    
    # Save as NetCDF using the turb module function
    turb.saveturbasnetcdf(period_data, output_path)
end

@info "Processing complete!"
@info "Output files saved to: $output_dir"