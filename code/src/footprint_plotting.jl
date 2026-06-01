"""
    footprint_world_file(image_file)

Return the expected QGIS world-file path next to an orthomosaic image.
"""
function footprint_world_file(image_file::AbstractString)
    stem, _ = splitext(String(image_file))
    world_file = stem * ".jgw"
    isfile(world_file) && return world_file

    uppercase_world_file = stem * ".JGW"
    isfile(uppercase_world_file) && return uppercase_world_file

    return world_file
end

"""
    footprint_meterperpxl(image_file, bgextend_m, bgextend_pxl)

Read meter-per-pixel values from a sibling `.jgw` file when it exists.
World-file x resolution maps to image columns; y resolution maps to rows.
If no world file is available, fall back to the previous TOML extent
calculation.
"""
function footprint_meterperpxl(image_file::Union{AbstractString, Nothing},
        bgextend_m::AbstractVector, bgextend_pxl::AbstractVector)

    if !isnothing(image_file)
        world_file = footprint_world_file(image_file)
        if isfile(world_file)
            values = Float64[]
            for line in eachline(world_file)
                stripped = strip(line)
                isempty(stripped) && continue
                push!(values, parse(Float64, stripped))
            end

            length(values) >= 4 || error("World file $(world_file) must contain at least four numeric rows.")
            return (
                meterperpxl_row=abs(values[4]),
                meterperpxl_col=abs(values[1]),
                source=:jgw,
                world_file=world_file,
            )
        end
    end

    return (
        meterperpxl_row=bgextend_m[1] / bgextend_pxl[1],
        meterperpxl_col=bgextend_m[2] / bgextend_pxl[2],
        source=:extent,
        world_file=nothing,
    )
end

function footprint_background_geometry(fluxloc::AbstractMatrix, bgextend_m::AbstractVector,
        bgextend_pxl::AbstractVector, figorigin::AbstractVector;
        image_file::Union{AbstractString, Nothing}=nothing)
    scale = footprint_meterperpxl(image_file, bgextend_m, bgextend_pxl)
    meterperpxl_row = scale.meterperpxl_row
    meterperpxl_col = scale.meterperpxl_col

    fluxloc_final = Array{Float64}(undef, size(fluxloc, 1), size(fluxloc, 2))
    fluxloc_final[:, 1] = (figorigin[1] .- fluxloc[:, 1]) .* meterperpxl_row
    fluxloc_final[:, 2] = (fluxloc[:, 2] .- figorigin[2]) .* meterperpxl_col

    bgextend_final = (
        -figorigin[2] * meterperpxl_col,
        (bgextend_pxl[2] - 1 - figorigin[2]) * meterperpxl_col,
        -(bgextend_pxl[1] - figorigin[1]) * meterperpxl_row,
        (figorigin[1] - 1) * meterperpxl_row,
    )

    return fluxloc_final, bgextend_final
end
