"""
    footprint_world_file(image_file)

Return an existing world-file path next to an orthomosaic image.
"""
function footprint_world_file(image_file::AbstractString)
    stem, _ = splitext(String(image_file))
    candidates = [
        stem * ".jgw",
        stem * ".JGW",
        stem * ".pgw",
        stem * ".PGW",
        stem * ".tfw",
        stem * ".TFW",
        stem * ".wld",
        stem * ".WLD",
    ]

    for world_file in candidates
        isfile(world_file) && return world_file
    end

    return first(candidates)
end

"""
    footprint_meterperpxl(image_file, bgextend_m, bgextend_pxl)

Read meter-per-pixel values from a sibling world file when it exists.
World-file x resolution maps to image columns; y resolution maps to rows.
If no world file is available, fall back to the previous TOML extent
calculation.
"""
function footprint_meterperpxl(image_file::Union{AbstractString, Nothing},
        bgextend_m::Union{AbstractVector, Nothing}=nothing,
        bgextend_pxl::Union{AbstractVector, Nothing}=nothing)

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
                source=:world_file,
                world_file=world_file,
            )
        end
    end

    if isnothing(bgextend_m) || isnothing(bgextend_pxl)
        image_label = isnothing(image_file) ? "image" : String(image_file)
        error("No world file found for $(image_label), and bgextend_m/bgextend_pxl fallback extents were not provided.")
    end

    return (
        meterperpxl_row=bgextend_m[1] / bgextend_pxl[1],
        meterperpxl_col=bgextend_m[2] / bgextend_pxl[2],
        source=:extent,
        world_file=nothing,
    )
end

function footprint_image_size(image_size)
    isnothing(image_size) && return nothing
    length(image_size) >= 2 || error("Footprint image size must include row and column dimensions.")
    return Float64(image_size[1]), Float64(image_size[2])
end

function footprint_background_size(scale, bgextend_pxl::Union{AbstractVector, Nothing}, image_size)
    if scale.source == :world_file
        size_from_image = footprint_image_size(image_size)
        isnothing(size_from_image) || return size_from_image
    end

    if !isnothing(bgextend_pxl)
        return Float64(bgextend_pxl[1]), Float64(bgextend_pxl[2])
    end

    size_from_image = footprint_image_size(image_size)
    isnothing(size_from_image) || return size_from_image

    error("Image size is required to build footprint background geometry when no bgextend_pxl fallback is provided.")
end

function footprint_background_geometry(fluxloc::AbstractMatrix,
        bgextend_m::Union{AbstractVector, Nothing},
        bgextend_pxl::Union{AbstractVector, Nothing}, figorigin::AbstractVector;
        image_file::Union{AbstractString, Nothing}=nothing, image_size=nothing)
    scale = footprint_meterperpxl(image_file, bgextend_m, bgextend_pxl)
    meterperpxl_row = scale.meterperpxl_row
    meterperpxl_col = scale.meterperpxl_col
    nrows, ncols = footprint_background_size(scale, bgextend_pxl, image_size)

    fluxloc_final = Array{Float64}(undef, size(fluxloc, 1), size(fluxloc, 2))
    fluxloc_final[:, 1] = (figorigin[1] .- fluxloc[:, 1]) .* meterperpxl_row
    fluxloc_final[:, 2] = (fluxloc[:, 2] .- figorigin[2]) .* meterperpxl_col

    bgextend_final = (
        -figorigin[2] * meterperpxl_col,
        (ncols - 1 - figorigin[2]) * meterperpxl_col,
        -(nrows - figorigin[1]) * meterperpxl_row,
        (figorigin[1] - 1) * meterperpxl_row,
    )

    return fluxloc_final, bgextend_final
end
