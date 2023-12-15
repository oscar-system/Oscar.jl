"""
    ColoredPartition

A colored partition is a set-partition where the points 
can have different colors. See Section 1.1 in [TW18](@cite).

It is represented by a set-partition `partition` and colors `color_upper_points` 
and `color_lower_points` for the upper and lower points.
"""
struct ColoredPartition <: AbstractPartition
    partition::SetPartition
    color_upper_points::Vector{Int}
    color_lower_points::Vector{Int}
end

function colored_partition(
    partition::SetPartition, 
    color_upper::Vector{Int}, 
    color_lower::Vector{Int})

    return ColoredPartition(partition, color_upper, color_lower)
end

function hash(p::ColoredPartition, h::UInt)

    return hash(p.partition, hash(p.color_upper_points, hash(p.color_lower_points, h)))
    
end

function ==(p::ColoredPartition, q::ColoredPartition)

    return p.partition == q.partition && 
        p.color_upper_points == q.color_upper_points && 
        p.color_lower_points == q.color_lower_points

end

function deepcopy_internal(p::ColoredPartition, stackdict::IdDict)
    if haskey(stackdict, p)
        return stackdict[p]
    end
    q = ColoredPartition(deepcopy_internal(p.partition, stackdict), 
                         deepcopy_internal(p.color_upper_points, stackdict), 
                         deepcopy_internal(p.color_lower_points, stackdict))
    stackdict[p] = q
    return q
end

function upper_points(p::ColoredPartition)
    return upper_points(p.partition)
end

function lower_points(p::ColoredPartition)
    return lower_points(p.partition)
end

function upper_colors(p::ColoredPartition)
    return p.upper_colors
end

function lower_colors(p::ColoredPartition)
    return p.lower_colors
end


"""
    tensor_product(p::ColoredPartition, q::ColoredPartition)

Return the tensor product of `p` and `q`.
"""
function tensor_product(p::ColoredPartition, q::ColoredPartition)

    return ColoredPartition(tensor_product(p.partition, q.partition), 
        vcat(p.color_upper_points, q.color_upper_points), 
        vcat(p.color_lower_points, q.color_lower_points))
end

"""
    involution(p::ColoredPartition)

Return the involution of `p`.
"""
function involution(p::ColoredPartition)

    return ColoredPartition(involution(p.partition), 
        p.color_lower_points, p.color_upper_points)

end

function is_composable(p::ColoredPartition, q::ColoredPartition)
    return p.color_upper_points == q.color_lower_points && 
        is_composable(p.partition, q.partition)
end

"""
    composition_loops(p::ColoredPartition, q::ColoredPartition)

Return the composition of `p` and `q` as well as the number of removed loops.
"""
function composition_loops(p::ColoredPartition, q::ColoredPartition)

    @req p.color_upper_points == q.color_lower_points "p upper and q 
        lower colors are different in composition"

    comp_loops = composition_loops(p.partition, q.partition)
    
    return (ColoredPartition(comp_loops[1], q.color_upper_points, p.color_lower_points), 
        comp_loops[2])

end

"""
    rotation(p::ColoredPartition)

Return the rotation of `p` in the direction given by `lr` and `tb`.

# Arguments
- `p`: input partition
- `lr`: whether left (true) or right (false)
- `tb`: whether top (true) or bottom (false)
"""
function rotation(p::ColoredPartition, lr::Bool, tb::Bool)

    if tb
        @req !isempty(upper_points(p)) "SetPartition has no top part"
    elseif !tb
        @req !isempty(lower_points(p)) "SetPartition has no bottom part"
    end

    ret = (deepcopy(upper_points(p)), deepcopy(lower_points(p)), 
        deepcopy(p.color_upper_points), deepcopy(p.color_lower_points))

    if lr
        if tb
            a = ret[1][1]
            splice!(ret[1], 1)
            pushfirst!(ret[2], a)

            a = ret[3][1]
            splice!(ret[3], 1)
            pushfirst!(ret[4], Int(!Bool(a)))
        else
            a = ret[2][1]
            splice!(ret[2], 1)
            pushfirst!(ret[1], a)

            a = ret[4][1]
            splice!(ret[4], 1)
            pushfirst!(ret[3], Int(!Bool(a)))
        end
    else
        if tb
            a = ret[1][end]
            pop!(ret[1])
            push!(ret[2], a)

            a = ret[3][end]
            pop!(ret[3])
            push!(ret[4], Int(!Bool(a)))
        else
            a = ret[2][end]
            pop!(ret[2])
            push!(ret[1], a)

            a = ret[4][end]
            pop!(ret[4])
            push!(ret[3], Int(!Bool(a)))
        end
    end
    return ColoredPartition(SetPartition(ret[1], ret[2]), ret[3], ret[4])
end
