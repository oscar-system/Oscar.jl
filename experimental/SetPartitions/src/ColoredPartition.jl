"""
    ColoredPartition

Initialize colored SetPartition object as ColoredPartition.
"""
struct ColoredPartition <: AbstractPartition
    partition::SetPartition
    color_upper_points::Vector
    color_lower_points::Vector
end

function colored_partition(
    partition::SetPartition, 
    color_upper::Vector, 
    color_lower::Vector)

    return ColoredPartition(partition, color_upper, color_lower)

end

function hash(p::ColoredPartition, h::UInt)

    hash(p.partition, hash(p.color_upper_points, hash(p.color_lower_points, h)))
    
end

function ==(p::ColoredPartition, q::ColoredPartition)

    p.partition == q.partition && 
        p.color_upper_points == q.color_upper_points && 
        p.color_lower_points == q.color_lower_points

end

function copy(p::ColoredPartition)
    return ColoredPartition(copy(p.partition), 
        copy(p.color_upper_points), 
        copy(p.color_lower_points))
end


"""
    tensor_product(p::ColoredPartition, q::ColoredPartition)

Apply tensor product of `p` and `q` and return result.
"""
function tensor_product(p::ColoredPartition, q::ColoredPartition)

    ColoredPartition(tensor_product(p.partition, q.partition), 
        vcat(p.color_upper_points, q.color_upper_points), 
        vcat(p.color_lower_points, q.color_lower_points))
end

"""
    involution(p::ColoredPartition)

Apply involution on `p` and return result.
"""
function involution(p::ColoredPartition)

    ColoredPartition(involution(p.partition), p.color_lower_points, p.color_upper_points)

end

"""
    composition_loops(p::ColoredPartition, q::ColoredPartition)

Apply composition between `p` and `q` and return tuple including the result
as well as the number of removed loops.
"""
function composition_loops(p::ColoredPartition, q::ColoredPartition)

    p.color_upper_points != q.color_lower_points ? 
        error("p upper and q lower colors are different in composition") : 

    comp_loops = composition_loops(p.partition, q.partition)
    
    (ColoredPartition(comp_loops[1], q.color_upper_points, p.color_lower_points), 
        comp_loops[2])

end

"""
    rotation(p::ColoredPartition)

Apply rotation on `p` and return result.

# Arguments
- `p`: Input partition
- `lr`: lr whether left (true) or right (false)
- `tb`: tb whether top (true) or bottom (false) rotation
"""
function rotation(p::ColoredPartition, lr::Bool, tb::Bool)

    if tb && isempty(upper_points(p))
        error("Partition has no top part.")
    elseif !tb && isempty(lower_points(p))
        error("Partition has no bottom part.")
    end

    ret = (copy(upper_points(p)), copy(lower_points(p)), 
        copy(p.color_upper_points), copy(p.color_lower_points))

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
    ColoredPartition(SetPartition(ret[1], ret[2]), ret[3], ret[4])
end
