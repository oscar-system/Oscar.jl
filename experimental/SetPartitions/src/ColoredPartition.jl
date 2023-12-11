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

    ColoredPartition(partition, color_upper, color_lower)

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
    ColoredPartition(copy(p.partition), 
        copy(p.color_upper_points), 
        copy(p.color_lower_points))
end


"""
    tensor_product(p::ColoredPartition, q::ColoredPartition)

Return the tensor product of `p` and `q`.
"""
function tensor_product(p::ColoredPartition, q::ColoredPartition)

    ColoredPartition(tensor_product(p.partition, q.partition), 
        vcat(p.color_upper_points, q.color_upper_points), 
        vcat(p.color_lower_points, q.color_lower_points))
end

"""
    involution(p::ColoredPartition)

Return the involution of `p`.
"""
function involution(p::ColoredPartition)

    ColoredPartition(involution(p.partition), p.color_lower_points, p.color_upper_points)

end

"""
    composition_loops(p::ColoredPartition, q::ColoredPartition)

Return the composition of `p` and `q` as well as the number of removed loops.
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

Return the rotation of `p` in the direction given by `lr` and `tb`.

# Arguments
- `p`: input partition
- `lr`: whether left (true) or right (false)
- `tb`: whether top (true) or bottom (false)
"""
function rotation(p::ColoredPartition, lr::Bool, tb::Bool)

    if tb && isempty(upper_points(p))
        error("Partition has no top part")
    elseif !tb && isempty(lower_points(p))
        error("Partition has no bottom part")
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
