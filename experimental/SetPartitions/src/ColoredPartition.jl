import Base.hash
import Base.==
import Base.copy

"""
ColoredPartition

Initialize colored Partition object

# Arguments
- `partition`: Partition object which is generalized with coloring `color_upper_points` and `color_lower_points`
- `color_upper_points`: color (in [0, 1]) of upper points
- `color_lower_points`: color (in [0, 1]) of lower points
"""
struct ColoredPartition <: AbstractPartition
    partition::Partition
    color_upper_points::Array{Int64, 1}
    color_lower_points::Array{Int64, 1}
end

function colored_partition(partition::Partition, color_upper::Array{Int64, 1}, color_lower::Array{Int64, 1})

    return ColoredPartition(partition, color_upper, color_lower)

end

function hash(p::ColoredPartition, h::UInt)

    hash(p.partition, hash(p.color_upper_points, hash(p.color_lower_points, h)))
    
end

function ==(p::ColoredPartition, q::ColoredPartition)

    p.partition == q.partition && p.color_upper_points == q.color_upper_points && p.color_lower_points == q.color_lower_points

end

function copy(p::ColoredPartition)
    return ColoredPartition(copy(p.partition), copy(p.color_upper_points), copy(p.color_lower_points))
end


"""
tensor_product(p::ColoredPartition, q::ColoredPartition)

This function applies on p tensor product with q (in O(n)).

# Arguments
- `p`: Input colored partition
- `q`: Second input colored partition

# Returns
- `p` tensor product `q`
"""
function tensor_product(p::ColoredPartition, q::ColoredPartition)

    ColoredPartition(tensor_product(p.partition, q.partition), vcat(p.color_upper_points, q.color_upper_points), vcat(p.color_lower_points, q.color_lower_points))
end

"""
involution(p::ColoredPartition)

This function applies an involution on `p` (in O(n) because normal_form, else O(1)).

# Arguments
- `p`: Input colored partition

# Returns
- involution of `p`
"""
function involution(p::ColoredPartition)

    ColoredPartition(involution(p.partition), p.color_lower_points, p.color_upper_points)

end

"""
composition_loops(p::ColoredPartition, q::ColoredPartition)

This function applies composition between p and q (in O(nlogn)).

# Arguments
- `p`: Input partition
- `q`: Second input partition

# Returns
- [`p` composition `q`, number of loops]
"""
function composition_loops(p::ColoredPartition, q::ColoredPartition)

    p.color_upper_points != q.color_lower_points ? error("p upper and q lower colors are different in composition") : 

    comp_loops = composition_loops(p.partition, q.partition)
    
    [ColoredPartition(comp_loops[1], q.color_upper_points, p.color_lower_points), comp_loops[2]]

end

"""
rotation(p::ColoredPartition)

This function applies a rotation on `p` (in O(n) because normal_form, else O(1)). 
Throws error if rotation not possible.

# Arguments
- `p`: Input colored partition
- `lr`: lr whether left (true) or right (false)
- `tb`: tb whether top (true) or bottom (false) rotation

# Returns
- rotation of `p`
"""
function rotation(p::ColoredPartition, lr::Bool, tb::Bool)

    if tb && isempty(get_upper_points(p))
        error("Partition has no top part.")
    elseif !tb && isempty(get_lower_points(p))
        error("Partition has no bottom part.")
    end

    ret::Array = [copy(get_upper_points(p)), copy(get_lower_points(p)), copy(p.color_upper_points), copy(p.color_lower_points)]

        if lr
            if tb
                a::Int = ret[1][1]
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
    ColoredPartition(Partition(ret[1], ret[2]), ret[3], ret[4])
end


