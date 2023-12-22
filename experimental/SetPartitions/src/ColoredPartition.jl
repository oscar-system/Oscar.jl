"""
    ColoredPartition

`ColoredPartition` represents a set-partition where the points 
can have different colors. See Section 1.1 in [TW18](@cite).
"""
struct ColoredPartition <: AbstractPartition
    partition::SetPartition
    color_upper_points::Vector{Int}
    color_lower_points::Vector{Int}

    function ColoredPartition(partition::SetPartition, 
                              color_upper_points::Vector{Int},
                              color_lower_points::Vector{Int})

        @req all(x -> x in (0, 1), color_upper_points) "upper coloring has to be binary in {0, 1}"
        @req num_upper_points(partition) == length(color_upper_points) "upper coloring format does not match upper points format"
        @req all(x -> x in (0, 1), color_lower_points) "lower coloring has to be binary in {0, 1}"
        @req num_lower_points(partition) == length(color_lower_points) "lower coloring format does not match and lower points format"

        return new(partition, color_upper_points, color_lower_points)
    end
end

"""
    colored_partition(partition::SetPartition, upper_colors::Vector, lower_colors::Vector)

Construct a `ColoredPartition` from `partition` with colors given by 
`upper_colors` and `lower_colors`.

# Examples
```jldoctest
julia> colored_partition(set_partition([2, 4], [4, 99]), [1, 0], [0, 1])
ColoredPartition(SetPartition([1, 2], [2, 3]), [1, 0], [0, 1])
```
"""
function colored_partition(partition::SetPartition, 
                           upper_colors::Vector, lower_colors::Vector)

    return ColoredPartition(partition, 
                            convert(Vector{Int}, upper_colors), 
                            convert(Vector{Int}, lower_colors))
end

"""
    colored_partition(upper_points::Vector, lower_points::Vector,
                      upper_colors::Vector, lower_colors::Vector))

Construct a `ColoredPartition` with underlying partition given by `upper_points` and 
`lower_points` and with colors `upper_colors` and `lower_colors`.

# Examples
```jldoctest
julia> colored_partition([2, 4], [4, 99], [1, 0], [0, 1])
ColoredPartition(SetPartition([1, 2], [2, 3]), [1, 0], [0, 1])
```
"""
function colored_partition(upper_points::Vector,lower_points::Vector,
                           upper_colors::Vector, lower_colors::Vector)

    return colored_partition(set_partition(upper_points, lower_points), 
                             upper_colors, lower_colors)
end

function hash(p::ColoredPartition, h::UInt)
    return hash(set_partition(p), hash(upper_colors(p), hash(lower_colors(p), h)))
end

function ==(p::ColoredPartition, q::ColoredPartition)
    return set_partition(p) == set_partition(q) && 
           upper_colors(p)  == upper_colors(q)  && 
           lower_colors(p)  == lower_colors(q)
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

"""
    set_partition(p::ColoredPartition)

Return the `SetPartition` part of `p`.

# Examples
```jldoctest
julia> set_partition(colored_partition([2, 4], [4, 99], [1, 0], [0, 1]))
SetPartition([1, 2], [2, 3])
```
"""
function set_partition(p::ColoredPartition)
    return p.partition
end

"""
    upper_points(p::ColoredPartition)

Return the upper points of `p`.

# Examples
```jldoctest
julia> upper_points(colored_partition([2, 4], [4, 99], [1, 0], [0, 1]))
[1, 2]
```
"""
function upper_points(p::ColoredPartition)
    return upper_points(set_partition(p))
end

"""
    lower_points(p::ColoredPartition)

Return the lower points of `p`.

# Examples
```jldoctest
julia> lower_points(colored_partition([2, 4], [4, 99], [1, 0], [0, 1]))
[2, 3]
```
"""
function lower_points(p::ColoredPartition)
    return lower_points(set_partition(p))
end

"""
    upper_colors(p::ColoredPartition)

Return the upper colors of `p`.

# Examples
```jldoctest
julia> upper_colors(colored_partition([2, 4], [4, 99], [1, 0], [0, 1]))
[1, 0]
```
"""
function upper_colors(p::ColoredPartition)
    return p.color_upper_points
end

"""
    lower_colors(p::ColoredPartition)

Return the lower colors of `p`.

# Examples
```jldoctest
julia> lower_colors(colored_partition([2, 4], [4, 99], [1, 0], [0, 1]))
[0, 1]
```
"""
function lower_colors(p::ColoredPartition)
    return p.color_lower_points
end

"""
    tensor_product(p::ColoredPartition, q::ColoredPartition)

Return the tensor product of `p` and `q`.

The tensor product of two colored partitions is given by their 
horizontal concatenation. See also Section 1.2 in [TW18](@cite) and 
`tensor_product(::SetPartition, ::SetPartition)`.

# Examples
```jldoctest
julia> tensor_product(colored_partition([1, 2], [2, 1], [1, 0], [1, 1]), colored_partition([1, 1], [1], [0, 1], [0]))
ColoredPartition(SetPartition([1, 2, 3, 3], [2, 1, 3]), [1, 0, 0, 1], [1, 1, 0])
```
"""
function tensor_product(p::ColoredPartition, q::ColoredPartition)
    return colored_partition(tensor_product(set_partition(p), set_partition(q)), 
                             vcat(upper_colors(p), upper_colors(q)), 
                             vcat(lower_colors(p), lower_colors(q)))
end

"""
    involution(p::ColoredPartition)

Return the involution of `p`.

The involution of a colored partition is obtained by swapping the upper 
and lower points and colors. See also Section 1.2 in [TW18](@cite) and 
`involution(::SetPartition)`.

# Examples
```jldoctest
julia> involution(colored_partition([1, 2, 3], [2, 1], [1, 1, 0], [0, 1]))
ColoredPartition(SetPartition([1, 2], [2, 1, 3]), [0, 1], [1, 1, 0])
```
"""
function involution(p::ColoredPartition)
    return colored_partition(involution(set_partition(p)), 
                             lower_colors(p), upper_colors(p))
end

"""
    is_composable(p::ColoredPartition, q::ColoredPartition)

Return whether `p` and `q` are composable, i.e. the upper colors of 
`p` equal the lower colors of `q`.

# Examples
```jldoctest
julia> is_composable(colored_partition([1, 2], [2, 1], [1, 1], [0, 1]), colored_partition([1, 2], [1, 1], [0, 1], [1, 1]))
true

julia> is_composable(colored_partition([1, 2], [2, 1], [0, 1], [1, 1]), colored_partition([1, 2], [1, 1], [0, 1], [1, 1]))
false
```
"""
function is_composable(p::ColoredPartition, q::ColoredPartition)
    return upper_colors(p) == lower_colors(q) && 
           is_composable(set_partition(p), set_partition(q))
end

"""
    compose_count_loops(p::ColoredPartition, q::ColoredPartition)

Return the composition of `p` and `q` as well as the number of removed loops.

The composition of two colored partitions is obtained by concatenating them vertically
and removing intermediate loops which are no longer connected to the top or bottom.
See also Section 1.2 in [TW18](@cite) and 
`compose_count_loops(::SetPartition, ::SetPartition)`.

The composition of `p` and `q` is only defined if the upper colors of 
`p` equals the lower colors of `q`. See also `is_composable(::ColoredPartition)`.

# Examples
```jldoctest
julia> compose_count_loops(colored_partition([1, 2], [2, 1], [1, 1], [0, 1]), colored_partition([1], [1, 1], [0], [1, 1]))
(ColoredPartition(SetPartition([1], [1, 1]), [0], [0, 1]), 0)

julia> compose_count_loops(colored_partition([1, 1], [2], [1, 1], [0]), colored_partition([1], [2, 2], [1], [1, 1]))
(ColoredPartition(SetPartition([1], [2]), [1], [0]), 1)

julia> compose_count_loops(colored_partition([1], [1, 2], [0], [1, 0]), colored_partition([1], [2, 2], [0], [0, 1]))
ERROR: ArgumentError: upper and lower colors are different
```
"""
function compose_count_loops(p::ColoredPartition, q::ColoredPartition)

    @req upper_colors(p) == lower_colors(q) "upper and lower colors are different"

    comp_loops = compose_count_loops(set_partition(p), set_partition(q))
    
    return (colored_partition(comp_loops[1], upper_colors(q), lower_colors(p)), 
            comp_loops[2])

end

"""
    rotate(p::ColoredPartition, lr::Bool, tb::Bool)

Rotate `p` in the direction given by `lr` and `tb`.

Rotating a colored partition moves the left- or right-most point of the upper points 
to the lower points or vice verca and inverts its color. 
See also Section 1.2 in [TW18](@cite) and `rotate(::SetPartition, ::Bool, ::Bool)`.

# Arguments
- `p`: input partition
- `lr`: rotating at the left (true) or at the right (false)
- `tb`: rotating from top to bottom (true) or from bottom to top (false)

# Examples
```jldoctest
julia> rotate(colored_partition([1, 2, 3], [2, 1], [1, 0, 1], [1, 1]), true, true)
ColoredPartition(SetPartition([1, 2], [3, 1, 3]), [0, 1], [0, 1, 1])

julia> rotate(colored_partition([1, 2, 3], [2, 1], [1, 0, 1], [1, 1]), true, false)
ColoredPartition(SetPartition([1, 2, 1, 3], [2]), [0, 1, 0, 1], [1])

julia> rotate(colored_partition([1, 2, 3], [2, 1], [1, 0, 1], [1, 1]), false, true)
ColoredPartition(SetPartition([1, 2], [2, 1, 3]), [1, 0], [1, 1, 0])

julia> rotate(colored_partition([1, 2, 3], [2, 1], [1, 0, 1], [1, 1]), false, false)
ColoredPartition(SetPartition([1, 2, 3, 1], [2]), [1, 0, 1, 0], [1])
```
"""
function rotate(p::ColoredPartition, lr::Bool, tb::Bool)

    if tb
        @req !isempty(upper_points(p)) "partition has no top part"
    elseif !tb
        @req !isempty(lower_points(p)) "partition has no bottom part"
    end

    ret = deepcopy((upper_points(p), lower_points(p), 
                    upper_colors(p), lower_colors(p)))

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
    return colored_partition(ret[1], ret[2], ret[3], ret[4])
end
