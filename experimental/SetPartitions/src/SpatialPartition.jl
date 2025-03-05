"""
    SpatialPartition

`SpatialPartition` is a set-partition where the points are 
arranged in multiple levels along the third dimension. 
See Section 2.2 in [CW16](@cite).
"""
struct SpatialPartition <: AbstractPartition
    partition::SetPartition
    levels::Int

    function SpatialPartition(partition::SetPartition, levels::Int)
        @req number_of_upper_points(partition) % levels == 0 "number of upper points not divisible by levels"
        @req number_of_lower_points(partition) % levels == 0 "number of lower points not divisible by levels"
        @req levels > 0 "levels needs to be greater than 0"

        return new(partition, levels)
    end
end

"""
    spatial_partition(partition::SetPartition, m::Int)

Construct a `SpatialPartition` on `m` levels from `partition`.

Note that number of upper and lower points of `partition` have to be a multiple of `m`. 
See Remark 2.4 in [CW16](@cite).
"""
function spatial_partition(partition::SetPartition, m::Int)
    return SpatialPartition(partition, m)
end

"""
    spatial_partition(upper_points::Vector, lower_points::Vector, m::Int)

Construct a `SpatialPartition` on `m` levels from a partition given by 
`upper_points` and `lower_points`.

Note that the length of `upper_points` and `lower_points` have to be a multiple of `m`. 
See Remark 2.4 in [CW16](@cite).
"""
function spatial_partition(upper_points::Vector, lower_points::Vector, m::Int)
    return spatial_partition(set_partition(upper_points, lower_points), m)
end

function hash(p::SpatialPartition, h::UInt)
    return hash(set_partition(p), hash(levels(p), h))
end

function ==(p::SpatialPartition, q::SpatialPartition)
    return set_partition(p) == set_partition(q) && levels(p) == levels(q)
end

function deepcopy_internal(p::SpatialPartition, stackdict::IdDict)
    if haskey(stackdict, p)
        return stackdict[p]
    end
    q = SpatialPartition(deepcopy_internal(p.partition, stackdict), 
                         deepcopy_internal(p.levels, stackdict))
    stackdict[p] = q
    return q
end

"""
    upper_points(p::SpatialPartition)

Return the upper points of `p`.

# Examples
```jldoctest
julia> upper_points(spatial_partition([2, 4], [4, 99], 2))
2-element Vector{Int64}:
 1
 2
```
"""
function upper_points(p::SpatialPartition)
    return upper_points(set_partition(p))
end

"""
    lower_points(p::SpatialPartition)

Return the lower points of `p`.

# Examples
```jldoctest
julia> lower_points(spatial_partition([2, 4], [4, 99], 2))
2-element Vector{Int64}:
 2
 3
```
"""
function lower_points(p::SpatialPartition)
    return lower_points(set_partition(p))
end

"""
    set_partition(p::SpatialPartition)

Return the `SetPartition` part of `p`.

# Examples
```jldoctest
julia> set_partition(spatial_partition([2, 4], [4, 99], 2))
SetPartition([1, 2], [2, 3])
```
"""
function set_partition(p::SpatialPartition)
    return p.partition
end

"""
    levels(p::SpatialPartition)

Return the number of levels of `p`.

# Examples
```jldoctest
julia> levels(spatial_partition([2, 4], [4, 99], 2))
2
```
"""
function levels(p::SpatialPartition)
    return p.levels
end


"""
    tensor_product(p::SpatialPartition, q::SpatialPartition)

Return the tensor product of `p` and `q`.

The tensor product of two spatial partitions is given by their horizontal concatenation.
See also Section 2.3 in [CW16](@cite) and `tensor_product(::SetPartition, ::SetPartition)`.

# Examples
```jldoctest
julia> tensor_product(spatial_partition([1, 2], [2, 1], 2), spatial_partition([1, 1], [1, 2], 2))
SpatialPartition(SetPartition([1, 2, 3, 3], [2, 1, 3, 4]), 2)
```
"""
function tensor_product(p::SpatialPartition, q::SpatialPartition)

    @req levels(p) == levels(q) "p and q have different levels"

    return spatial_partition(tensor_product(set_partition(p), set_partition(q)), levels(p))
end

"""
    involution(p::SpatialPartition)

Return the involution of `p`.

The involution of a spatial partition is obtained by swapping the upper 
and lower points. See also Section 2.3 in [CW16](@cite) and `involution(::SetPartition)`.

# Examples
```jldoctest
julia> involution(spatial_partition([1, 2], [2, 1, 3, 3], 2))
SpatialPartition(SetPartition([1, 2, 3, 3], [2, 1]), 2)
```
"""
function involution(p::SpatialPartition)
    return spatial_partition(involution(set_partition(p)), levels(p))
end

"""
    is_composable(p::SpatialPartition, q::SpatialPartition)

Return whether `p` and `q` are composable, i.e. they have the same 
number of levels and the number of upper points of `p` equal the 
number of lower points of `q`.

# Examples
```jldoctest
julia> is_composable(spatial_partition([1, 2], [2, 1], 2), spatial_partition([1, 2], [1, 1], 2))
true

julia> is_composable(spatial_partition([1, 2], [2, 1], 2), spatial_partition([1, 2], [1, 1], 1))
false
```
"""
function is_composable(p::SpatialPartition, q::SpatialPartition)
    return levels(p) == levels(q) && is_composable(set_partition(p), set_partition(q))
end

"""
    compose_count_loops(p::SpatialPartition, q::SpatialPartition)

Return the composition of `p` and `q` as well as the number of removed loops.

The composition of two spatial partitions is obtained by concatenating them vertically
and removing intermediate loops which are no longer connected to the top or bottom.
See also Section 2.3 in [CW16](@cite) and 
`compose_count_loops(::SetPartition, ::SetPartition)`.

The composition of `p` and `q` is only defined if they have the same 
number of levels and the number of upper points of `p` equal the 
number of lower points of `q`. See also `is_composable(::SpatialPartition)`.

# Examples
```jldoctest
julia> compose_count_loops(spatial_partition([1, 2], [2, 1], 2), spatial_partition([1, 2], [1, 1], 2))
(SpatialPartition(SetPartition([1, 2], [1, 1]), 2), 0)

julia> compose_count_loops(spatial_partition([1, 1], [2, 2], 2), spatial_partition([1, 1], [2, 2], 2))
(SpatialPartition(SetPartition([1, 1], [2, 2]), 2), 1)

julia> compose_count_loops(spatial_partition([1, 2], [2, 1], 2), spatial_partition([1, 2], [1, 1], 1))
ERROR: ArgumentError: p and q have different levels
[...]
```
"""
function compose_count_loops(p::SpatialPartition, q::SpatialPartition)

    @req is_composable(p, q) "p and q have different levels"

    comp_loops = compose_count_loops(set_partition(p), set_partition(q))

    return (spatial_partition(comp_loops[1], levels(p)), comp_loops[2])

end
