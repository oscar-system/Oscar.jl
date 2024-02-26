"""
    PartitionedPermutation

The type of partitioned permutations. Fieldnames are
- p::Perm{Int} - a permutation
- V::SetPartition - a partition
- check::Bool = true
If the permutation has length `n`, then the partition must have `n` upper points and 0 lower points. 
Further, if `W` is the partition given by the cycles of `p`, then `W` must be dominated by `V` in the 
sense that every block of `W` is contained in one block of `V`. There is one inner constructer of PartitionedPermutation:
- PartitionedPermutation(_p::Perm{Int}, _V::Vector{Int}) constructs the partitioned permutation where the partition is given by the vector _V.
If the optional flag `check` is set to `false`, then the constructor skips the validation of the requirements mentioned above.

# Examples
```jldoctest
julia> PartitionedPermutation(Perm([2, 1, 3]), [1, 1, 2])
PartitionedPermutation((1,2), SetPartition([1, 1, 2], Int64[]))
```
"""
struct PartitionedPermutation
    p::Perm{Int}
    V::SetPartition
    check::Bool

    function PartitionedPermutation(_p::Perm{Int}, _V::Vector{Int}, _check::Bool=true)
        __V = SetPartition(_V, Int[])
        if _check
            @req parent(_p).n == length(_V) "permutation and partition must have the same length"
            @req cycle_partition(_p) <= __V "permutation must be dominated by partition"
        end
        new(_p, __V)
    end
end

function partitioned_permutation(p::Perm{Int}, V::Vector{Int}, check::Bool=true)
    return PartitionedPermutation(p, V, check)
end

function ==(pp_1::PartitionedPermutation, pp_2::PartitionedPermutation)
    return pp_1.p == pp_2.p && pp_1.V == pp_2.V
end

function hash(pp::PartitionedPermutation, h::UInt)
    return hash(pp.p, hash(pp.V, h))
end

function deepcopy_internal(pp::PartitionedPermutation, stackdict::IdDict)
    if haskey(stackdict, pp)
        return stackdict[pp]
    end
    q = PartitionedPermutation(deepcopy_internal(pp.p, stackdict), 
                     deepcopy_internal(pp.V, stackdict))
    stackdict[pp] = q
    return q
end

"""
    length(pp::PartitionedPermutation)

Return the length of a partitioned permutation, i.e. the size of the underlying set.

# Examples
```jldoctest
julia> length(PartitionedPermutation(Perm([2, 1]), [1, 1]))
2
```
"""
function length(pp::PartitionedPermutation)
    return parent(pp.p).n
end

"""
    length2(pp::PartitionedPermutation)

Return the adjusted length of a partitioned permutation as described in [CITE] as `|(V, pi)|`
for a partition `V` and a permutation `pi`.

# Examples
```jldoctest
julia> length2(PartitionedPermutation(Perm([2, 1]), [1, 1]))
1
```
"""
function length2(pp::PartitionedPermutation)
    p = pp.p
    V = pp.V
    return parent(p).n - (2*number_of_blocks(V) - length(cycles(p)))
end
