export
    isomorphism, size_groundset,
    bases, nonbases, circuits, hyperplanes, flats, cyclic_flats, closure, 
    rank, nullity, 
    fundamental_circuit, fundamental_cocircuit,
    spanning_sets, independent_sets,
    cobases, cocircuits, cohyperplanes, corank,
    is_clutter,
    is_regular, is_binary, is_ternary,
    n_connected_components, connected_components, is_connected,
    loops, coloops, is_loopless, is_coloopless, is_simple, direct_sum_components,
    connectivity_function, is_vertical_k_separation, is_k_separation,
    vertical_connectivity, girth, tutte_connectivity,
    tutte_polynomial, characteristic_polynomial, charpoly, reduced_characteristic_polynomial,
    revlex_basis_encoding

################################################################################
##  Properties and basic functions
################################################################################


@doc Markdown.doc"""
    isomorphism(M::matroid, gs::AbstractVector)

Return a matroid isomorphic to `M` on the groundset `gs`.

# Example
```jldoctest
julia> isomorphism(fano_matroid(), [2,3,4,1,5,6,7])
Matroid of rank 3 on 7 elements
```
"""
function isomorphism(M::Matroid, gs::GroundsetType)
    if length(M.groundset)!=length(gs)
        error("Sets of different size")
    end
    gs2num = create_gs2num(gs)
    return Matroid(M.pm_matroid,gs,gs2num)
end


@doc Markdown.doc"""
    size_groundset(M::matroid)

Return the size of the ground set of the matroid `M`.

# Example
```jldoctest
julia> size_groundset(fano_matroid())
7
```
"""
size_groundset(M::Matroid) = length(M.groundset)

@doc Markdown.doc"""
    bases(M::matroid)

Return the list of bases of the matroid `M`.

# Example
```jldoctest
julia> bases(uniform_matroid(2, 3))
3-element Vector{Vector{Int64}}:
 [1, 2]
 [1, 3]
 [2, 3]
```
"""
bases(M::Matroid) = [[M.groundset[i+1] for i in sort(collect(C))] for C in Vector{Set{Int}}(M.pm_matroid.BASES)]

@doc Markdown.doc"""
    nonbases(M::matroid)

Return the list of nonbases of the matroid `M`.

# Example
```jldoctest
julia> nonbases(fano_matroid())
7-element Vector{Vector{Int64}}:
 [1, 6, 7]
 [2, 5, 7]
 [3, 4, 7]
 [3, 5, 6]
 [2, 4, 6]
 [1, 4, 5]
 [1, 2, 3]
```
"""
function nonbases(M::Matroid)
    Bases = bases(M)
    return filter(set -> set ∉ Bases, Oscar.Hecke.subsets(Vector(1:length(M.groundset)),rank(M)))
end 


@doc Markdown.doc"""
    circuits(M::matroid)

Return the list of circuits of the matroid `M`.

# Example
```jldoctest
julia> circuits(uniform_matroid(2, 4))
4-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 2, 4]
 [1, 3, 4]
 [2, 3, 4]
```
"""
circuits(M::Matroid) = [[M.groundset[i+1] for i in sort(collect(C))] for C in Vector{Set{Int}}(M.pm_matroid.CIRCUITS)]

@doc Markdown.doc"""
    hyperplanes(M::matroid)

Return the list of hyperplanes of the matroid `M`.

# Example
```jldoctest
julia> hyperplanes(fano_matroid())
7-element Vector{Vector{Int64}}:
 [3, 5, 6]
 [3, 4, 7]
 [2, 5, 7]
 [2, 4, 6]
 [1, 6, 7]
 [1, 4, 5]
 [1, 2, 3]
```
"""
hyperplanes(M::Matroid) = [[M.groundset[i+1] for i in sort(collect(C))] for C in Vector{Set{Int}}(M.pm_matroid.MATROID_HYPERPLANES)]

@doc Markdown.doc"""
    flats(M::matroid, r::Union{Int,Nothing})

Return the list of flats of the matroid `M`.
By default all flats are returned.
One may specify a rank `r` as the second parameter in which case only the flats of rank `r` are returned.

# Example
```jldoctest
julia> M = fano_matroid()
Matroid of rank 3 on 7 elements

julia> flats(M)
16-element Vector{Vector}:
 Any[]
 [1]
 [2]
 [3]
 [4]
 [5]
 [6]
 [7]
 [1, 2, 3]
 [1, 4, 5]
 [1, 6, 7]
 [2, 4, 6]
 [2, 5, 7]
 [3, 5, 6]
 [3, 4, 7]
 [1, 2, 3, 4, 5, 6, 7]

julia> flats(M, 2)
7-element Vector{Vector}:
 [1, 2, 3]
 [1, 4, 5]
 [1, 6, 7]
 [2, 4, 6]
 [2, 5, 7]
 [3, 5, 6]
 [3, 4, 7]
```
"""
function flats(M::Matroid, r::Union{Int,Nothing}=nothing)
    num_flats = M.pm_matroid.LATTICE_OF_FLATS.N_NODES
    pm_flats = M.pm_matroid.LATTICE_OF_FLATS.FACES
    jl_flats = [Vector{Int}(Polymake.to_one_based_indexing(Polymake._get_entry(pm_flats, i))) for i in 0:(num_flats-1)]
    if M.pm_matroid.LATTICE_OF_FLATS.TOP_NODE==0
        jl_flats = reverse(jl_flats)
    end
    matroid_flats = [[M.groundset[i] for i in flat] for flat in jl_flats]
    if r!=nothing
        if r<0 || r>rank(M)
            error("The specified rank needs to be between 0 and the rank of the matroid.")
        end
        matroid_flats = filter(flat -> rank(M,flat)==r, matroid_flats)
    end
    return matroid_flats
end

@doc Markdown.doc"""
    cyclic_flats(M::matroid, r::Union{Int,Nothing})

Return the list of cyclic flats of the matroid `M`.
These are the flats that are the union of cycles.
See Section 2.1 in Oxl11 (@cite).

By default all cylic flats are returned.
One may specify a rank `r` as the second parameter.
In this case only the cylic flats of this rank are returned.

# Example
```jldoctest
julia> M = fano_matroid()
Matroid of rank 3 on 7 elements

julia> cyclic_flats(M)
9-element Vector{Vector}:
 Any[]
 [1, 2, 3]
 [1, 4, 5]
 [1, 6, 7]
 [2, 4, 6]
 [2, 5, 7]
 [3, 5, 6]
 [3, 4, 7]
 [1, 2, 3, 4, 5, 6, 7]

julia> cyclic_flats(M, 2)
7-element Vector{Vector}:
 [1, 2, 3]
 [1, 4, 5]
 [1, 6, 7]
 [2, 4, 6]
 [2, 5, 7]
 [3, 5, 6]
 [3, 4, 7]
```
"""
function cyclic_flats(M::Matroid, r::Union{Int,Nothing}=nothing)
    num_flats = M.pm_matroid.LATTICE_OF_CYCLIC_FLATS.N_NODES
    pm_flats = M.pm_matroid.LATTICE_OF_CYCLIC_FLATS.FACES
    jl_flats = [Vector{Int}(Polymake.to_one_based_indexing(Polymake._get_entry(pm_flats, i))) for i in 0:(num_flats-1)]
    if M.pm_matroid.LATTICE_OF_FLATS.TOP_NODE==0
        jl_flats = reverse(jl_flats)
    end
    matroid_flats = [[M.groundset[i] for i in flat] for flat in jl_flats]
    if r!=nothing
        if r<0 || r>rank(M)
            error("The specified rank needs to be between 0 and the rank of the matroid.")
        end
        matroid_flats = filter(flat -> rank(M,flat)==r, matroid_flats)
    end
    return matroid_flats
end

@doc Markdown.doc"""
    closure(M::matroid, set::Vector)

Return the closure of `set` in the matroid `M`.

# Example
```jldoctest
julia> closure(fano_matroid(), [1,2])
3-element Vector{Int64}:
 1
 2
 3
```
"""
function closure(M::Matroid, set::GroundsetType)
    cl = M.groundset
    for flat in flats(M)
        if issubset(set,flat) && issubset(flat,cl)
            cl = flat
        end
    end
    return cl
end

@doc Markdown.doc"""
    rank(M::matroid)

Return the rank of the matroid `M`.

# Example
```jldoctest
julia> rank(fano_matroid())
3
```
"""
rank(M::Matroid) = M.pm_matroid.RANK::Int

@doc Markdown.doc"""
    rank(M::matroid, set::Vector)

Return the rank of `set` in the matroid `M`.

# Example
```jldoctest
julia> M = fano_matroid();

julia> rank(M, [1,2,3])
2
```
"""
function rank(M::Matroid, set::GroundsetType)
    if length(set)==0
        return 0
    else
        return Polymake.matroid.rank( M.pm_matroid, Set([M.gs2num[i]-1 for i in set]) )
    end
end

@doc Markdown.doc"""
    nullity(M::matroid, set::Vector)

Return the nullity of `set` in the matroid `M`.
This is defined to be |set| - rk(set).

# Example
```jldoctest
julia> M = fano_matroid();

julia> nullity(M, [1,2,3])
1
```
"""
nullity(M::Matroid, set::GroundsetType) = length(set)-Polymake.matroid.rank( M.pm_matroid, Set([M.gs2num[i]-1 for i in set]) )


@doc Markdown.doc"""
    fundamental_circuit(M::matroid, basis::Vector, elem::ElementType)

Return the unique circuit contained in the union of `basis` and `elem` of the matroid `M`.
See Section 1.2 of Oxl11 (@cite).
Note that `elem` needs to be in the complement of the `basis` in this case.


# Example
```jldoctest
julia> M = fano_matroid();

julia> fundamental_circuit(M, [1,2,4], 7)
4-element Vector{Int64}:
 1
 2
 4
 7

julia> fundamental_circuit(M, [1,2,4], 3)
3-element Vector{Int64}:
 1
 2
 3
```
"""
function fundamental_circuit(M::Matroid, basis::GroundsetType, elem::ElementType)
    if !(basis in bases(M))
        error("The set is not a basis of M")
    end
    if !(elem in M.groundset)
        error("The element is not in the groundset")
    end
    if elem in basis
        error("The element is in the basis, but has to be in the complement")
    end
    return circuits(restriction(M,[elem; basis]))[1]
end

@doc Markdown.doc"""
    fundamental_cocircuit(M::matroid, cobasis::Vector, elem::ElementType)

Return the unique circuit of the dual matroid of `M` in the union of the complement of `basis` and `elem`.
See Section 2.1 of Oxl11 (@cite).
Note that `elem` needs to be an element of the `basis` in this case.


# Example
```jldoctest
julia> fundamental_cocircuit(fano_matroid(), [1,2,4], 4)
4-element Vector{Int64}:
 4
 5
 6
 7
```
"""
fundamental_cocircuit(M::Matroid, basis::GroundsetType, elem::ElementType) = fundamental_circuit(dual_matroid(M), setdiff(M.groundset, basis), elem)

@doc Markdown.doc"""
    independet_sets(M::matroid)

Return the list of independent sets of the matroid `M`.
These are all subsets of the bases.

# Example
```jldoctest
julia> independent_sets(uniform_matroid(2, 3))
7-element Vector{Vector{Int64}}:
 []
 [1]
 [2]
 [3]
 [1, 3]
 [2, 3]
 [1, 2]
```
"""
function independent_sets(M::Matroid)
    pm_bases = Vector{Set{Int}}(M.pm_matroid.BASES)
    n = length(M.groundset)
    gs = M.groundset
    sets = Vector{Vector{Int64}}()
    push!(sets,[])
    for k in 1:rank(M)
        for set in Oscar.Hecke.subsets(Vector(0:n-1),k)
            for B in pm_bases
                if issubset(set,B)
                    push!( sets, [gs[i+1] for i in set])
                    break
                end
            end
        end
    end
    return sets
end

@doc Markdown.doc"""
    spanning_sets(M::matroid)

Return the list of spanning sets of the matroid `M`.
These are all sets containing a basis.

# Example
```jldoctest
julia> spanning_sets(uniform_matroid(2, 3))
4-element Vector{Vector{Int64}}:
 [1, 2]
 [1, 3]
 [2, 3]
 [1, 2, 3]
```
"""
function spanning_sets(M::Matroid)
    # To avoid code duplication we use that spanning sets are the complements of independent sets in the dual matroid.
    coindependent_sets = independent_sets(dual_matroid(M))
    span_sets = [filter(k -> !(k in set), 1:size_groundset(M)) for set in coindependent_sets]
    return reverse(span_sets)
end

@doc Markdown.doc"""
    cobases(M::matroid)

Return the bases of the dual matroid of `M`.
See Section 2 in Oxl11 (@cite).

# Example
```jldoctest
julia> cobases(uniform_matroid(2, 3))
3-element Vector{Vector{Int64}}:
 [3]
 [2]
 [1]
```
"""
cobases(M::Matroid) = bases(dual_matroid(M))

@doc Markdown.doc"""
    cocircuits(M::matroid)

Return the circuits of the dual matroid of `M`.
See Section 2 in Oxl11 (@cite).

# Example
```jldoctest
julia> cocircuits(uniform_matroid(2, 5))
5-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
 [1, 2, 3, 5]
 [1, 2, 4, 5]
 [1, 3, 4, 5]
 [2, 3, 4, 5]
```
"""
cocircuits(M::Matroid) = circuits(dual_matroid(M))

@doc Markdown.doc"""
    cohyperplanes(M::matroid)

Return the hyperplanes of the dual matroid of `M`.
See Section 2 in Oxl11 (@cite).

# Example
```jldoctest
julia> cohyperplanes(fano_matroid())
14-element Vector{Vector{Int64}}:
 [4, 5, 6, 7]
 [2, 3, 6, 7]
 [2, 3, 4, 5]
 [1, 3, 5, 7]
 [1, 3, 4, 6]
 [1, 2, 5, 6]
 [1, 2, 4, 7]
 [3, 5, 6]
 [3, 4, 7]
 [2, 5, 7]
 [2, 4, 6]
 [1, 6, 7]
 [1, 4, 5]
 [1, 2, 3]
```
"""
cohyperplanes(M::Matroid) = hyperplanes(dual_matroid(M))

@doc Markdown.doc"""
    corank(M::matroid, set::Vector)

Return the rank of `set` in the dual matroid of `M`.

# Example
```jldoctest
julia> corank(fano_matroid(), [1,2,3])
4
```
"""
corank(M::Matroid, set::GroundsetType) = length(set)-rank(M, set) + rank(M, setdiff(M.groundset, set))

@doc Markdown.doc"""
    is_clutter(sets)

Checks if the collection of subsets `sets` is a clutter.
A collection of subsets is a clutter if none of the sets is a proper subset of another.
See Section 2.1 in Oxl11 (@cite).

# Example
```jldoctest
julia> is_clutter([[1,2], [1,2,3]])
false

julia> is_clutter(circuits(fano_matroid()))
true
```
"""
function is_clutter(sets::AbstractVector{T}) where T <: GroundsetType
    for A in sets
        for B in sets
            if A!=B && (issubset(A,B) || issubset(B,A))
                return false
            end
        end
    end
    return true
end

@doc Markdown.doc"""
    is_regular(M::Matroid)

Checks if the matroid `M` is regular, that is representable over every field.
See Section 6.6 in Oxl11 (@cite).

# Example
```jldoctest
julia> is_regular(uniform_matroid(2, 3))
true

julia> is_regular(fano_matroid())
false
```
"""
is_regular(M::Matroid) = M.pm_matroid.REGULAR

@doc Markdown.doc"""
    is_binary(M::Matroid)

Checks if the matroid `M` is binary, that is representable over the finite field F2.
See Section 6.5 in Oxl11 (@cite).

# Example
```jldoctest
julia> is_binary(uniform_matroid(2, 4))
false

julia> is_binary(fano_matroid())
true
```
"""
is_binary(M::Matroid) = M.pm_matroid.BINARY

@doc Markdown.doc"""
    is_ternary(M::Matroid)

Checks if the matroid `M` is ternary, that is representable over the finite field F3.
See Section 4.1 in Oxl11 (@cite).

# Example
```jldoctest
julia> is_ternary(uniform_matroid(2, 4))
true

julia> is_ternary(fano_matroid())
false
```
"""
is_ternary(M::Matroid) = M.pm_matroid.TERNARY

@doc Markdown.doc"""
    n_connected_components(M::Matroid)

Return the number of connected components of `M`.
See Section 4.1 in Oxl11 (@cite).

# Example
```jldoctest
julia> n_connected_components(fano_matroid())
1

julia> n_connected_components(uniform_matroid(3, 3))
3
```
"""
n_connected_components(M::Matroid) = length(M.pm_matroid.CONNECTED_COMPONENTS)

@doc Markdown.doc"""
    connected_components(M::Matroid)

Return the connected components of `M`. The function returns a partition of the ground set where each part corresponds to one connected component. 
See Section 4.1 in Oxl11 (@cite).

# Example
```jldoctest
julia> connected_components(fano_matroid())
1-element Vector{Vector{Int64}}:
 [1, 2, 3, 4, 5, 6, 7]

julia> connected_components(uniform_matroid(3, 3))
3-element Vector{Vector{Int64}}:
 [1]
 [2]
 [3]
```
"""
connected_components(M::Matroid) = [[M.groundset[i+1] for i in comp] for comp in M.pm_matroid.CONNECTED_COMPONENTS]

@doc Markdown.doc"""
    is_connected(M::Matroid)

Check if the matroid `M` is connected, that is has one connected component
See Section 4.1 in Oxl11 (@cite).

# Example
```jldoctest
julia> is_connected(fano_matroid())
true

julia> is_connected(uniform_matroid(3, 3))
false
```
"""
is_connected(M::Matroid) = M.pm_matroid.CONNECTED

@doc Markdown.doc"""
    loops(M::Matroid)

Return the loops of `M`. A loop is an element of the ground set that is not contained in any basis.

# Example
```jldoctest
julia> loops(matroid_from_bases([[1,2]], 4))
2-element Vector{Int64}:
 3
 4

julia> loops(fano_matroid())
Any[]
```
"""
loops(M::Matroid) = [M.groundset[i+1] for i in M.pm_matroid.LOOPS]

@doc Markdown.doc"""
    coloops(M::Matroid)

Return the coloops of `M`. A coloop is an element of the ground set that is contained in every basis.

# Example
```jldoctest
julia> coloops(matroid_from_bases([[1,2]], 4))
2-element Vector{Int64}:
 1
 2

julia> coloops(fano_matroid())
Any[]
```
"""
coloops(M::Matroid) = [M.groundset[i+1] for i in M.pm_matroid.DUAL.LOOPS]

@doc Markdown.doc"""
    is_loopless(M::Matroid)

Check if `M` has a loop. Return `true` if `M` does not have a loop.
See also `loops`.

# Example
```jldoctest
julia> is_loopless(matroid_from_bases([[1,2]], 4))
false

julia> is_loopless(fano_matroid())
true
```
"""
is_loopless(M::Matroid) = length(M.pm_matroid.LOOPS)==0 ? true : false

@doc Markdown.doc"""
    is_coloopless(M::Matroid)

Check if `M` has a coloop. Return `true` if `M` does not have a coloop.
See also `coloops`.

# Example
```jldoctest
julia> is_coloopless(matroid_from_bases([[1,2]], 4))
false

julia> is_coloopless(fano_matroid())
true
```
"""
is_coloopless(M::Matroid) = length(M.pm_matroid.DUAL.LOOPS)==0 ? true : false

@doc Markdown.doc"""
    is_simple(M::Matroid)

Check if `M` has is simple. A matroid is simple if it doesn't have loops and doesn't have parallel elements.
Return `true` if `M` is simple.
See also `loops`.

# Example
```jldoctest
julia> is_simple(matroid_from_bases([[1,2]], 4))
false

julia> is_simple(fano_matroid())
true
```
"""
is_simple(M::Matroid) = M.pm_matroid.SIMPLE

@doc Markdown.doc"""
    direct_sum_components(M::Matroid)

Return the connected components of `M` as a list of matroids.
See Section 4.1 in Oxl11 (@cite).

# Example
```jldoctest
julia> direct_sum_components(fano_matroid())
1-element Vector{Matroid}:
 Matroid of rank 3 on 7 elements

julia> direct_sum_components(uniform_matroid(3, 3))
3-element Vector{Matroid}:
 Matroid of rank 1 on 1 elements
 Matroid of rank 1 on 1 elements
 Matroid of rank 1 on 1 elements
```
"""
function direct_sum_components(M::Matroid)
    res = Vector{Matroid}()
    for set in connected_components(M)
        res = [res; restriction(M,set)]
    end
    return res
end

@doc Markdown.doc"""
    connectivity_function(M::Matroid, set::Vector)

Return the value of the connectivity function of `set` in the matroid `M`.
See Section 8.1 in Oxl11 (@cite).

# Example
```jldoctest
julia> connectivity_function(fano_matroid(), [1,2,4])
3

```
"""
function connectivity_function(M::Matroid, set::GroundsetType)
    return rank(M,set) + rank(M,setdiff(Set(M.groundset), set)) - rank(M)
end

@doc Markdown.doc"""
    is_vertical_k_separation(M::Matroid, k::Int, set::Vector)

Check if `set` together with its complement defines a `k` separation in `M`
See Section 8.6 in Oxl11 (@cite).

# Example
```jldoctest
julia> is_vertical_k_separation(fano_matroid(), 2, [1,2,4])
false

```
"""
function is_vertical_k_separation(M::Matroid,k::IntegerUnion, set::GroundsetType) 
    return k<=rank(M,set) && k<=rank(M,setdiff( Set(M.groundset), set)) && k> connectivity_function(M,set)
end

@doc Markdown.doc"""
    is_k_separation(M::Matroid, k::Int, set::Vector)

Check if `set` together with its complement defines a `k` separation in `M`
See Section 8.1 in Oxl11 (@cite).

# Example
```jldoctest
julia> is_k_separation(fano_matroid(), 2, [1,2,4])
false

```
"""
function is_k_separation(M::Matroid,k::IntegerUnion, set::GroundsetType)
    return k<=length(set) && k<=length(setdiff( Set(M.groundset), set)) && k> connectivity_function(M,set)
end

@doc Markdown.doc"""
    vertical_connectivity(M::Matroid)

If 'M' has two disjoint cocircuits, its vertical connectivity is defined to be least positive integer k such that `M` has a vertical k separation.
Otherwise its vertical connectivity is defined to be the rank of `M`.
See Section 8.6 in Oxl11 (@cite).
# Example
```jldoctest
julia> vertical_connectivity(fano_matroid())
3

```
"""
function vertical_connectivity(M::Matroid)
    gs = Set(M.groundset)
    res = rank(M)
    for set in cocircuits(M)
        comp = setdiff(gs,set)
        if length(set)==0 || length(set)>length(comp)
            continue
        end
        rk_set = rank(M,set)
        rk_comp = rank(M,comp)
        k = minimum([rk_set,rk_comp])
        if k<res && k>rk_set+rk_comp-rank(M)
            res = k
        end
    end
    return res
end

@doc Markdown.doc"""
    girth(M::Matroid, set::Vector)

Return the girth of `set` in the matroid `M`.
This is the size of the smalles circuit contained in `set` and infintie otherwise.
See Section 8.6 in Oxl11 (@cite).

# Example
```jldoctest
julia> girth(fano_matroid(), [1,2,3,4])
3

```
"""
girth(M::Matroid, set::GroundsetType=M.groundset) = minimum([inf; [issubset(C,set) ? length(C) : inf for C in circuits(M)]])

@doc Markdown.doc"""
    tutte_connectivity(M::Matroid)

The Tutte connectivity of `M` is the least integer k such that `M` has a k separation. It can be infinte if no k separation exists.
See Section 8.6 in Oxl11 (@cite).

# Example
```jldoctest
julia> tutte_connectivity(fano_matroid())
3

julia> tutte_connectivity(uniform_matroid(2,4))
PosInf()

```
"""
function tutte_connectivity(M::Matroid)
    r = M.pm_matroid.RANK
    n = M.pm_matroid.N_ELEMENTS
    #if M is uniform, apply Cor. 8.6.3 otherwise Thm. 8.6.4
    if M.pm_matroid.N_BASES==binomial(n,r)
        if n>=2r+2 
            return r+1 
        elseif n<=2r-2 
            n-r+1
        else
            return inf
        end
    end
    return minimum([vertical_connectivity(M),girth(M)])
end

@doc Markdown.doc"""
    tutte_polynomial(M::Matroid)

Return the Tutte polynomial of `M`. This is polynomial in the variables x and y with integral coefficients.
See Section 15.3 in Oxl11 (@cite).

# Example
```jldoctest
julia> tutte_polynomial(fano_matroid())
x^3 + 4*x^2 + 7*x*y + 3*x + y^4 + 3*y^3 + 6*y^2 + 3*y

```
"""
function tutte_polynomial(M::Matroid)
    R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
    poly = M.pm_matroid.TUTTE_POLYNOMIAL
    exp = Polymake.monomials_as_matrix(poly)
    return R(Vector{Int}(Polymake.coefficients_as_vector(poly)),[[exp[i,1],exp[i,2]] for i in 1:size(exp)[1]])
end

@doc Markdown.doc"""
    characteristic_polynomial(M::Matroid)

Return the characteristic polynomial of `M`. This is polynomial in the variable q with integral coefficients.
It is computed as an evaluation of the Tutte polynmomial.
See Section 15.2 in Oxl11 (@cite).

# Example
```jldoctest
julia> characteristic_polynomial(fano_matroid())
q^3 - 7*q^2 + 14*q - 8

```
"""
function characteristic_polynomial(M::Matroid)
    R, q = PolynomialRing(ZZ, 'q')
    return (-1)^M.pm_matroid.RANK*tutte_polynomial(M)(1-q,0)
end
charpoly(M::Matroid) = characteristic_polynomial(M)

@doc Markdown.doc"""
    reduced_characteristic_polynomial(M::Matroid)

Return the reduced characteristic polynomial of `M`. This is the quotient of the characteristic polynomial by (q-1).
See Section 15.2 in Oxl11 (@cite).

# Example
```jldoctest
julia> reduced_characteristic_polynomial(fano_matroid())
q^2 - 6*q + 8

```
"""
function reduced_characteristic_polynomial(M::Matroid)
    R, q = PolynomialRing(ZZ, 'q')
    p = characteristic_polynomial(M)
    c = Vector{Int}(undef,degree(p))
    s = 0
    for i in 1:degree(p)
        s-= coeff(p,i-1)
        c[i] = s
    end
    return R(c)
end

# This function compares two sets A and B in reverse lexicographic order.
# It assumes that both sets are of the same length and ordered.
# It returns true if A is less than B in this order
function revlex_order(A::AbstractVector{Int64}, B::AbstractVector{Int64})
    @assert length(A) == length(B)
    r = length(A)
    if A[end] < B[end]
        return true
    elseif A[end] < B[end]
        return false
    elseif r > 1
        return revlex_order(view(A,1:(r-1)),view(B,1:(r-1)))
    else
        return false
    end
end

#This functions computes the matrix of characteristic vectors of all r element subsets of [n] in revlex order.
function revlex_bases_matrix(r::Int64,n::Int64)
    all_bases = Oscar.Hecke.subsets(Vector(1:n),r)
    sort!(all_bases, lt=revlex_order)
    M = zeros(Int64,length(all_bases),n)
    for i in 1:length(all_bases)
        for j in 1:r
            M[i,all_bases[i][j]] = 1
        end
    end
    # Add columns of ones for polymake
    M = hcat(ones(Int64,length(all_bases)),M)
    return M
end

@doc Markdown.doc"""
Computes the ``revlex basis encoding`` and the ``minimal revlex basis encoding`` among isomorphic matroids 

# Examples
To get the revlex basis encoding of the fano matroid and to preduce a matrod form the encoding write:
```jldoctest
julia> string1, string2 = revlex_basis_encoding(fano_matroid())
("0******0******0***0******0*0**0****", "000****0**0**0***0*****************")

julia> matroid_from_revlex_basis_encoding(string2, 3, 7)
Matroid of rank 3 on 7 elements
```
"""
function revlex_basis_encoding(M::Matroid)
	rvlx = M.pm_matroid.REVLEX_BASIS_ENCODING
	indicies = findall(x->x=='*', rvlx)
	v = zeros(Int,length(rvlx))
	[v[i]=1 for i in indicies]
	hy = Polymake.polytope.hypersimplex(rank(M),size_groundset(M), group=1);
	Polymake.Shell.pair = Polymake.group.lex_minimal(hy.GROUP.VERTICES_ACTION, v)
	Polymake.shell_execute(raw"""$min_v = $pair->first;""")
	return  rvlx, String( [Polymake.Shell.min_v[i]==1 ? '*' : '0' for i in 1:length(rvlx)] )
end
