################################################################################
##  Properties and basic functions
################################################################################


@doc raw"""
    isomorphic_matroid(M::Matroid, gs::GroundsetType)

Return a matroid isomorphic to `M` on the groundset `gs`.

# Examples
```jldoctest
julia> isomorphic_matroid(fano_matroid(), [2,3,4,1,5,6,7])
Matroid of rank 3 on 7 elements
```
"""
function isomorphic_matroid(M::Matroid, gs::GroundsetType)
    if length(M.groundset)!=length(gs) || length(M.groundset)!=length(Set(gs))
        error("Sets of different size")
    end
    gs2num = create_gs2num(gs)
    return Matroid(pm_object(M),gs,gs2num)
end


@doc raw"""
    length(M::Matroid)

Return the size of the ground set of the matroid `M`.

# Examples
```jldoctest
julia> length(fano_matroid())
7
```
"""
length(M::Matroid) = length(M.groundset)::Int


@doc raw"""
    bases([::Type{Int},] M::Matroid)

Return the list of bases of the matroid `M`.
If `Int` is passed as a first argument then the bases will be returned as
indices instead of ground set elements.

# Examples
```jldoctest
julia> bases(uniform_matroid(2, 3))
3-element Vector{Vector{Int64}}:
 [1, 2]
 [1, 3]
 [2, 3]
```
"""
bases(M::Matroid{T}) where T = _indices_to_gs(bases(Int, M), M.groundset)::Vector{Vector{T}}

bases(::Type{Int}, M::Matroid) = _pmset_to_indices(pm_object(M).BASES)::Vector{Vector{Int}}

@doc raw"""
    nonbases(M::Matroid)

Return the list of nonbases of the matroid `M`.

# Examples
```jldoctest
julia> nonbases(fano_matroid())
7-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 4, 5]
 [1, 6, 7]
 [2, 4, 6]
 [2, 5, 7]
 [3, 4, 7]
 [3, 5, 6]
```
"""
nonbases(M::Matroid{T}) where T = _indices_to_gs(nonbases(Int, M), M.groundset)::Vector{Vector{T}}

nonbases(::Type{Int}, M::Matroid) = _pmset_to_indices(pm_object(M).NON_BASES)::Vector{Vector{Int}}


@doc raw"""
    circuits(M::Matroid)

Return the list of circuits of the matroid `M`.

# Examples
```jldoctest
julia> circuits(uniform_matroid(2, 4))
4-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 2, 4]
 [1, 3, 4]
 [2, 3, 4]
```
"""
circuits(M::Matroid{T}) where T = _property_to_gs(M, :CIRCUITS)::Vector{Vector{T}}

@doc raw"""
    hyperplanes(M::Matroid)

Return the list of hyperplanes of the matroid `M`.

# Examples
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
hyperplanes(M::Matroid{T}) where T = _property_to_gs(M, :MATROID_HYPERPLANES)::Vector{Vector{T}}

@doc raw"""
    flats(M::Matroid, [r::Int])

Return the list of flats of the matroid `M`.
By default all flats are returned.
One may specify a rank `r` as the second parameter in which case only the flats of rank `r` are returned.

# Examples
```jldoctest
julia> M = fano_matroid()
Matroid of rank 3 on 7 elements

julia> flats(M)
16-element Vector{Vector{Int64}}:
 []
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
7-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 4, 5]
 [1, 6, 7]
 [2, 4, 6]
 [2, 5, 7]
 [3, 5, 6]
 [3, 4, 7]
```
"""
flats(M::Matroid, r::Union{Int,Nothing}=nothing) = flats_impl(M::Matroid, r::Union{Int,Nothing}, pm_object(M).LATTICE_OF_FLATS.N_NODES,  pm_object(M).LATTICE_OF_FLATS.FACES)

function flats_impl(M::Matroid{T}, r::Union{Int,Nothing}, num_flats::Int, pm_flats) where T
    jl_flats = [Vector{Int}(Polymake.to_one_based_indexing(Polymake._get_entry(pm_flats, i))) for i in 0:(num_flats-1)]
    if pm_object(M).LATTICE_OF_FLATS.TOP_NODE==0
        jl_flats = reverse(jl_flats)
    end
    matroid_flats = _indices_to_gs(jl_flats, M.groundset)
    if r !== nothing
        if r<0 || r>rank(M)
            error("The specified rank needs to be between 0 and the rank of the matroid.")
        end
        matroid_flats = filter(flat -> rank(M,flat)==r, matroid_flats)
    end
    matroid_flats = Vector{Vector{T}}(matroid_flats)
    return matroid_flats
end

@doc raw"""
    cyclic_flats(M::Matroid, [r::Int])

Return the list of cyclic flats of the matroid `M`.
These are the flats that are the union of cycles.
See Section 2.1 in [Oxl11](@cite).

By default all cyclic flats are returned.
One may specify a rank `r` as the second parameter.
In this case only the cyclic flats of this rank are returned.

# Examples
```jldoctest
julia> M = fano_matroid()
Matroid of rank 3 on 7 elements

julia> cyclic_flats(M)
9-element Vector{Vector{Int64}}:
 []
 [1, 2, 3]
 [1, 4, 5]
 [1, 6, 7]
 [2, 4, 6]
 [2, 5, 7]
 [3, 5, 6]
 [3, 4, 7]
 [1, 2, 3, 4, 5, 6, 7]

julia> cyclic_flats(M, 2)
7-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 4, 5]
 [1, 6, 7]
 [2, 4, 6]
 [2, 5, 7]
 [3, 5, 6]
 [3, 4, 7]
```
"""
cyclic_flats(M::Matroid{T}, r::Union{Int,Nothing}=nothing) where T = flats_impl(M, r, pm_object(M).LATTICE_OF_CYCLIC_FLATS.N_NODES,  pm_object(M).LATTICE_OF_CYCLIC_FLATS.FACES)

@doc raw"""
    closure(M::Matroid, set::GroundsetType)

Return the closure of `set` in the matroid `M`.

# Examples
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

@doc raw"""
    rank(M::Matroid)

Return the rank of the matroid `M`.

# Examples
```jldoctest
julia> rank(fano_matroid())
3
```
"""
rank(M::Matroid) = pm_object(M).RANK::Int

@doc raw"""
    rank(M::Matroid, set::GroundsetType)

Return the rank of `set` in the matroid `M`.

# Examples
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
        return Polymake.matroid.rank(pm_object(M), _gs_to_pmindices(set, M.gs2num; type=Set))::Int
    end
end

@doc raw"""
    nullity(M::Matroid, set::GroundsetType)

Return the nullity of `set` in the matroid `M`.
This is defined to be `|set| - rk(set)`.

# Examples
```jldoctest
julia> M = fano_matroid();

julia> nullity(M, [1,2,3])
1
```
"""
nullity(M::Matroid, set::GroundsetType) = length(set)-Polymake.matroid.rank(pm_object(M), _gs_to_pmindices(set, M.gs2num; type=Set))::Int


@doc raw"""
    fundamental_circuit(M::Matroid, basis::GroundsetType, elem::ElementType)

Return the unique circuit contained in the union of `basis` and `elem` of the matroid `M`.
See Section 1.2 of [Oxl11](@cite).
Note that `elem` needs to be in the complement of the `basis` in this case.


# Examples
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
    if !(Set(basis) in Set.(bases(M)))
        error("The set is not a basis of M")
    end
    if !(elem in M.groundset)
        error("The element is not in the groundset")
    end
    if elem in basis
        error("The element is in the basis, but has to be in the complement")
    end
    return circuits(restriction(M,[elem; collect(basis)]))[1]
end

@doc raw"""
    fundamental_cocircuit(M::Matroid, cobasis::GroundsetType, elem::ElementType)

Return the unique circuit of the dual matroid of `M` in the union of the complement of `basis` and `elem`.
See Section 2.1 of [Oxl11](@cite).
Note that `elem` needs to be an element of the `basis` in this case.


# Examples
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

@doc raw"""
    independent_sets(M::Matroid)

Return the list of independent sets of the matroid `M`.
These are all subsets of the bases.

# Examples
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
    pm_bases = bases(Int, M)
    gs = matroid_groundset(M)
    elt = eltype(gs)
    n = length(gs)
    sets = Vector{Vector{elt}}()
    push!(sets,elt[])
    for k in 1:rank(M)
        for set in Oscar.Hecke.subsets(Vector(1:n),k)
            for B in pm_bases
                if issubset(set,B)
                    push!(sets, gs[collect(set)])
                    break
                end
            end
        end
    end
    return sets
end

@doc raw"""
    spanning_sets(M::Matroid)

Return the list of spanning sets of the matroid `M`.
These are all sets containing a basis.

# Examples
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
    span_sets = [filter(!in(set), matroid_groundset(M)) for set in coindependent_sets]
    return reverse(span_sets)
end

@doc raw"""
    cobases(M::Matroid)

Return the bases of the dual matroid of `M`.
See Section 2 in [Oxl11](@cite).

# Examples
```jldoctest
julia> cobases(uniform_matroid(2, 3))
3-element Vector{Vector{Int64}}:
 [3]
 [2]
 [1]
```
"""
cobases(M::Matroid) = bases(dual_matroid(M))

@doc raw"""
    cocircuits(M::Matroid)

Return the circuits of the dual matroid of `M`.
See Section 2 in [Oxl11](@cite).

# Examples
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

@doc raw"""
    cohyperplanes(M::Matroid)

Return the hyperplanes of the dual matroid of `M`.
See Section 2 in [Oxl11](@cite).

# Examples
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

@doc raw"""
    corank(M::Matroid, set::GroundsetType)

Return the rank of `set` in the dual matroid of `M`.

# Examples
```jldoctest
julia> corank(fano_matroid(), [1,2,3])
3
```
"""
corank(M::Matroid, set::GroundsetType) = length(set)-rank(M, M.groundset) + rank(M, setdiff(M.groundset, set))

@doc raw"""
    is_clutter(sets::AbstractVector{T}) where T <: GroundsetType

Check if the collection of subsets `sets` is a clutter.
A collection of subsets is a clutter if none of the sets is a proper subset of another.
See Section 2.1 in [Oxl11](@cite).

# Examples
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

@doc raw"""
    is_regular(M::Matroid)

Check if the matroid `M` is regular, that is representable over every field.
See Section 6.6 in [Oxl11](@cite).

# Examples
```jldoctest
julia> is_regular(uniform_matroid(2, 3))
true

julia> is_regular(fano_matroid())
false
```
"""
function is_regular(M::Matroid)
    #if statement to avoid bug in polymake (TODO remove when fixed)
    if rank(M)==length(M)
        return true
    end
    return pm_object(M).REGULAR::Bool
end

@doc raw"""
    is_binary(M::Matroid)

Check whether the matroid `M` is binary, that is representable over the finite field `F_2`.
See Section 6.5 in [Oxl11](@cite).

# Examples
```jldoctest
julia> is_binary(uniform_matroid(2, 4))
false

julia> is_binary(fano_matroid())
true
```
"""
is_binary(M::Matroid) = pm_object(M).BINARY::Bool

@doc raw"""
    is_ternary(M::Matroid)

Check if the matroid `M` is ternary, that is representable over the finite field `F_3`.
See Section 4.1 in [Oxl11](@cite).

# Examples
```jldoctest
julia> is_ternary(uniform_matroid(2, 4))
true

julia> is_ternary(fano_matroid())
false
```
"""
function is_ternary(M::Matroid)
    #if statement to avoid bug in polymake (TODO remove when fixed)
    if rank(M)==length(M)
        return true
    end
    return pm_object(M).TERNARY::Bool
end

@doc raw"""
    n_connected_components(M::Matroid)

Return the number of connected components of `M`.
See Section 4.1 in [Oxl11](@cite).

# Examples
```jldoctest
julia> n_connected_components(fano_matroid())
1

julia> n_connected_components(uniform_matroid(3, 3))
3
```
"""
n_connected_components(M::Matroid) = length(pm_object(M).CONNECTED_COMPONENTS)::Int

@doc raw"""
    connected_components(M::Matroid)

Return the connected components of `M`. The function returns a partition of the ground set where each part corresponds to one connected component. 
See Section 4.1 in [Oxl11](@cite).

# Examples
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
connected_components(M::Matroid{T}) where T = _property_to_gs(M, :CONNECTED_COMPONENTS)::Vector{Vector{T}}

@doc raw"""
    is_connected(M::Matroid)

Check if the matroid `M` is connected, that is has one connected component
See Section 4.1 in [Oxl11](@cite).

# Examples
```jldoctest
julia> is_connected(fano_matroid())
true

julia> is_connected(uniform_matroid(3, 3))
false
```
"""
is_connected(M::Matroid) = pm_object(M).CONNECTED::Bool

@doc raw"""
    loops(M::Matroid)

Return the loops of `M`. A loop is an element of the ground set that is not contained in any basis.

# Examples
```jldoctest
julia> loops(matroid_from_bases([[1,2]], 4))
2-element Vector{Int64}:
 3
 4

julia> loops(fano_matroid())
Int64[]
```
"""
loops(M::Matroid{T}) where T = _property_to_gs(M, :LOOPS)::Vector{T}

@doc raw"""
    coloops(M::Matroid)

Return the coloops of `M`. A coloop is an element of the ground set that is contained in every basis.

# Examples
```jldoctest
julia> coloops(matroid_from_bases([[1,2]], 4))
2-element Vector{Int64}:
 1
 2

julia> coloops(fano_matroid())
Int64[]
```
"""
coloops(M::Matroid{T}) where T = _property_to_gs(M, Symbol("DUAL.LOOPS"))::Vector{T}

@doc raw"""
    is_loopless(M::Matroid)

Check if `M` has a loop. Return `true` if `M` does not have a loop.
See also `loops`.

# Examples
```jldoctest
julia> is_loopless(matroid_from_bases([[1,2]], 4))
false

julia> is_loopless(fano_matroid())
true
```
"""
is_loopless(M::Matroid) = length(pm_object(M).LOOPS)==0 ? true : false

@doc raw"""
    is_coloopless(M::Matroid)

Check if `M` has a coloop. Return `true` if `M` does not have a coloop.
See also `coloops`.

# Examples
```jldoctest
julia> is_coloopless(matroid_from_bases([[1,2]], 4))
false

julia> is_coloopless(fano_matroid())
true
```
"""
is_coloopless(M::Matroid) = length(pm_object(M).DUAL.LOOPS)==0 ? true : false

@doc raw"""
    is_simple(M::Matroid)

Check if `M` has is simple. A matroid is simple if it doesn't have loops and doesn't have parallel elements.
Return `true` if `M` is simple.
See also `loops`.

# Examples
```jldoctest
julia> is_simple(matroid_from_bases([[1,2]], 4))
false

julia> is_simple(fano_matroid())
true
```
"""
is_simple(M::Matroid) = pm_object(M).SIMPLE::Bool

@doc raw"""
    direct_sum_components(M::Matroid)

Return the connected components of `M` as a list of matroids.
See Section 4.1 in [Oxl11](@cite).

# Examples
```jldoctest
julia> direct_sum_components(fano_matroid())
1-element Vector{Matroid}:
 Matroid of rank 3 on 7 elements

julia> direct_sum_components(uniform_matroid(3, 3))
3-element Vector{Matroid}:
 Matroid of rank 1 on 1 element
 Matroid of rank 1 on 1 element
 Matroid of rank 1 on 1 element
```
"""
function direct_sum_components(M::Matroid)
    res = Vector{Matroid}()
    for set in connected_components(M)
        res = [res; restriction(M,set)]
    end
    return res
end

@doc raw"""
    connectivity_function(M::Matroid, set::GroundsetType)

Return the value of the connectivity function of `set` in the matroid `M`.
See Section 8.1 in [Oxl11](@cite).

# Examples
```jldoctest
julia> connectivity_function(fano_matroid(), [1,2,4])
3

```
"""
function connectivity_function(M::Matroid, set::GroundsetType)
    return rank(M,set) + rank(M,setdiff(Set(M.groundset), set)) - rank(M)
end

@doc raw"""
    is_vertical_k_separation(M::Matroid, k::IntegerUnion, set::GroundsetType)

Check if `set` together with its complement defines a `k` separation in `M`
See Section 8.6 in [Oxl11](@cite).

# Examples
```jldoctest
julia> is_vertical_k_separation(fano_matroid(), 2, [1,2,4])
false

```
"""
function is_vertical_k_separation(M::Matroid,k::IntegerUnion, set::GroundsetType) 
    return k<=rank(M,set) && k<=rank(M,setdiff( Set(M.groundset), set)) && k> connectivity_function(M,set)
end

@doc raw"""
    is_k_separation(M::Matroid, k::IntegerUnion, set::GroundsetType)

Check if `set` together with its complement defines a `k` separation in `M`
See Section 8.1 in [Oxl11](@cite).

# Examples
```jldoctest
julia> is_k_separation(fano_matroid(), 2, [1,2,4])
false

```
"""
function is_k_separation(M::Matroid,k::IntegerUnion, set::GroundsetType)
    return k<=length(set) && k<=length(setdiff( Set(M.groundset), set)) && k> connectivity_function(M,set)
end

@doc raw"""
    vertical_connectivity(M::Matroid)

If 'M' has two disjoint cocircuits, its vertical connectivity is defined to be least positive integer k such that `M` has a vertical k separation.
Otherwise its vertical connectivity is defined to be the rank of `M`.
See Section 8.6 in [Oxl11](@cite).
# Examples
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

@doc raw"""
    girth(M::Matroid, set::GroundsetType)

Return the girth of `set` in the matroid `M`.
This is the size of the smallest circuit contained in `set` and infinite otherwise.
See Section 8.6 in [Oxl11](@cite).

# Examples
```jldoctest
julia> girth(fano_matroid(), [1,2,3,4])
3

```
"""
girth(M::Matroid, set::GroundsetType=M.groundset) = minimum(issubset(C,set) ? length(C) : inf for C in circuits(M); init=inf)

@doc raw"""
    tutte_connectivity(M::Matroid)

The Tutte connectivity of `M` is the least integer k such that `M` has a k separation. It can be infinite if no k separation exists.
See Section 8.6 in [Oxl11](@cite).

# Examples
```jldoctest
julia> tutte_connectivity(fano_matroid())
3

julia> tutte_connectivity(uniform_matroid(2,4))
infinity

```
"""
function tutte_connectivity(M::Matroid)
    r = rank(M)
    n = length(M)
    #if M is uniform, apply Cor. 8.6.3 otherwise Thm. 8.6.4
    if pm_object(M).N_BASES==binomial(n,r)
        if n>=2r+2 
            return r+1 
        elseif n<=2r-2 
            n-r+1
        else
            return inf
        end
    end
    return min(vertical_connectivity(M), girth(M))
end

@doc raw"""
    tutte_polynomial(M::Matroid)
    tutte_polynomial(M::Matroid; parent::ZZMPolyRing)
    tutte_polynomial(parent::ZZMPolyRing, M::Matroid)

Return the Tutte polynomial of `M`. This is polynomial in the variables x and y with integral coefficients.
See Section 15.3 in [Oxl11](@cite).

# Examples
```jldoctest
julia> tutte_polynomial(fano_matroid())
x^3 + 4*x^2 + 7*x*y + 3*x + y^4 + 3*y^3 + 6*y^2 + 3*y

```
"""
function tutte_polynomial(M::Matroid;
           parent::ZZMPolyRing = polynomial_ring(ZZ, [:x, :y]; cached = false)[1])
  @assert ngens(parent) >= 2
  poly = pm_object(M).TUTTE_POLYNOMIAL
  coeffs = Vector{ZZRingElem}(Polymake.coefficients_as_vector(poly))
  exp = Polymake.monomials_as_matrix(poly)
  ev = [[exp[i, 1],exp[i, 2]] for i in 1:size(exp)[1]]::Vector{Vector{Int}}
  return parent(coeffs, ev)
end

tutte_polynomial(R::ZZMPolyRing, M::Matroid) = tutte_polynomial(M, parent = R)

@doc raw"""
    characteristic_polynomial(M::Matroid)
    characteristic_polynomial(M::Matroid; parent::ZZPolyRing)
    characteristic_polynomial(parent::ZZPolyRing, M::Matroid)

Return the characteristic polynomial of `M`. This is polynomial in the variable q with integral coefficients.
It is computed as an evaluation of the Tutte polynmomial.
See Section 15.2 in [Oxl11](@cite).

# Examples
```jldoctest
julia> characteristic_polynomial(fano_matroid())
q^3 - 7*q^2 + 14*q - 8

```
"""
function characteristic_polynomial(M::Matroid;
           parent::ZZPolyRing = polynomial_ring(ZZ, :q; cached = false)[1])
  return (-1)^rank(M) * tutte_polynomial(M)(1 - gen(parent), 0)::ZZPolyRingElem
end

characteristic_polynomial(R::ZZPolyRing, M::Matroid) = characteristic_polynomial(M, parent = R)

@doc raw"""
    reduced_characteristic_polynomial(M::Matroid)
    reduced_characteristic_polynomial(M::Matroid; parent::ZZPolyRing)
    reduced_characteristic_polynomial(parent::ZZPolyRing, M::Matroid)

Return the reduced characteristic polynomial of `M`. This is the quotient of the characteristic polynomial by (q-1).
See Section 15.2 in [Oxl11](@cite).

# Examples
```jldoctest
julia> reduced_characteristic_polynomial(fano_matroid())
q^2 - 6*q + 8

```
"""
function reduced_characteristic_polynomial(M::Matroid;
           parent::ZZPolyRing = polynomial_ring(ZZ, :q; cached = false)[1])
  p = characteristic_polynomial(M, parent = parent)
  c = Vector{ZZRingElem}(undef, degree(p))
  s = ZZ(0)
  for i in 1:degree(p)
    s -= coeff(p, i - 1)
    c[i] = s
  end
  return parent(c)
end

reduced_characteristic_polynomial(R::ZZPolyRing, M::Matroid) = reduced_characteristic_polynomial(M, parent = R)

# This function compares two sets A and B in reverse lexicographic order.
# It assumes that both sets are of the same length and ordered.
# It returns true if A is less than B in this order
function revlex_order(A::AbstractVector{Int64}, B::AbstractVector{Int64})
    @assert length(A) == length(B)
    r = length(A)
    for i in 0:r-1
        if A[r-i] < B[r-i]
            return true
        elseif A[r-i] > B[r-i]
            return false
        end
    end
    return false
end

#This functions computes the matrix of characteristic vectors of all r element subsets of [n] in revlex order.
function revlex_bases_matrix(r::Int64,n::Int64)
    all_bases = subsets(Vector(1:n),r)
    sort!(all_bases, lt=revlex_order)
    M = zeros(Int64,length(all_bases),n)
    for i in 1:length(all_bases), j in 1:r
        M[i,all_bases[i][j]] = 1
    end
    # Add columns of ones for polymake
    M = hcat(ones(Int64,length(all_bases)),M)
    return M
end

_revlex_basis_to_vector(s::AbstractString) = Int[x=='*' ? 1 : 0 for x in s]
_revlex_basis_from_vector(v::AbstractVector{<:Integer}) = join(isone(x) ? '*' : '0' for x in v)

@doc raw"""
    revlex_basis_encoding(M::Matroid)

Compute the revlex basis encoding of the matroid M.

# Examples
To get the revlex basis encoding of the fano matroid and to produce a matrod form the encoding write:
```jldoctest
julia> str = revlex_basis_encoding(fano_matroid())
"0******0******0***0******0*0**0****"

julia> matroid_from_revlex_basis_encoding(str, 3, 7)
Matroid of rank 3 on 7 elements

```
"""
function revlex_basis_encoding(M::Matroid)
    return String(pm_object(M).REVLEX_BASIS_ENCODING)
end

@doc raw"""
    min_revlex_basis_encoding(M::Matroid)

Compute the minimal revlex basis encoding among isomorphic matroids.

# Examples
To get the minimal revlex basis encoding of the fano matroid write:
```jldoctest
julia> str = min_revlex_basis_encoding(fano_matroid())
"0******0******0***0******0*0**0****"

```
"""
function min_revlex_basis_encoding(M::Matroid)
    rvlx = pm_object(M).REVLEX_BASIS_ENCODING
    v = _revlex_basis_to_vector(rvlx)
    n = length(M)
    A = revlex_bases_matrix(rank(M),n)
    if n<2
      gens = [collect(0:n-1)]
    else
      gens = [[[1, 0]; collect(2:n-1)], [[n-1]; collect(0:n-2)]]
    end
    vertperm = Polymake.group.induced_permutations(gens, A)
    action = Polymake.group.PermutationAction(GENERATORS=vertperm)
    all_elem = Polymake.group.all_group_elements(action)
    st = Polymake.SwitchTable(all_elem)
    min_rvlx = first(Polymake.lex_minimize_vector(st, convert(Polymake.PolymakeType, v)))
    return _revlex_basis_from_vector(min_rvlx)
end

@doc raw"""
    matroid_hex(M::Matroid)

Store a matroid as a string of hex characters. The first part of the string is
"r" followed by the rank of the matroid. This is followed by "n" and the
number of elements. The rest of the string is the revlex basis encoding. The
encoding is done by converting the basis encoding to a vector of bits and then
to a string of characters. The bits are padded to a multiple of 4 and then
converted to hex characters.

# Examples
To get the hex encoding of the fano matroid write:
```jldoctest
julia> matroid_hex(fano_matroid())
"r3n7_3f7eefd6f"

```
"""
function matroid_hex(M::Matroid)
  rvlx = min_revlex_basis_encoding(M)
  r,n = rank(M), length(M) 
  v = zeros(Int, 4*ceil(Int, length(rvlx)/4))
  v[length(v)-length(rvlx)+1:end] = _revlex_basis_to_vector(rvlx)

  v = reshape(v,4,:)
  v = [string(parse(Int, join(v[:, j]), base=2), base=16) for j in 1:size(v)[2]]

  return "r$(r)n$(n)_" * join(v)
end

@doc raw"""
    matroid_from_matroid_hex(str::AbstractString)

Return a matroid from a string of hex characters.

# Examples
To retrieve the fano matroid from its hex encoding write:

```jldoctest
julia> matroid_from_matroid_hex("r3n7_3f7eefd6f")
Matroid of rank 3 on 7 elements

```
"""
function matroid_from_matroid_hex(str::AbstractString)
  @req occursin(r"^r\d+n\d+_[0-9a-f]+$", str) "Invalid hex encoding"

  sep = split(str, "_")
  (r,n) = parse.(Int,split(sep[1][2:end],"n"))

  v = [digits(parse(Int, x, base=16), base=2, pad=4) |> reverse for x in sep[2]]
  v = foldl(append!, v)
  v = v[(length(v)-binomial(n, r)+1):end]

  return matroid_from_revlex_basis_encoding(_revlex_basis_from_vector(v), r, n)
end

@doc raw"""
    is_isomorphic(M1::Matroid, M2::Matroid)

Check if the matroid `M1` is isomorphic to the matroid `M2` under the action of the symmetric group that acts on their groundsets.

# Examples
To compare two matrods write:
```jldoctest
julia> H = [[1,2,4],[2,3,5],[1,3,6],[3,4,7],[1,5,7],[2,6,7],[4,5,6]];

julia> M = matroid_from_hyperplanes(H,7);

julia> is_isomorphic(M,fano_matroid())
true

```
"""
function is_isomorphic(M1::Matroid, M2::Matroid)
    if length(M1) != length(M2)
        return false
    end
    return Polymake.matroid.is_isomorphic_to(M1.pm_matroid, M2.pm_matroid)::Bool
end

@doc raw"""
    is_minor(M::Matroid, N::Matroid)

Check if the matroid `M` is isomorphic to a minor of the matroid `N`.

# Examples
```jldoctest
julia> is_minor(direct_sum(uniform_matroid(0,1), uniform_matroid(2,2)), fano_matroid())
false

julia> is_minor(direct_sum(uniform_matroid(0,1), uniform_matroid(2,2)), parallel_extension(uniform_matroid(3,4), 1, 5))
true
```
"""
function is_minor(Minor::Matroid, M::Matroid)
    n = length(M)
    c = rank(M)-rank(Minor)
    d = n-length(Minor)-c
    nB = length(bases(Minor))
    if(c<0 || d<0 || length(bases(M))< nB)
        return false
    end
    basesM = bases(Int, M)
    basesMinor = bases(Int, Minor)

    for set_C in Oscar.Hecke.subsets(Vector(1:n), c) # set to contract
        bases_C = filter(B->issubset(set_C,B), basesM)
        for set_R in Oscar.Hecke.subsets(setdiff(Set(1:n), set_C), n-c-d) # set t restrict on
            bases_R = filter(B->issubset(B,union(set_R, set_C)), bases_C)
            if length(bases_R) == nB
                I = Polymake.IncidenceMatrix(nB, n, bases_R)[1:nB, collect(set_R)]
                Iminor = Polymake.IncidenceMatrix(nB, length(Minor), basesMinor)
                if !isnothing(Polymake.graph.find_row_col_permutation(I,Iminor))
                    return true
                end
            end
        end
    end
    return false
end

@doc raw"""
    matroid_base_polytope(M::Matroid) 

The base polytope of the matroid `M`.  

# Examples
```jldoctest
julia> D = matroid_base_polytope(uniform_matroid(2,4));

julia> vertices(D)
6-element SubObjectIterator{PointVector{QQFieldElem}}:
 [1, 1, 0, 0]
 [1, 0, 1, 0]
 [1, 0, 0, 1]
 [0, 1, 1, 0]
 [0, 1, 0, 1]
 [0, 0, 1, 1]
```
"""
function matroid_base_polytope(M::Matroid)
    n = length(matroid_groundset(M))
    M = isomorphic_matroid(M, [i for i in 1:n])
    Delta_verts = hcat([indicator_vector(x, n) for x in bases(M)]...)
    return convex_hull(Delta_verts') 
end


function indicator_vector(S::Vector{Int}, n::Int)
    return map(x -> x in S ? 1 : 0 , 1:n)
end
