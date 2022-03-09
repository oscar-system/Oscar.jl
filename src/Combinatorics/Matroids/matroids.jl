export
    Matroid, matroid_groundset, 
    matroid_from_revlex_basis_encoding,
    matroid_from_bases, matroid_from_nonbases, matroid_from_circuits, matroid_from_hyperplanes,
    matroid_from_matrix_columns, matroid_from_matrix_rows,
    cycle_matroid, bond_matroid, cocycle_matroid,
    dual_matroid, direct_sum, restriction, deletion, contraction, minor,
    principal_extension, free_extension, series_extension, parallel_extension,
    uniform_matroid, fano_matroid, non_fano_matroid, pappus_matroid, non_pappus_matroid, vamos_matroid,
    all_subsets_matroid, projective_plane, projective_geometry, affine_geometry


include("representations.jl")

################################################################################
##  Constructing
################################################################################
ElementType = Union{IntegerUnion,Char,String}

struct Matroid
    pm_matroid::Polymake.BigObject
    groundset::Vector{ElementType} # groundset of the matroid 
    gs2num::Dict{Any, IntegerUnion}# dictionary to map the groundset to the integers from 1 to its size
end

pm_object(M::Matroid) = M.pm_matroid

function Base.show(io::IO, M::Matroid)
    r = rank(M)
    n = length(M.groundset)
    print(io, "Matroid of rank $(r) on $(n) elements")
end

"""
Matroid(M)
Construct a `matroid` from a ``polymake matroid M`` on the default ground set `{1,...,n}`.
"""
function Matroid(pm_matroid::Polymake.BigObjectAllocated)
    gs2num = Dict{Any,IntegerUnion}()
    [gs2num[i] = i for i in 1:pm_matroid.N_ELEMENTS]
    return Matroid(pm_matroid, 1:pm_matroid.N_ELEMENTS, gs2num)
end

"""
matroid_from_revlex_encoding(M, r, n)
Construct a `matroid` from a ``polymake matroid M`` of rank `r` on the default ground set `{1,...,n}`.
"""
function matroid_from_revlex_basis_encoding(rvlx::String, r::IntegerUnion, n::IntegerUnion)
    if match(r"[^*0]",rvlx)!=nothing
	    error("The revlex encoding uses only `*` and `0`")
    end
    if length(rvlx)!= binomial(n,r)
	    error("The length of the string does not match the rank and number of elements")
    end
    return Matroid(Polymake.matroid.Matroid(N_ELEMENTS=n, RANK=r, REVLEX_BASIS_ENCODING=rvlx))
end

@doc Markdown.doc"""
    matroid_from_bases(B, [n, E])

# Arguments
- `B::AbstractVector`: The set of bases of the matroid.
- `n::InterUnion`: The size of the ground set. The ground set will be `{1,..n}` in this case.
- `E::AbstractVector`: An explicit ground set passed as vector.

Construct a `matroid` with bases `B` on the ground set `E` (which can be the empty set).
The set `B` is a non-empty collection of subsets of the ground set `E` satisfying an exchange property,
and the default value for `E` is the set `{1,..n}` for a non-negative value `n`.

See Section 1.2 of Oxl11 (@cite).

# Examples
To construct a rank two matroid with five bases on four elements you can write:
```jldoctest
julia> B = [[1,2],[1,3],[1,4],[2,3],[2,4]];

julia> M = matroid_from_bases(B,4)
Matroid of rank 2 on 4 elements
```

# Examples
To construct the same matroid on the four elements 1,2,i,j you may write:
```jldoctest
julia> M = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j'])
Matroid of rank 2 on 4 elements
```
"""
matroid_from_bases(bases::Union{AbstractVector{<:AbstractVector{<:IntegerUnion}}, AbstractVector{<:AbstractSet{<:IntegerUnion}}}, nelements::IntegerUnion; check::Bool=true) = matroid_from_bases(bases,Vector(1:nelements);check=check)

function matroid_from_bases(bases::Union{AbstractVector{<:AbstractVector}, AbstractVector{<:AbstractSet}}, groundset::AbstractVector; check::Bool=true)
    if check && length(groundset)!=length(Set(groundset))
        error("Input is not a valid groundset of a matroid")
    end
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in groundset
        gs2num[elem] = i
        i+=1
    end
    pm_bases = [[gs2num[i]-1 for i in B] for B in bases]
    M = Polymake.matroid.Matroid(BASES=pm_bases,N_ELEMENTS=length(groundset))
    if check && !Polymake.matroid.check_basis_exchange_axiom(M.BASES)
        error("Input is not a collection of bases")
    end
    return Matroid(M,groundset,gs2num)
end

@doc Markdown.doc"""
    matroid_from_nonbases(N, [n, E])

# Arguments
- `N::AbstractVector`: The set of nonbases of the matroid.
- `n::InterUnion`: The size of the ground set. The ground set will be `{1,..n}` in this case.
- `E::AbstractVector`: An explicit ground set passed as vector.

Construct a `matroid` with nonbases `N` on the ground set `E` (which can be the empty set).
That means that the matroid has as bases all subsets of the size `|N[1]|` of the ground set that are not in `N`.
The set `N` can't be empty in this function.
The described complement of `N` needs to be a non-empty collection of subsets of the ground set `E` satisfying an exchange property,
and the default value for `E` is the set `{1,..n}` for a non-negative value `n`.

See Section 1.2 of Oxl11 (@cite).

# Examples
To construct the Fano matroid you may write:
```jldoctest
julia> H = [[1,2,4],[2,3,5],[1,3,6],[3,4,7],[1,5,7],[2,6,7],[4,5,6]];

julia> M = matroid_from_nonbases(H,7)
Matroid of rank 3 on 7 elements

```
"""
matroid_from_nonbases(nonbases::Union{AbstractVector{<:AbstractVector{<:IntegerUnion}}, AbstractVector{<:AbstractSet{<:IntegerUnion}}}, nelements::IntegerUnion; check::Bool=true) = matroid_from_nonbases(nonbases,Vector(1:nelements); check)

function matroid_from_nonbases(nonbases::Union{AbstractVector{<:AbstractVector}, AbstractVector{<:AbstractSet}},groundset::AbstractVector; check::Bool=true)
    if check && length(groundset)!=length(Set(groundset))
        error("Input is not a valid groundset of a matroid")
    end
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in groundset
        gs2num[elem] = i
        i+=1
    end
    if length(nonbases)==0
        error("The collection of nonbases should not be empty.")
    end

    rk = length(nonbases[1])
    Bases = filter(set -> set ∉ nonbases, Oscar.Hecke.subsets(Vector(1:length(groundset)),rk))
    pm_bases = [[gs2num[i]-1 for i in B] for B in Bases]
    M = Polymake.matroid.Matroid(BASES=pm_bases,N_ELEMENTS=length(groundset))
    if check && !Polymake.matroid.check_basis_exchange_axiom(M.BASES)
        error("Input is not a collection of nonbases")
    end
    return Matroid(M,groundset,gs2num)
end

@doc Markdown.doc"""
    matroid_from_circuits(C, [n, E])

# Arguments
- `C::AbstractVector`: The set of circuits of the matroid.
- `n::InterUnion`: The size of the ground set. The ground set will be `{1,..n}` in this case.
- `E::AbstractVector`: An explicit ground set passed as vector.

A matroid with circuits `C` on the ground set `E` (which can be the empty set).
The set `C` is a collection of subsets of the ground set `E` satisfying an exchange property,
and the default value for `E` is the set `{1,..n}` for a non-negative value `n`. 

See Section 1.1 of Oxl11 (@cite).

# Examples
To construct a rank two matroid with five bases on four elements by its circuits you may write:
```jldoctest
julia> C = [[1,2,3],[1,2,4],[3,4]];

julia> M = matroid_from_circuits(C,4)
Matroid of rank 2 on 4 elements
```

# Examples
To construct the same matroid on the ground set `{1,2,i,j}` you may write:
```jldoctest
julia> C = [[1,2,'j'],[1,2,'i'],['i','j']];

julia> M = matroid_from_circuits(C,[1,2,'i','j'])
Matroid of rank 2 on 4 elements
```
"""
matroid_from_circuits(circuits::Union{AbstractVector{<:AbstractVector{<:IntegerUnion}}, AbstractVector{<:AbstractSet{<:IntegerUnion}}}, nelements::IntegerUnion) = matroid_from_circuits(circuits,Vector(1:nelements))

function matroid_from_circuits(circuits::Union{AbstractVector{<:AbstractVector}, AbstractVector{<:AbstractSet}},groundset::AbstractVector; check::Bool=true)
    if check && length(groundset)!=length(Set(groundset))
        error("Input is not a valid groundset of a matroid")
    end
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in groundset
        gs2num[elem] = i
        i+=1
    end
    pm_circuits = [[gs2num[i]-1 for i in C] for C in circuits]
    M = Polymake.matroid.Matroid(CIRCUITS=pm_circuits,N_ELEMENTS=length(groundset))
    #TODO check_circuit_exchange_axiom (requires an update of polymake)
    #if check && !Polymake.matroid.check_circuit_exchange_axiom(M.CIRCUITS)
    #   error("Input is not a collection of circuits")
    #end
    return Matroid(M,groundset,gs2num)
end

@doc Markdown.doc"""
    matroid_from_hyperplanes(H, [n, E])

# Arguments
- `H::AbstractVector`: The set of hyperplanes of the matroid.
- `n::InterUnion`: The size of the ground set. The ground set will be `{1,..n}` in this case.
- `E::AbstractVector`: An explicit ground set passed as vector.

A matroid with hyperplanes `H` on the ground set `E` (which can be the empty set).
A hyperplane is a flat of rank `r-1`.
The set `H` is a collection of subsets of the ground set `E` satisfying an exchange property,
and the default value for `E` is the set `{1,..n}` for a non-negative value `n`. 

See Section 1.4 of Oxl11 (@cite).

# Examples
To construct the Fano matroid you may write:
```jldoctest
julia> H = [[1,2,4],[2,3,5],[1,3,6],[3,4,7],[1,5,7],[2,6,7],[4,5,6]];

julia> M = matroid_from_hyperplanes(H,7)
Matroid of rank 3 on 7 elements

```
"""
matroid_from_hyperplanes(hyperplanes::Union{AbstractVector{<:AbstractVector{<:IntegerUnion}}, AbstractVector{<:AbstractSet{<:IntegerUnion}}}, nelements::IntegerUnion) = matroid_from_hyperplanes(hyperplanes,Vector(1:nelements))

function matroid_from_hyperplanes(hyperplanes::Union{AbstractVector{<:AbstractVector}, AbstractVector{<:AbstractSet}},groundset::AbstractVector; check::Bool=true)
    if check && length(groundset)!=length(Set(groundset))
        error("Input is not a valid groundset of a matroid")
    end
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in groundset
        gs2num[elem] = i
        i+=1
    end
    pm_hyperplanes = [[gs2num[i]-1 for i in H] for H in hyperplanes]
    M = Polymake.matroid.Matroid(MATROID_HYPERPLANES=pm_hyperplanes,N_ELEMENTS=length(groundset))
    #TODO implement a check if these are actually the hyperplanes of a matroid
    return Matroid(M,groundset,gs2num)
end

@doc Markdown.doc"""
    matroid_from_matrix_columns(A::MatrixElem)

A matroid represented by the column vectors of a matrix `A`.

See Section 1.1 of Oxl11 (@cite).

# Examples
To construct the vector matroid (a.k.a linear matroid) of the matrix `A` over the field with two elements write:
```jldoctest
julia> A = matrix(GF(2),[[1,0,1,1],[0,1,1,1]]);

julia> M = matroid_from_matrix_columns(A)
Matroid of rank 2 on 4 elements
```
"""
function matroid_from_matrix_columns(A::MatrixElem; check=true)
    rk = rank(A)
    nr = nrows(A)
    bases = Vector{Vector{IntegerUnion}}()
    # We use a temporary matrix which we will rewrite to avoid creating new matrices in every iteration.
    tmp_mat = A[:,1:rk]
    for set in subsets(Vector(1:ncols(A)),rk)
        for is in 1:rk
            for ir in 1:nr
                tmp_mat[ir,is] = A[ir,set[is]]
            end
        end
        if rank(tmp_mat)==rk
            push!(bases, set);
        end
    end
    return matroid_from_bases(bases, ncols(A); check=check)
end

@doc Markdown.doc"""
    matroid_from_matrix_columns(A::MatrixElem)

A matroid represented by the row vectors of a matrix.

See Section 1.1 of Oxl11 (@cite).

# Examples
To construct the linear matroid of the rows of the matrix `A` over the field with two elements write:
```jldoctest
julia> A = matrix(GF(2),[[1,0],[0,1],[1,1],[1,1]]);

julia> M = matroid_from_matrix_rows(A)
Matroid of rank 2 on 4 elements
```
"""
matroid_from_matrix_rows(A::MatrixElem, ; check=true) = matroid_from_matrix_columns(transpose(A); check=check)

@doc Markdown.doc"""
    cycle_matroid(g::Oscar.Graphs.Graph)

The cycle matroid of a graph `g`.

See Section 1.1 of Oxl11 (@cite).

# Examples
To construct the cycle matroid of the complete graph of 4 vertices write:
```jldoctest
julia> g = Oscar.Graphs.complete_graph(4);

julia> M = cycle_matroid(g)
Matroid of rank 3 on 6 elements
```
"""
function cycle_matroid(g::Oscar.Graphs.Graph)
    pm_Graph = Polymake.graph.Graph(ADJACENCY=g.pm_graph)
    M = Polymake.matroid.matroid_from_graph(pm_Graph)
    n = Oscar.Graphs.ne(g)
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in 1:n
        gs2num[elem] = i
        i+=1
    end
    return Matroid(M,1:n,gs2num)
end

@doc Markdown.doc"""
    bond_matroid(g::Oscar.Graphs.Graph)

The `bond matroid` or `cocylce matroid` of a graph `g` which is the dual of a cycle matroid, e.g, cographic.

See Section 2.3 of Oxl11 (@cite).

# Examples
To construct the bond or cocycle matroid of the complete graph of 4 vertices write:
```jldoctest
julia> g = Oscar.Graphs.complete_graph(4);

julia> M = bond_matroid(g)
Matroid of rank 3 on 6 elements
```

or equivalently

```jldoctest
julia> g = Oscar.Graphs.complete_graph(4);

julia> M = cocycle_matroid(g)
Matroid of rank 3 on 6 elements
```
"""
bond_matroid(g::Oscar.Graphs.Graph) = dual_matroid(cycle_matroid(g))

@doc Markdown.doc"""
See bond_matroid
"""
cocycle_matroid(g::Oscar.Graphs.Graph) = bond_matroid(g::Oscar.Graphs.Graph)


@doc Markdown.doc"""
    dual_matroid(M::Matroid)

The `dual matroid` of a given matroid `M`.

See page 65 and Sectrion 2 in Oxl11 (@cite).

# Examples
To construct the dual of the fano matroid write:
```jldoctest
julia> M = dual_matroid(fano_matroid())
Matroid of rank 4 on 7 elements
```
"""
dual_matroid(M::Matroid) = Matroid(M.pm_matroid.DUAL,M.groundset,M.gs2num)

@doc Markdown.doc"""
    groundset(M::Matroid)

The ground set `E` of a matroid `M`.

To obtain the ground set of the fano matroid type:
# Example
```jldoctest
julia> groundset(fano_matroid())
7-element Vector{Int64}:
 1
 2
 3
 4
 5
 6
 7
```
"""
matroid_groundset(M::Matroid) = M.groundset

@doc Markdown.doc"""
    direct_sum(M::Matroid, N::Matroid)
The `direct sum` of the matroids `M` and `N`.
Optionally one can also pass a vector of matroids.

See Section 4.2 of Oxl11 (@cite).

To obtain the direct sum of the fano and a uniform matroid type:
# Example
```jldoctest
julia> direct_sum(fano_matroid(), uniform_matroid(2,4))
Matroid of rank 5 on 11 elements
```

To take the sum of three uniform matroids use:
# Example
```jldoctest
julia> matroids = Vector([uniform_matroid(2,4), uniform_matroid(1,3), uniform_matroid(3,4)]);

julia> M = direct_sum(matroids)
Matroid of rank 6 on 11 elements
```
"""
function direct_sum(M::Matroid, N::Matroid)
    gsN = N.groundset
    while any(in(M.groundset), gsN)
        gsN = Vector{Any}(copy(gsN))
        for i in 1:length(gsN)
            gsN[i] = string(gsN[i],'\'')
        end
    end
    new_gs2num = Dict{Any,IntegerUnion}()
    i = length(M.groundset)+1
    for elem in gsN
        new_gs2num[elem] = i
        i+=1
    end
    gs2num = merge(M.gs2num,new_gs2num)
    return Matroid(Polymake.matroid.direct_sum(M.pm_matroid,N.pm_matroid),[M.groundset;gsN],gs2num)
end

direct_sum(comp::Vector{Matroid}) = foldl(direct_sum, comp)

@doc Markdown.doc"""
    deletion(M, [S, e])

# Arguments
- `M::Matroid`: A matroid `M`.
- `S::Union{AbstractVector,Set}`: A subset `S` of the ground set of `M`.
- `e::ElementType`: An element `e` of the ground set of `M`.

The `deletion M\S` of an element `e` or a subset `S` of the ground set `E` of the matroid `M`.

See Section 3 of Oxl11 (@cite).

# Example
```jldoctest
julia> M = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j']);

julia> N = deletion(M,'i')
Matroid of rank 2 on 3 elements
```

# Example
```jldoctest
julia> M = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j']);

julia> N = deletion(M,['i','j'])
Matroid of rank 2 on 2 elements

julia> groundset(N)
2-element Vector{Any}:
 1
 2
```
"""
function deletion(M::Matroid,set::Union{AbstractVector, Set})
    sort_set = Vector(undef,length(M.groundset)-length(set))
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in M.groundset
        if length(findall(x->x==elem, set))==0
            sort_set[i]=elem
            gs2num[elem] = i
            i+=1
        end
    end
    pm_del = Polymake.matroid.deletion(M.pm_matroid, Set([M.gs2num[i]-1 for i in set]))
    return Matroid(pm_del, sort_set, gs2num)
end

deletion(M::Matroid,elem::ElementType) = deletion(M,Vector([elem]))

@doc Markdown.doc"""
    restriction(M, S)

# Arguments
- `M::Matroid`: A matroid `M`.
- `S::Union{AbstractVector,Set}`: A subset `S` of the ground set of `M`.

The `restriction M|S` on a subset `S` of the ground set `E` of the matroid `M`.

See Section 3 of Oxl11 (@cite).

# Example
```jldoctest
julia> M = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j']);

julia> N = restriction(M,[1,2])
Matroid of rank 2 on 2 elements

julia> groundset(N)
2-element Vector{Int64}:
 1
 2
```
"""
function restriction(M::Matroid, set::Union{AbstractVector, Set})
    deleted_elems = filter(x -> x ∉ set, M.groundset)
    return deletion(M, deleted_elems)
end

@doc Markdown.doc"""
    contraction(M, [S, e])

# Arguments
- `M::Matroid`: A matroid `M`.
- `S::Union{AbstractVector,Set}`: A subset `S` of the ground set of `M`.
- `e::ElementType`: An element `e` of the ground set of `M`.

The `contraction M/S` of an element or a subset `S` of the ground set `E` of the matroid `M`.

See Section 3 of Oxl11 (@cite).

# Example
```jldoctest
julia> M = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j']);

julia> N = contraction(M,'i')
Matroid of rank 1 on 3 elements
```

# Example
```jldoctest
julia> M = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j']);

julia> N = contraction(M,['i','j'])
Matroid of rank 1 on 2 elements

julia> groundset(N)
2-element Vector{Any}:
 1
 2
```
"""
function contraction(M::Matroid,set::Union{AbstractVector, Set})
    # We use that the contraction by set is the dual of the deletion of set in the dual matroid. 
    return dual_matroid(deletion(dual_matroid(M), set))
end

contraction(M::Matroid,elem::ElementType) = contraction(M,Vector([elem]))

@doc Markdown.doc"""
    minor(M::Matroid, set_del::Union{AbstractVector, Set}, set_cont::Union{AbstractVector, Set})

The 'minor M\S/T` of disjoint subsets  `S` and `T` of the ground set `E` of the matroid `M`.

See also ``contraction`` and ``deletion``. You can find more in Section 3 of Oxl11 (@cite).

# Example
```jldoctest
julia> M = fano_matroid();

julia> S = [1,2,3];

julia> T = [4];

julia>  N = minor(M,S,T) 
Matroid of rank 2 on 3 elements
```
"""
function minor(M::Matroid, set_del::Union{AbstractVector, Set}, set_cont::Union{AbstractVector, Set})
    if any(in(set_del), set_cont)
        error("The two sets are not disjoined, which is required")
    end
    return contraction(deletion(M, set_del), set_cont)
end

@doc Markdown.doc"""
    principal_extension(M::Matroid, F::Union{AbstractVector,Set}, e::ElementType)

The `principal extension M +_F e` of a matroid `M` where the element `e` is freely added to the flat `F`.

See Section 7.2 of Oxl11 (@cite).

# Example
To add `4` freely to the flat `{1,2}` of the uniform matroid U_{3,4} do
```jldoctest
julia> M = uniform_matroid(3,4);

julia> N = principal_extension(M,[1,2],5)
Matroid of rank 3 on 5 elements
```
"""
function principal_extension(M::Matroid, set::Union{AbstractVector,Set}, elem::ElementType)
    if elem in M.groundset
        error("The element you are about to add is already contained in the ground set")
    end
    gs2num = copy(M.gs2num)
    gs2num[elem] = length(M.groundset)
    return Matroid(Polymake.matroid.principal_extension(M.pm_matroid,Set(set)),[M.groundset;elem],gs2num)
end

@doc Markdown.doc"""
    free_extension(M::Matroid, e::ElementType)
The `free extension M +_E e` of a matroid `M` where the element `e`.

See ``principal_extension`` and Section 7.2 of Oxl11 (@cite).

# Example
To add `4` freely to the uniform matroid U_{3,4} do
```jldoctest
julia> M = uniform_matroid(3,4);

julia>  N = free_extension(M,5)
Matroid of rank 3 on 5 elements
```
"""
free_extension(M::Matroid, elem::ElementType) = principal_extension(M, M.groundset, elem)

@doc Markdown.doc"""
    series_extension(M::Matroid, f::ElementType, e::ElementType)

The `series extension` of a matroid `M` where the element `e` is added in series to `f`.

This is actually a coextension see also Section 7.2 of Oxl11 (@cite).

# Example
To add `e` in series to `1` in the uniform matroid U_{3,4} do
```jldoctest
julia> M = uniform_matroid(1,4);

julia> N = series_extension(M,1,'e')
Matroid of rank 2 on 5 elements

julia> cocircuits(N)[1]
2-element Vector{Any}:
 2
  'e': ASCII/Unicode U+0065 (category Ll: Letter, lowercase)
```
"""
function series_extension(M::Matroid, old::ElementType, new::ElementType)
	return dual_matroid(principal_extension(dual_matroid(M),[old], new))
end


@doc Markdown.doc"""
    parallel_extension(M::Matroid, f::ElementType, e::ElementType)

The `parallel extension M +_{f} e` of a matroid `M` where the element `e` is added parallel to `f`.

See Section 7.2 of Oxl11 (@cite).

# Example
To add `e` parallel to `1` in the uniform matroid U_{3,4} do
```jldoctest
julia> M = uniform_matroid(3,4);

julia> N = parallel_extension(M,1,'e')
Matroid of rank 3 on 5 elements

julia> circuits(N)[1]
2-element Vector{Any}:
 2
  'e': ASCII/Unicode U+0065 (category Ll: Letter, lowercase)
```
"""
function parallel_extension(M::Matroid, old::ElementType, new::ElementType)
	if !(old in M.groundset)
		error("The element ".old." is not in the ground set")
	end
	return principal_extension(M,closure(M,[old]), new)
end


"""
    uniform_matroid(r,n)

Construct the uniform matroid of rank `r` on the `n` elements `{1,...,n}`.
"""
uniform_matroid(r::IntegerUnion,n::IntegerUnion) = Matroid(Polymake.matroid.uniform_matroid(r,n))

"""
    fano_matroid()

Construct the fano_matroid.
"""
fano_matroid() = matroid_from_matrix_rows(matrix(GF(2),[[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]]))

"""
    non_fano_matroid()

Construct the non-fano matroid.
"""
non_fano_matroid() = matroid_from_matrix_rows(matrix(QQ,[[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]]))

"""
    non_pappus_matroid()

Construct the non-Pappus matroid.
"""
non_pappus_matroid() = Matroid(Polymake.matroid.non_pappus_matroid())

"""
    pappus_matroid()

Construct the non-Pappus matroid.
"""
pappus_matroid() = Matroid(Polymake.matroid.pappus_matroid())

"""
    vamos_matroid()

Construct the Vamos matroid.
"""
vamos_matroid() = Matroid(Polymake.matroid.vamos_matroid())

"""
    all_subsets_matroid(r)

Construct the all-subsets-matroid of rank r, a.k.a. the matroid underlying the resonance arrangement.
"""
function all_subsets_matroid(r::Int)
    M = []
    for i in 1:2^r-1
        M = vcat(M,digits(i, base=2, pad=r))
    end
    M = convert(Array{Int64, 2}, reshape(M, r, 2^r-1))
    return matroid_from_matrix_columns(matrix(QQ, M); check=false)
end

@doc Markdown.doc"""
    projective_plane(q::Int)

The projective plane of order `q`.
Note that this only works for prime numbers `q` for now.

See Section 6.1 of Oxl11 (@cite).

# Example
```jldoctest
julia> M = projective_plane(3)
Matroid of rank 3 on 13 elements

```
"""
function projective_plane(q::Int)
    if !isprime(q)
        error("Only works for prime q at the moment.")
    end
    return Matroid(Polymake.matroid.projective_plane(q))
end


@doc Markdown.doc"""
    projective_geometry(r::Int, q::Int)

The projective geometry of order `q` and rank `r`.
Note that this only works for prime numbers `q` for now.

See Section 6.1 of Oxl11 (@cite).
Warning: Unlike in the book of Oxley, `r` is the actual rank of the matroid.

# Example
```jldoctest
julia> M = projective_geometry(3, 3)
Matroid of rank 3 on 13 elements

```
"""
function projective_geometry(r::Int, q::Int)
    if !isprime(q)
        error("q is not a prime.")
    end
    if r<3
        error("The rank should be at least 3")
    elseif r==3
        return projective_plane(q)
    end
    M=[]
    n=Int((q^r-1)/(q-1))
    for i in 1:(q^r-1)
        new_column = digits(i, base=q, pad=r)
        if new_column[findfirst(k->k!=0, new_column)]==1
            M = vcat(M, new_column)
        end
    end
    M = convert(Array{Int64,2},reshape(M, r, n))
    return matroid_from_matrix_columns(matrix(GF(q), M); check=false)
end

@doc Markdown.doc"""
    affine_geometry(r::Int, q::Int)

The affine geometry of order `q` and rank `r`.
Note that this only works for prime numbers `q` for now.

See Section 6.1 of Oxl11 (@cite).
Warning: Unlike in the book of Oxley, `r` is the actual rank of the matroid.

# Example
```jldoctest
julia> M = affine_geometry(3, 3)
Matroid of rank 3 on 9 elements
```
"""
function affine_geometry(r::Int, q::Int)
    if !isprime(q)
        error("q is not a prime.")
    end
    if r<3
        error("The rank should be at least 3")
    end

    PG = projective_geometry(r, q; check=false)
    return restriction(PG, Vector(2:(length(PG.groundset)-q)))
end
