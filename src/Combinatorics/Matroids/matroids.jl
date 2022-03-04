export
    Matroid, groundset,
    matroid_from_bases, matroid_from_circuits,
    matroid_from_matrix_columns, matroid_from_matrix_rows,
    cycle_matroid, bond_matroid, cocycle_matroid,
    dual_matroid, direct_sum, restriction, deletion, contraction, minor,
    principal_extension,
    uniform_matroid, fano_matroid, non_fano_matroid

################################################################################
##  Constructing
################################################################################

struct Matroid
    pm_matroid::Polymake.BigObject
    groundset::Vector # groundset of the matroid 
    gs2num::Dict{Any, IntegerUnion}# dictionary to map the groundset to the integers from 1 to its size
end

pm_object(M::Matroid) = M.pm_matroid

function Base.show(io::IO, M::Matroid)
    r = rank(M)
    n = length(M.groundset)
    print(io, "Matroid of rank $(r) on $(n) elements")
end


@doc Markdown.doc"""
Construct a `matroid` with bases `B` on the ground set `E` (which can be the empty set).
The set `B` is a non-empty collection of subsets of the ground set `E` satisfying an exchange property,
and the default value for `E` is the set `{1,..n}` for a non-negative value `n`.

See Section 1.2 of Oxl11 (@cite)

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
julia> M = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j']);
Matroid of rank 2 on 4 elements
```
"""
matroid_from_bases(bases::Union{AbstractVector{<:AbstractVector{<:IntegerUnion}}, AbstractVector{<:AbstractSet{<:IntegerUnion}}}, nelements::IntegerUnion) = matroid_from_bases(bases,Vector(1:nelements))

function matroid_from_bases(bases::Union{AbstractVector{<:AbstractVector}, AbstractVector{<:AbstractSet}},groundset::AbstractVector; check::Bool=true)
    if check && size(groundset)[1]!=length(Set(groundset))
        error("Input is not a valid groundset of a matroid")
    end
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in groundset
        gs2num[elem] = i
        i+=1
    end
    pm_bases = [[gs2num[i]-1 for i in B] for B in bases]
    M = Polymake.matroid.Matroid(BASES=pm_bases,N_ELEMENTS=size(groundset)[1])
    if check && !Polymake.matroid.check_basis_exchange_axiom(M.BASES)
        error("Input is not a collection of bases")
    end
    return Matroid(M,groundset,gs2num)
end

@doc Markdown.doc"""
A matroid with circuits `C` on the ground set `E` (which can be the empty set).
The set `C` is a collection of subsets of the ground set `E` satisfying an exchange property,
and the default value for `E` is the set `{1,..n}` for a non-negative value `n`. 

See Section 1.1 of Oxl11 (@cite)

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

julia> M = matroid_from_circuits(C,4)
Matroid of rank 2 on 4 elements
```
"""
matroid_from_circuits(circuits::Union{AbstractVector{<:AbstractVector{<:IntegerUnion}}, AbstractVector{<:AbstractSet{<:IntegerUnion}}}, nelements::IntegerUnion) = matroid_from_circuits(circuits,Vector(1:nelements))

function matroid_from_circuits(circuits::Union{AbstractVector{<:AbstractVector}, AbstractVector{<:AbstractSet}},groundset::AbstractVector; check::Bool=true)
    if check && size(groundset)[1]!=length(Set(groundset))
        error("Input is not a valid groundset of a matroid")
    end
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in groundset
        gs2num[elem] = i
        i+=1
    end
    pm_circuits = [[gs2num[i]-1 for i in C] for C in circuits]
    M = Polymake.matroid.Matroid(CIRCUITS=pm_circuits,N_ELEMENTS=size(groundset)[1])
    #TODO check_circuit_exchange_axiom (requires an update of polymake)
    #if check && !Polymake.matroid.check_circuit_exchange_axiom(M.CIRCUITS)
    #   error("Input is not a collection of circuits")
    #end
    return Matroid(M,groundset,gs2num)
end

@doc Markdown.doc"""
A matroid represented by the column vectors of a matrix.

See Section 1.1 of Oxl11 (@cite)

# Examples
To construct the vector matroid (a.k.a linear matroid) of the matrix `A` over the field with two elements write:
```jldoctest
julia> A = matrix(GF(2),[[1,0,1,1],[0,1,1,1]]);

julia> M = matroid_from_matrix_columns(A)
Matroid of rank 2 on 4 elements
```
"""
function matroid_from_matrix_columns(A::MatrixElem)
    rk = rank(A)
    bases = Vector{Vector{IntegerUnion}}()
    for set in Oscar.Hecke.subsets(Vector(1:ncols(A)),rk)
        if rank(A[:,set])==rk
            push!(bases, set);
        end
    end
    return matroid_from_bases(bases,ncols(A))
end

@doc Markdown.doc"""
A matroid represented by the row vectors of a matrix.

See Section 1.1 of Oxl11 (@cite)

# Examples
To construct the linear matroid of the rows of the matrix `A` over the field with two elements write:
```jldoctest
julia> A = matrix(GF(2),[[1,0],[0,1],[1,1],[1,1]]);

julia> M = matroid_from_matrix_rows(A)
Matroid of rank 2 on 4 elements
```
"""
matroid_from_matrix_rows(A::MatrixElem) = matroid_from_matrix_columns(transpose(A))

@doc Markdown.doc"""
The cycle matroid of a graph.

See Section 1.1 of Oxl11 (@cite)

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
The `bond matroid` or `cocylce matroid` of a graph which is the dual of a cycle matroid, e.g, cographic.

See Section 2.3 of Oxl11 (@cite)

# Examples
To construct the bond or cocycle matroid of the complete graph of 4 vertices write:
```jldoctest
julia> g = Oscar.Graphs.complete_graph(4);

julia> M = bond_matroid(g)
Matroid of rank 3 on 6 elements
```

or equivalently

```jldoctest
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
The `dual matroid` of a given matroid.

See page 65 and Sectrion 2 in Oxl11 (@cite)

# Examples
To construct the dual of the fano matroid write:
```jldoctest
julia> M = dual_matroid(fano_matroid())
Matroid of rank 4 on 7 elements
```
"""
dual_matroid(M::Matroid) = Matroid(M.pm_matroid.DUAL,M.groundset,M.gs2num)

@doc Markdown.doc"""
The ground set `E` of a matroid.

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
groundset(M::Matroid) = M.groundset

@doc Markdown.doc"""
The `direct sum` of matroids.

See Section 4.2 of Oxl11 (@cite)

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
        for i in 1:size(gsN)[1]
            gsN[i] = string(gsN[i],'\'')
        end
    end
    new_gs2num = Dict{Any,IntegerUnion}()
    i = size(M.groundset)[1]+1
    for elem in gsN
        new_gs2num[elem] = i
        i+=1
    end
    gs2num = merge(M.gs2num,new_gs2num)
    return Matroid(Polymake.matroid.direct_sum(M.pm_matroid,N.pm_matroid),[M.groundset;gsN],gs2num)
end

direct_sum(comp::Vector{Matroid}) = foldl(direct_sum, comp)

@doc Markdown.doc"""
The `deletion M\S` of an element or a subset `S` of the ground set `E` of the matroid `M`.

See Section 3 of Oxl11 (@cite)

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
    sort_set = Vector(undef,size(M.groundset)[1]-size(set)[1])
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in M.groundset
        if size(findall(x->x==elem, set))[1]==0
            sort_set[i]=elem
            gs2num[elem] = i
            i+=1
        end
    end
    pm_del = Polymake.matroid.deletion(M.pm_matroid, Set([M.gs2num[i]-1 for i in set]))
    return Matroid(pm_del, sort_set, gs2num)
end

deletion(M::Matroid,elem::Union{IntegerUnion,Char,String}) = deletion(M,Vector([elem]))

@doc Markdown.doc"""
The `restriction M|S` on a subset `S` of the ground set `E` of the matroid `M`.

See Section 3 of Oxl11 (@cite)

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
function restriction(M::Matroid,set::Union{AbstractVector, Set})
    sort_set = copy(set)
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in M.groundset
        if size(findall(x->x==elem, set))[1]>0
            sort_set[i]=elem
            gs2num[elem] = i
            i+=1
        end
    end
    pm_complement = setdiff(Set(0:size(M.groundset)[1]-1),Set([M.gs2num[i]-1 for i in set]))
    pm_rest = Polymake.matroid.deletion(M.pm_matroid, pm_complement)
    return Matroid(pm_rest, sort_set, gs2num)
end

@doc Markdown.doc"""
The `contraction M/S` of an element or a subset `S` of the ground set `E` of the matroid `M`.

See Section 3 of Oxl11 (@cite)

# Example
```jldoctest
julia> M = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j']);

julia> N = contraction(M,'i')
Matroid of rank 1 on 3 elements
```

# Example
```jldoctest
julia> M = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j']);

julia> N = deletion(M,['i','j'])
Matroid of rank 0 on 2 elements

julia> groundset(N)
2-element Vector{Any}:
 1
 2
```
"""
function contraction(M::Matroid,set::Union{AbstractVector, Set})
    sort_set = Vector(undef,size(M.groundset)[1]-size(set)[1])
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in M.groundset
        if size(findall(x->x==elem, set))[1]==0
            sort_set[i]=elem
            gs2num[elem] = i
            i+=1
        end
    end
    pm_contr = Polymake.matroid.contraction(M.pm_matroid, Set([M.gs2num[i]-1 for i in set]))
    return Matroid(pm_contr, sort_set, gs2num)
end

contraction(M::Matroid,elem::Union{IntegerUnion,Char,String}) = contraction(M,Vector([elem]))

@doc Markdown.doc"""
The 'minor M\S/T` of disjoint subsets  `S` and `T` of the ground set `E` of the matroid `M`.

See also ``contraction`` and ``deletion``. You can find more in Section 3 of Oxl11 (@cite)

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
The `principal extension M +_F e` of a matroid `M` where the element `e` is freely added to the flat `F`.

See Section 7.2 of Oxl11 (@cite)

# Example
To add `4` freely to the flat `{1,2}` of the uniform matroid U_{2,3} do
```jldoctest
julia> M = uniform_matroid(3,4);

julia>  N = principal_extension(M,[1,2],5)
Matroid of rank 3 on 5 elements
```
"""
function principal_extension(M::Matroid, set::Union{AbstractVector,Set}, elem::Union{IntegerUnion,Char,String})
    if issubset([elem],M.groundset)
        error("The element you are about to add is already contained in the ground set")
    end
    gs2num = copy(M.gs2num)
    gs2num[elem] = length(M.groundset)
    return Matroid(Polymake.matroid.principal_extension(M.pm_matroid,Set(set)),[M.groundset;elem],gs2num)
end


"""
    uniform_matroid(r,n)

Construct the uniform matroid of rank `r` on the `n` elements `{1,...,n}`.
"""
function uniform_matroid(r::IntegerUnion,n::IntegerUnion)
    gs2num = Dict{Any,IntegerUnion}()
    i = 1
    for elem in 1:n
        gs2num[elem] = i
        i+=1
    end
    return Matroid(Polymake.matroid.uniform_matroid(r,n),1:n,gs2num)
end

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
