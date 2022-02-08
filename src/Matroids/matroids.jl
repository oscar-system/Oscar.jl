export
	Matroid, groundset,
	matroid_from_bases, matroid_from_circuits,
	matroid_from_matrix, matroid_from_reduced_matrix,
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
    gs2num::Dict{Any,Int}# dictionary to map the groundset to the integers from 1 to its size
end

pm_object(M::Matroid) = M.pm_matroid

@doc Markdown.doc"""
A matroid with bases `B` on the ground set ``{1,2,...`n`}`` (which can be the empty set) where `n` is a non-negative integer and `B` is a collection of subsets of the ground set satisfying a exchange property.

# Examples
To construct a rank two matroid with five bases on the four elements 1,2,3,4 you can write:
```jldoctest
julia> B = [[1,2],[1,3],[1,4],[2,3],[2,4]];

julia> M = matroid_from_bases(B,4);
A matroid on 4 elements with 5 bases.
```

# Examples
To construct the same matroid on the four elements 2,1,i,j you may write:
```jldoctest
julia> M = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[2,1,'i','j']);
```
"""
matroid_from_bases(bases::Union{AbstractVector{<:AbstractVector{<:Base.Integer}}, AbstractVector{<:AbstractSet{<:Base.Integer}}}, nelements::Int) = matroid_from_bases(bases,Vector(1:nelements))

function matroid_from_bases(bases::Union{AbstractVector{<:AbstractVector}, AbstractVector{<:AbstractSet}},groundset::AbstractVector, check_input::Bool=true)
	if( check_input && size(groundset)[1]!=length(Set(groundset)) )
		throw("Input is not a valid groundset of a matroid")
	end
	gs2num = Dict{Any,Int}()
	i = 1
	for elem in groundset
		gs2num[elem] = i
		i+=1
	end
	pm_bases = [[gs2num[i]-1 for i in B] for B in bases]
    	M = Polymake.matroid.Matroid(BASES=pm_bases,N_ELEMENTS=size(groundset)[1])
	if(check_input && !Polymake.matroid.check_basis_exchange_axiom(M.BASES))
		throw("Input is not a collection of bases")
	end
    	Matroid(M,groundset,gs2num)
end

@doc Markdown.doc"""
A matroid with circuits `C` on the ground set ``{1,2,...`n`}`` (which can be the empty set) where `n` is a non-negative integer and `C` is a collection of subsets of the ground set.

# Examples
To construct a rank two matroid with five bases on the four elements 1,2,3,4 you can write:
```jldoctest
julia> B = [[1,2,3],[1,2,4],[3,4]];

julia> M = matroid_from_circuits(B,4);
A matroid on 4 elements with 5 bases.
```
# Examples
To construct the same matroid on the four elements 2,1,i,j you may write:
```jldoctest
julia> M = matroid_from_bases([[1,2,'i'],[1,2,'j'],['i','j']],[2,1,'i','j']);
```
"""
matroid_from_circuits(circuits::Union{AbstractVector{<:AbstractVector{<:Base.Integer}}, AbstractVector{<:AbstractSet{<:Base.Integer}}}, nelements::Int) = matroid_from_circuits(circuits,Vector(1:nelements))

function matroid_from_circuits(circuits::Union{AbstractVector{<:AbstractVector}, AbstractVector{<:AbstractSet}},groundset::AbstractVector, check_input::Bool=true)
	if( check_input && size(groundset)[1]!=length(Set(groundset)) )
		throw("Input is not a valid groundset of a matroid")
	end
	gs2num = Dict{Any,Int}()
	i = 1
	for elem in groundset
		gs2num[elem] = i
		i+=1
	end
	pm_circuits = [[gs2num[i]-1 for i in C] for C in circuits]
    	M = Polymake.matroid.Matroid(CIRCUITS=pm_circuits,N_ELEMENTS=size(groundset)[1])
	#TODO check_circuit_exchange_axiom
	#if(check_input && !Polymake.matroid.check_circuit_exchange_axiom(M.CIRCUITS))
	#	throw("Input is not a collection of circuits")
	#end
    	Matroid(M,groundset,gs2num)
end

@doc Markdown.doc"""
A matroid represented by the column vectors of a matrix.

# Examples
To construct the linear matroid of the matrix `A` over the field with two elements write:
```jldoctest
julia> A = [[1,0],[0,1],[1,1],[1,1]]

julia> M = matroid_from_matrix(matrix(GF(2),A))
```

or to obtain the same result
```jldoctest
julia> A = [[1,1],[1,1]]

julia> M = matroid_from_reduced_matrix(matrix(GF(2),A))
```
"""
matroid_from_reduced_matrix(A::MatrixElem) = matroid_from_matrix(hcat(identity_matrix(A.base_ring,nrows(A)),A))

function matroid_from_matrix(A::MatrixElem)
	rk = rank(A)
	bases = Vector{Vector{Int}}()
	for set in Oscar.Hecke.subsets(Vector(1:ncols(A)),rk)
		if(rank(A[:,set])==rk)
			push!(bases, set);
		end
	end
	matroid_from_bases(bases,ncols(A))
end

@doc Markdown.doc"""
The cycle matroid of a graph.

# Examples
To construct the cycle matroid of the complete graph of 4 vertices write:
```jldoctest
julia> g = Oscar.Graphs.complete_graph(4)

julia> M = cycle_matroid(g)
```
"""
function cycle_matroid(g::Oscar.Graphs.Graph)
	pm_Graph = Polymake.graph.Graph(ADJACENCY=g.pm_graph)
	M = Polymake.matroid.matroid_from_graph(pm_Graph)
	n = Oscar.Graphs.ne(g)
	gs2num = Dict{Any,Int}()
	i = 1
	for elem in 1:n
		gs2num[elem] = i
		i+=1
	end
	Matroid(M,1:n,gs2num)
end

@doc Markdown.doc"""
The bond or cocylce matroid of a graph which is the dual of a cycle matroid, e.g, cographic.

# Examples
To construct the bond or cocycle matroid of the complete graph of 4 vertices write:
```jldoctest
julia> g = Oscar.Graphs.complete_graph(4)

julia> M = bond_matroid(g)
```

or
```jldoctest
julia> M = cocycle_matroid(g)
```
"""
bond_matroid(g::Oscar.Graphs.Graph) = dual_matroid(cycle_matroid(g))
cocycle_matroid(g::Oscar.Graphs.Graph) = bond_matroid(g::Oscar.Graphs.Graph)


@doc Markdown.doc"""
The dual matroid of a given matroid.

# Examples
To construct the dual of the fano matroid write:
```jldoctest
julia> M = dual_matroid(fano_matroid())
```
"""
dual_matroid(M::Matroid) = Matroid(M.pm_matroid.DUAL,M.groundset,M.gs2num)

@doc Markdown.doc"""
The ground set of a matroid.

To obtain the ground set of the fano matroid type:
# Example
```jldoctest
julia> ground_set(fano())
```
"""
groundset(M::Matroid) = M.groundset

@doc Markdown.doc"""
The direct sum of matroids.

To obtain the direct sum of the fano and non-fano matroid type:
# Example
```jldoctest
julia> direct_sum(fano_matroid(),non_fano_matroid())
```

or to take the sum of more matroids use:
# Example
```jldoctest
julia> matroids = Vector(uniform_matroid(2,4), uniform_matroid(1,3), uniform_matroid(3,4))

julia> direct_sum(matroids)
```
"""
function direct_sum(M::Matroid, N::Matroid)
	gsN = N.groundset
	while(size( intersect(M.groundset,gsN) )[1]>0)
		gsN = Vector{Any}(copy(gsN))
		for i in 1:size(gsN)[1]
			gsN[i] = string(gsN[i],'\'')
		end
	end
	new_gs2num = Dict{Any,Int}()
	i = size(M.groundset)[1]+1
	for elem in gsN
		new_gs2num[elem] = i
		i+=1
	end
	gs2num = merge(M.gs2num,new_gs2num)
	Matroid(Polymake.matroid.direct_sum(M.pm_matroid,N.pm_matroid),[M.groundset;gsN],gs2num)
end

function direct_sum(comp::Vector{Matroid})
	M = comp[1]
	for i in 2:size(comp)[1]
		M = direct_sum(M,comp[i])
	end
	return M
end

@doc Markdown.doc"""
The deletion of an element or a subset of the ground set.

# Example
```jldoctest
julia> 
```
"""
function deletion(M::Matroid,set::Union{AbstractVector, Set})
	sort_set = Vector(undef,size(M.groundset)[1]-size(set)[1])
	gs2num = Dict{Any,Int}()
	i = 1
	for elem in M.groundset
		if(size(findall(x->x==elem, set))[1]==0)
			sort_set[i]=elem
			gs2num[elem] = i
			i+=1
		end
	end
	pm_del = Polymake.matroid.deletion(M.pm_matroid, Set([M.gs2num[i]-1 for i in set]))
	Matroid(pm_del, sort_set, gs2num)
end

deletion(M::Matroid,elem::Union{Int,Char,String}) = deletion(M,Vector([elem]))

function restriction(M::Matroid,set::Union{AbstractVector, Set})
	sort_set = copy(set)
	gs2num = Dict{Any,Int}()
	i = 1
	for elem in M.groundset
		if(size(findall(x->x==elem, set))[1]>0)
			sort_set[i]=elem
			gs2num[elem] = i
			i+=1
		end
	end
	pm_complement = setdiff(Set(0:size(M.groundset)[1]-1),Set([M.gs2num[i]-1 for i in set]))
	pm_rest = Polymake.matroid.deletion(M.pm_matroid, pm_complement)
	Matroid(pm_rest, sort_set, gs2num)
end

function contraction(M::Matroid,set::Union{AbstractVector, Set})
	sort_set = Vector(undef,size(M.groundset)[1]-size(set)[1])
	gs2num = Dict{Any,Int}()
	i = 1
	for elem in M.groundset
		if(size(findall(x->x==elem, set))[1]==0)
			sort_set[i]=elem
			gs2num[elem] = i
			i+=1
		end
	end
	pm_contr = Polymake.matroid.contraction(M.pm_matroid, Set([M.gs2num[i]-1 for i in set]))
	Matroid(pm_contr, sort_set, gs2num)
end

contraction(M::Matroid,elem::Union{Int,Char,String}) = contraction(M,Vector([elem]))


function minor(M::Matroid,set_del::Union{AbstractVector, Set},set_cont::Union{AbstractVector, Set})
	if( length(intersect(set_del,set_cont)>0) )
		throw("The two sets are not disjoined, which is required")
	end
	return contraction(deletion(M, set_del), set_contr)
end

function principal_extension(M::Matroid, set::Union{AbstractVector,Set}, elem::Union{Int,Char,String})
	if(issubset([elem],M.groundset))
		throw("The element you are about to add is already contained in the ground set")
	end
	gs2num = copy(M.gs2num)
	gs2num[elem] = length(M.groundset)
	Matroid(Polymake.matroid.principal_extension(M.pm_matroid,Set(set)),[M.groundset;elem],gs2num)
end


function uniform_matroid(r::Int,n::Int)
	gs2num = Dict{Any,Int}()
	i = 1
	for elem in 1:n
		gs2num[elem] = i
		i+=1
	end
	Matroid(Polymake.matroid.uniform_matroid(r,n),1:n,gs2num)
end

fano_matroid() = Matroid(Polymake.matroid.fano_matroid(),0:6, Dict{Any,Int}(1=>0, 2=>1, 3=>2, 4=>3, 5=>4, 6=>5, 7=>6))

non_fano_matroid() = Matroid(Polymake.matroid.non_fano_matroid(),0:6, Dict{Any,Int}(1=>0, 2=>1, 3=>2, 4=>3, 5=>4, 6=>5, 7=>6))

