export
	isomorphism, size_groundset,
	bases, circuits, hyperplanes, flats, cyclic_flats, closure, 
	rank, nullity, 
	fundamental_circuit, fundamental_cocircuit,
	spanning_sets, independent_sets,
	cobases, cocircuits, cohyperplanes, corank,
	is_clutter, get_blocker,
	is_regular, is_graphic, is_cographic,
	is_binary, is_ternary,
	n_connected_components, connected_components, is_connected,
	loops, coloops, is_loopless, is_coloopless, direct_sum_components,
	connectivity_function, is_vertical_k_separation,
	vertical_connectivity, girth, tutte_connectivity,
	principal_extension, free_extension, series_extension,
	tutte_polynomial, characteristic_polynomial, reduced_characteristic_polynomial,
	is_quotient

################################################################################
##  Properties and basic functions
################################################################################

function isomorphism(M::Matroid, gs::AbstractVector)
	if size(M.groundset)[1]!=size(gs)[1]
		throw("Sets of different size")
	end
	gs2num = Dict{Any,Int}()
	i = 1
	for elem in gs
		gs2num[elem] = i
		i+=1
	end
	return Matroid(M.pm_matroid,gs,gs2num)
end

size_groundset(M::Matroid) = length(M.groundset)[1]

@doc Markdown.doc"""
    bases(M::matroid)

Return the list of bases of the matroid `M`.

# Example
```jldoctest
julia> bases(fano())

```
"""
bases(M::Matroid) = [[M.groundset[i+1] for i in sort(collect(C))] for C in Vector{Set{Int}}(M.pm_matroid.BASES)]

@doc Markdown.doc"""
    circuits(M::matroid)

Return the list of circuits of the matroid `M`.

# Example
```jldoctest
julia> 

```
"""
circuits(M::Matroid) = [[M.groundset[i+1] for i in sort(collect(C))] for C in Vector{Set{Int}}(M.pm_matroid.CIRCUITS)]

hyperplanes(M::Matroid) = [[M.groundset[i+1] for i in sort(collect(C))] for C in Vector{Set{Int}}(M.pm_matroid.MATROID_HYPERPLANES)]


function flats(M::Matroid)
	pm_IncidenceMatrix = Polymake.@pm common.convert_to{IncidenceMatrix}(M.pm_matroid.LATTICE_OF_FLATS.FACES)
	indices = M.pm_matroid.LATTICE_OF_FLATS.TOP_NODE!=0 ? (1:Polymake.size(pm_IncidenceMatrix,1)) : reverse(1:Polymake.size(pm_IncidenceMatrix,1))
	return [[M.groundset[i] for i in flat] for flat in [Vector{Int}(Polymake.row(pm_IncidenceMatrix, i)) for i in indices]]
end

function cyclic_flats(M::Matroid)
	pm_IncidenceMatrix = Polymake.@pm common.convert_to{IncidenceMatrix}(M.pm_matroid.LATTICE_OF_CYCLIC_FLATS.FACES)
	indices = M.pm_matroid.LATTICE_OF_CYCLIC_FLATS.TOP_NODE==0 ? (1:Polymake.size(pm_IncidenceMatrix,1)) : reverse(1:Polymake.size(pm_IncidenceMatrix,1))
	return [[M.groundset[i] for i in flat] for flat in [Vector{Int}(Polymake.row(pm_IncidenceMatrix, i)) for i in indices]]
end

function closure(M::Matroid,set::Union{AbstractVector,Set})
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
julia> M = 
julia> groundset(M)
```
"""
rank(M::Matroid) = M.pm_matroid.RANK::Int


@doc Markdown.doc"""
    rank(M::matroid,set::Vector)

Return the rank of `set` in the matroid `M`.

# Example
```jldoctest
julia> 

```
"""
rank(M::Matroid, set::Union{AbstractVector,Set}) = Polymake.matroid.rank( M.pm_matroid, Set([M.gs2num[i]-1 for i in set]) )

nullity(M::Matroid, set::Union{AbstractVector,Set}) = size(set)[1]-Polymake.matroid.rank( M.pm_matroid, Set([M.gs2num[i]-1 for i in set]) )

function fundamental_circuit(M::Matroid, basis::AbstractVector, elem::Union{Int,Char,String})
	if !(basis in bases(M))
		throw("The set is not a basis of M")
	end
	if !(elem in M.groundset)
		throw("The element is not in the groundset")
	end
	if elem in basis
		throw("The element is in the basis, but has to be in the complement")
	end
	return circuits(restriction(M,[elem; basis]))[1]
end

fundamental_cocircuit(M::Matroid, basis::AbstractVector, elem::Union{Int,Char,String}) = fundamental_circuit(dual_matroid(M),basis,elem)

@doc Markdown.doc"""
    independet_sets(M::matroid)

Return the list of independent sets of the matroid `M`.

# Example
```jldoctest
julia> 

```
"""
function spanning_sets(M::Matroid)
	pm_bases = Vector{Set{Int}}(M.pm_matroid.BASES)
	n = size(M.groundset)[1]
	gs = M.groundset
	sets = Vector()
	for k in rank(M):n
		for set in Oscar.Hecke.subsets(Vector(0:n-1),k)
			for B in pm_bases
				if issubset(B,set)
					push!(sets, [Vector([gs[i+1] for i in set])])
					break
				end
			end
		end
	end
	return sets
end

function independent_sets(M::Matroid)
	pm_bases = Vector{Set{Int}}(M.pm_matroid.BASES)
	n = size(M.groundset)[1]
	gs = M.groundset
	sets = Vector(undef,1)
	for k in 1:rank(M)
		for set in Oscar.Hecke.subsets(Vector(0:n-1),k)
			for B in pm_bases
				if issubset(set,B)
					push!( sets, [Vector([gs[i+1] for i in set])])
					break
				end
			end
		end
	end
	sets[1]=[]
	return sets
end

cobases(M::Matroid) = bases( dual_matroid(M) )

cocircuits(M::Matroid) = circuits( dual_matroid(M) )

cohyperplanes(M::Matroid) = hyperplanes( dual_matroid(M) )

corank(M::Matroid, set::Vector) = size(set)[1]-rank(M, set) + rank(M, setdiff(M.groundset,set))

function is_clutter(sets::Union{AbstractVector{<:AbstractVector}, AbstractVector{<:AbstractSet}})
	for A in sets
		for B in sets
			if issubset(A,B) ⊻ issubset(B,A)
				return false
			end
		end
	end
	return true
end

function get_blocker(clutter::Union{AbstractVector{<:AbstractVector}, AbstractVector{<:AbstractSet}}, groundset::Vector)
	#TODO
end

is_regular(M::Matroid) = M.pm_matroid.REGULAR

is_graphic(M::Matroid) = true #TODO

is_cographic(M::Matroid) = is_graphic(dual_matroid(M))

is_binary(M::Matroid) = M.pm_matroid.BINARY

is_ternary(M::Matroid) = M.pm_matroid.TERNARY

n_connected_components(M::Matroid) = size(M.pm_matroid.CONNECTED_COMPONENTS)[1]

connected_components(M::Matroid) = [[M.groundset[i+1] for i in comp] for comp in M.pm_matroid.CONNECTED_COMPONENTS]
is_connected(M::Matroid) = M.pm_matroid.CONNECTED

loops(M::Matroid) = [M.groundset[i+1] for i in M.pm_matroid.LOOPS]

coloops(M::Matroid) = [M.groundset[i+1] for i in M.pm_matroid.DUAL.LOOPS]

is_loopless(M::Matroid) = length(M.pm_matroid.LOOPS)==0 ? true : false

is_coloopless(M::Matroid) = length(M.pm_matroid.DUAL.LOOPS)==0 ? true : false

function direct_sum_components(M::Matroid)
	res = Vector{Matroid}()
	for set in connected_components(M)
		res = [res; restriction(M,set)]
	end
	return res
end

connectivity_function(M::Matroid, set::Union{AbstractVector, Set}) = rank(M,set) + rank(M,setdiff( Set(M.groundset), set)) - rank(M)

is_vertical_k_separation(M::Matroid,k::Int, set::Union{AbstractVector,Set}) = k<=rank(M,set) && k<=rank(M,setdiff( Set(M.groundset), set)) && k> connectivity_function(M,set)

is_k_separation(M::Matroid,k::Int, set::Union{AbstractVector,Set}) = k<=size(set)[1] && k<=size(setdiff( Set(M.groundset), set))[1] && k> connectivity_function(M,set)

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

girth(M::Matroid, set::Vector=M.groundset) = minimum([inf; [issubset(C,set) ? size(C)[1] : inf for C in circuits(M)]])


function tutte_connectivity(M::Matroid)
	r = M.pm_matroid.RANK
	n = M.pm_matroid.N_ELEMENTS
	#if M is uniform, apply Cor. 8.6.3 otherwise Thm. 8.6.4
	if n>2r && M.pm_matroid.N_BASES==binomial(n,r)
		return n>=2r-2 ? r+1 : inf
	end
	return minimum([vertical_connectivity(M),girth(M)])
end

principal_extension()=0
free_extension()=0
series_extension()=0

function tutte_polynomial(M::Matroid)
	R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
	poly = M.pm_matroid.TUTTE_POLYNOMIAL
	exp = Polymake.monomials_as_matrix(poly)
	return R(Vector{Int}(Polymake.coefficients_as_vector(poly)),[[exp[i,1],exp[i,2]] for i in 1:size(exp)[1]])
end

function characteristic_polynomial(M::Matroid)
	R, q = PolynomialRing(ZZ, 'q')
	return (-1)^M.pm_matroid.RANK*tutte_polynomial(M)(1-q,0)
end

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

is_quotient() = 0
