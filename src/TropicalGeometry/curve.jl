###
# Tropical curves in Oscar
# ========================
###

using RecipesBase
###
# 1. Definition
# -------------
# M = typeof(min) or typeof(max):
#   min or max convention, affecting initial ideals
# EMB = true or false:
#   embedded or abstract tropical curves
#   embedded tropical variety = graph embedded in euclidean space with weighted edges and vertices
#   abstract tropical variety = graph with enumerated vertices with weighted edges and vertices
###

@attributes mutable struct TropicalCurve{M,EMB} <: TropicalVarietySupertype{M,EMB}
    polyhedralComplex::PolyhedralComplex
    function TropicalCurve{M,EMB}(Sigma::PolyhedralComplex) where {M,EMB}
        if EMB
            if dim(Sigma)!=1
                error("TropicalCurve: input polyhedral complex not one-dimensional")
            end
        end
        return new{M,EMB}(Sigma)
    end
end
export TropicalCurve

function pm_object(T::TropicalCurve{M, EMB}) where {M, EMB}
    if has_attribute(T,:polymake_bigobject)
        return get_attribute(T,:polymake_bigobject)
    end
    error("pm_object(T::TropicalCurve): no polymake bigobject attributed")
end


###
# 2. Basic constructors
# ---------------------
###

@doc Markdown.doc"""
    TropicalCurve(PC::PolyhedralComplex)

Construct a tropical curve from a polyhedral complex.
If the curve is embedded, vertices must are points in $\mathbb R^n$.
If the curve is abstract, the polyhedral complex is empty, vertices must be 1, ..., n,
and the graph is given as attribute.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4]])
3×4 IncidenceMatrix
[1, 2]
[1, 3]
[1, 4]


julia> VR = [0 0; 1 0; -1 0; 0 1]
4×2 Matrix{Int64}:
  0  0
  1  0
 -1  0
  0  1

julia> PC = PolyhedralComplex{fmpq}(IM, VR)
A polyhedral complex in ambient dimension 2

julia> TC = TropicalCurve(PC)
A min tropical curve in 2-dimensional Euclidean space

julia> abs_TC = TropicalCurve(IM)
An abstract min tropical curve
```
"""
function TropicalCurve(PC::PolyhedralComplex, M::Union{typeof(min),typeof(max)}=min)
   @assert dim(PC)==1 "The polyhedral complex is not of dimenion 1."
   return TropicalCurve{M, true}(PC)
end

function TropicalCurve(graph::IncidenceMatrix, M::Union{typeof(min),typeof(max)}=min)
    # Columns correspond to nodes
    # Rows correpons to edges
    empty = PolyhedralComplex(Polymake.fan.PolyhedralComplex())
    result = TropicalCurve{M, false}(empty)
    set_attribute!(result, :graph, graph)
    return result
end

@doc Markdown.doc"""
    graph(tc::TropicalCurve)

Return the graph of an abstract tropical curve `tc`.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve(IM)
An abstract min tropical curve

julia> graph(tc)
6×4 IncidenceMatrix
[1, 2]
[1, 3]
[1, 4]
[2, 3]
[2, 4]
[3, 4]
```
"""
function graph(tc::TropicalCurve)
    if !has_attribute(tc, :graph)
        throw(ArgumentError("No graph attached"))
    end
    return get_attribute(tc, :graph)
end

@doc Markdown.doc"""
    n_nodes(tc::TropicalCurve)

Return the number of nodes of an abstract tropical curve `tc`.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve(IM)
An abstract min tropical curve

julia> n_nodes(tc)
4
```
"""
function n_nodes(tc::TropicalCurve)
    G = graph(tc)
    return Polymake.ncols(G)
end


function Base.show(io::IO, tc::TropicalCurve{M,EMB}) where {M,EMB}
    if EMB
        print(io, "A "*string(M)*" tropical curve in $(ambient_dim(tc))-dimensional Euclidean space")
    else
        print(io, "An abstract "*string(M)*" tropical curve")
    end
end

@doc Markdown.doc"""
    DivisorOnTropicalCurve(tc::TropicalCurve, coeffs::Vector{Int})

Construct a divisor with coefficients `coeffs` on an abstract tropical curve `tc`.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve(IM)
An abstract min tropical curve

julia> coeffs = [0, 1, 1, 1];

julia> dtc = DivisorOnTropicalCurve(tc,coeffs)
DivisorOnTropicalCurve{min, false}(An abstract min tropical curve, [0, 1, 1, 1])
```
"""
struct DivisorOnTropicalCurve{M, EMB}
    base_curve::TropicalCurve{M, EMB}
    coefficients::Vector{Int}
    function DivisorOnTropicalCurve(tc::TropicalCurve{M, EMB}, coeffs::Vector{Int}) where {M,EMB}
        if EMB
            error("Not implemented yet")
        else
            if n_nodes(tc) != length(coeffs)
                throw(ArgumentError("Wrong number coefficients"))
            end
            return new{M, EMB}(tc, coeffs)
        end
    end
end

base_curve(dtc::DivisorOnTropicalCurve) = dtc.base_curve

###
# 3.Basic properties
#
@doc Markdown.doc"""
    coefficients(dtc::DivisorOnTropicalCurve)

Construct a divisor `dtc` with coefficients `coeffs` on an abstract tropical curve.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve(IM)
An abstract min tropical curve

julia> coeffs = [0, 1, 1, 1];

julia> dtc = DivisorOnTropicalCurve(tc,coeffs)
DivisorOnTropicalCurve{min, false}(An abstract min tropical curve, [0, 1, 1, 1])

julia> coefficients(dtc)
4-element Vector{Int64}:
 0
 1
 1
 1
```
"""
coefficients(dtc::DivisorOnTropicalCurve) = dtc.coefficients

@doc Markdown.doc"""
   degree(dtc::DivisorOnTropicalCurve)

Compute the degree of  a divisor `dtc` on an abstract tropical curve.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve(IM)
An abstract min tropical curve

julia> coeffs = [0, 1, 1, 1];

julia> dtc = DivisorOnTropicalCurve(tc,coeffs)
DivisorOnTropicalCurve{min, false}(An abstract min tropical curve, [0, 1, 1, 1])

julia> degree(dtc)
3
```
"""
degree(dtc::DivisorOnTropicalCurve) = sum(coefficients(dtc))

@doc Markdown.doc"""
    is_effective(dtc::DivisorOnTropicalCurve)

Check whether a divisor `dtc` on an abstract tropical curve is effective.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve(IM)
An abstract min tropical curve

julia> coeffs = [0, 1, 1, 1];

julia> dtc = DivisorOnTropicalCurve(tc,coeffs)
DivisorOnTropicalCurve{min, false}(An abstract min tropical curve, [0, 1, 1, 1])

julia> is_effective(dtc)
true
```
"""
is_effective(dtc::DivisorOnTropicalCurve) = all(e -> e>=0, coefficients(dtc))


@doc Markdown.doc"""
   chip_firing_move(dtc::DivisorOnTropicalCurve, position::Int)

Given a divisor `dtc` and vertex labelled `position`, compute the linearly equivalent divisor obtained by a chip firing move from the given vertex `position`.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve(IM)
An abstract min tropical curve

julia> coeffs = [0, 1, 1, 1];

julia> dtc = DivisorOnTropicalCurve(tc,coeffs)
DivisorOnTropicalCurve{min, false}(An abstract min tropical curve, [0, 1, 1, 1])

julia> chip_firing_move(dtc,1)
DivisorOnTropicalCurve{min, false}(An abstract min tropical curve, [-3, 2, 2, 2])
```
"""
function chip_firing_move(dtc::DivisorOnTropicalCurve, position::Int)
    G = graph(base_curve(dtc))
    newcoeffs = Vector{Int}(coefficients(dtc))
    for i in 1:Polymake.nrows(G)
        row = Polymake.row(G, i)
        if position in row
            newcoeffs[position] -= 1
            for i in row
                if i != position
                    newcoeffs[i] += 1
                end
            end
        end
    end
    return DivisorOnTropicalCurve(base_curve(dtc), newcoeffs)
end

### The function computes the outdegree of a vertex v  with respect to a given subset W  of vertices.
### This is the number of vertices not in W adjecent to v. 1,
function outdegree(tc::TropicalCurve, W::Set{Int}, v::Int)
    G = graph(tc)
    m = Polymake.nrows(G) #number of edges of tc
    @assert v in W "Vertex number $v not in $W"
    deg = 0 #outdeg
    for i in 1:m
        row = Polymake.row(G,i)
        if v in row
            for j in row
                if !(j in W)
                    deg = deg +1
                end
            end
        end
    end
    return deg
end

@doc Markdown.doc"""
   v_reduced(dtc::DivisorOnTropicalCurve, vertex::Int)

Given a divisor `dtc` and vertex labelled `vertex`, compute the unique divisor reduced with repspect to `vertex`
as defined in [BN07](@cite).
The divisor `dtc` must have positive coefficients apart from `vertex`.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve(IM)
An abstract min tropical curve

julia> coeffs = [0, 1, 1, 1];

julia> dtc = DivisorOnTropicalCurve(tc,coeffs)
DivisorOnTropicalCurve{min, false}(An abstract min tropical curve, [0, 1, 1, 1])

julia> v_reduced(dtc,1)
DivisorOnTropicalCurve{min, false}(An abstract min tropical curve, [3, 0, 0, 0])
```
"""
function v_reduced(dtc::DivisorOnTropicalCurve, vertex::Int)
    tc = base_curve(dtc)
    G = graph(base_curve(dtc))
    n = Polymake.ncols(G)
    newcoeff = Vector{Int}(coefficients(dtc))
    S = Set{Int}(1:n)
    W = setdiff(S,vertex)
    @assert all(j -> newcoeff[j]>=0, W) "Divisor not effective outside the vertex number $vertex"
    while !(isempty(W))
	w0 = vertex
	if all(w -> newcoeff[w] >= outdegree(tc,W,w),W)
	    coeffdelta = zeros(Int, n)
	    delta = DivisorOnTropicalCurve(tc,coeffdelta)
	    for j in W
	        delta = chip_firing_move(delta,j)
	    end
	    newcoeff = newcoeff + coefficients(delta)
	else
	    for w in W
	        if newcoeff[w] < outdegree(tc,W,w)
		    w0 = w
		    break
	        end
	    end
	    W = setdiff(W,w0)
        end
    end
    reduced = DivisorOnTropicalCurve(tc, newcoeff)
    return reduced
end

@doc Markdown.doc"""
   is_linearly_equivalent(dtc1::DivisorOnTropicalCurve, dtc2::DivisorOnTropicalCurve)

Given two effective divisors `dtc1` and `dtc2` on the same tropical curve, check whether they are linearly equivalent.

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve(IM)
An abstract min tropical curve

julia> coeffs1 = [0, 1, 1, 1];

julia> dtc1 = DivisorOnTropicalCurve(tc,coeffs1)
DivisorOnTropicalCurve{min, false}(An abstract min tropical curve, [0, 1, 1, 1])

julia> coeffs2 = [3,0,0,0];

julia> dtc2 = DivisorOnTropicalCurve(tc,coeffs2)
DivisorOnTropicalCurve{min, false}(An abstract min tropical curve, [3, 0, 0, 0])

julia> is_linearly_equivalent(dtc1, dtc2)
true
```
"""
function is_linearly_equivalent(dtc1::DivisorOnTropicalCurve{M, EMB}, dtc2::DivisorOnTropicalCurve{M, EMB}) where {M, EMB}
    @assert is_effective(dtc1) "The divisor $dtc1 is not effective"
    @assert is_effective(dtc2) "The divisor $dtc2 is not effective"
    @assert base_curve(dtc1) === base_curve(dtc2) "The input curve needs to be the same for both divisors"
    v = 1
    reduced1 = v_reduced(dtc1,v)
    reduced2 = v_reduced(dtc2,v)
    coeff1 = coefficients(reduced1)
    coeff2 = coefficients(reduced2)
    return coeff1 == coeff2
end

export DivisorOnTropicalCurve,
    n_nodes,
    degree,
    coefficients,
    is_effective,
    graph,
    base_curve,
    chip_firing_move,
    v_reduced,
    is_linearly_equivalent,
    structure_tropical_jacobian,
    visualize

###
# 4. More properties
# -------------------
###

@doc Markdown.doc"""
    structure_tropical_jacobian(TC::TropicalCurve)

Compute the elementary divisors $n_i$ of the Laplacian matrix of the tropical curve `TC`.
The tropical Jacobian is then isomorphic to $\prod (Z/(n_i)Z)$.

# Examples
```jldoctest
julia> cg = Oscar.Graphs.complete_graph(5);

julia> IM1=IncidenceMatrix([[Oscar.Graphs.src(e), Oscar.Graphs.dst(e)] for e in Oscar.Graphs.edges(cg)])
10×5 IncidenceMatrix
[1, 2]
[1, 3]
[2, 3]
[1, 4]
[2, 4]
[3, 4]
[1, 5]
[2, 5]
[3, 5]
[4, 5]

julia> TC1 = TropicalCurve(IM1)
An abstract min tropical curve

julia> structure_tropical_jacobian(TC1)
(General) abelian group with relation matrix
[1 0 0 0; 0 5 0 0; 0 0 5 0; 0 0 0 5]

julia> cg2 = Oscar.Graphs.complete_graph(3);

julia> IM2=IncidenceMatrix([[Oscar.Graphs.src(e), Oscar.Graphs.dst(e)] for e in Oscar.Graphs.edges(cg2)])
3×3 IncidenceMatrix
[1, 2]
[1, 3]
[2, 3]

julia> TC2 = TropicalCurve(IM2)
An abstract min tropical curve

julia> structure_tropical_jacobian(TC2)
(General) abelian group with relation matrix
[1 0; 0 3]

julia> IM3 = IncidenceMatrix([[1,2],[2,3],[3,4],[4,5],[1,5]])
5×5 IncidenceMatrix
[1, 2]
[2, 3]
[3, 4]
[4, 5]
[1, 5]

julia> TC3=TropicalCurve(IM3)
An abstract min tropical curve

julia> G = structure_tropical_jacobian(TC3)
(General) abelian group with relation matrix
[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 5]
```
"""
function structure_tropical_jacobian(TC::TropicalCurve)
    gg=Graphs.Graph{Graphs.Undirected}(n_nodes(TC))
    IM = graph(TC)
    for i in 1:Polymake.nrows(IM)
        row = Vector{Int}(Polymake.row(IM,i))
        Graphs.add_edge!(gg, row[1],row[2])
    end
    lap = Polymake.graph.laplacian(Oscar.pm_object(gg))
    L = Polymake.@convert_to Matrix{Int} lap
    LL = matrix(ZZ, L)
    ED = elementary_divisors(LL)[1:nrows(LL)-1]
    G = abelian_group(ED)
    return G
end

@recipe function visualize(tc::TropicalCurve{M,EMB}) where {M,EMB}
    @assert EMB "Tropical curve is abstract."
    PC = tc.polyhedralComplex
    MaxPoly= maximal_polyhedra(PC)
    list_vertices = Vector{Complex{Float64}}()
    for P in MaxPoly
        V = vertices(P)
        R = rays(P)
        #V = [Vector{Rational{Int}}(x) for x in V]
        V = [Vector{Float64}(Vector{Rational}(Vector{Polymake.Rational}(x))) for x in V]
        #R = [Vector{Rational{Int}}(x) for x in R]
        R = [Vector{Float64}(Vector{Rational}(Vector{Polymake.Rational}(x))) for x in R]
        if length(V)==2
            append!(list_vertices,V[1][1]+V[1][2]*im,V[2][1]+V[2][2]*im,Inf+0*im)
        else
            B = V[1]+R[1]
            append!(list_vertices,V[1][1]+V[1][2]*im,B[1]+B[2]*im,Inf+0*im)
        end
    end
    return list_vertices #(list_vertices,legend=false,axis=false,xlabel="",ylabel="")
end
