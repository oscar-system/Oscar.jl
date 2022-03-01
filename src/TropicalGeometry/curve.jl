###
# Tropical curves in Oscar
# ========================
###



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

function pm_object(T::TropicalCurve)
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
    TropicalCurve{M, EMB}()

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


julia>  VR = [0 0; 1 0; -1 0; 0 1]
4×2 Matrix{Int64}:
  0  0
  1  0
 -1  0
  0  1

julia> PC = PolyhedralComplex(IM, vr)
A polyhedral complex in ambient dimension 2

julia> TC = TropicalCurve{min}(PC)
A tropical curve in 2-dimensional Euclidean space

julia> abs_TC = TropicalCurve{min}(IM)
An abstract tropical curve
```
"""
function TropicalCurve{M}(PC::PolyhedralComplex) where {M}
   @assert dim(PC)==1 "The polyhedral complex is not of dimenion 1."   
  return TropicalCurve{M, true}(PC)
end


function TropicalCurve{M}(graph::IncidenceMatrix) where {M}
    # Columns correspond to nodes
    # Rows correpons to edges
    empty = PolyhedralComplex(Polymake.fan.PolyhedralComplex())
    result = TropicalCurve{M, false}(empty)
    set_attribute!(result, :graph, graph)    
    return result
end

@doc Markdown.doc"""
    graph(tc::TropicalCurve{M, EMB})

Returns the graph of an abstract tropical curve `tc`. 
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve{min}(IM)
An abstract tropical curve

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
function graph(tc::TropicalCurve{M, EMB}) where {M, EMB}
    if !has_attribute(tc, :graph)
        throw(ArgumentError("No graph attached"))
    end
    return get_attribute(tc, :graph)
end

@doc Markdown.doc"""
    n_nodes(tc::TropicalCurve{M, EMB})      

Returns the number of nodes of an abstract tropical curve `tc`.
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve{min}(IM)
An abstract tropical curve

julia> n_nodes(tc)
4
```
"""
function n_nodes(tc::TropicalCurve{M, EMB}) where {M, EMB}
    G = graph(tc)
    return Polymake.ncols(G)
end


function Base.show(io::IO, tc::TropicalCurve{M, EMB}) where {M, EMB}
    if EMB
        print(io, "A tropical curve in $(ambient_dim(tc))-dimensional Euclidean space")
    else
        print(io, "An abstract tropical curve")
    end
end

@doc Markdown.doc"""
    DivisorOnTropicalCurve(tc::TropicalCurve{M, EMB}, coeffs::Vector{Int})       

Construct a divisor with coefficients `coeffs` on an abstract tropical curve `tc`.
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve{min}(IM)
An abstract tropical curve

julia> coeffs = [0, 1, 1, 1];

julia> dtc = DivisorOnTropicalCurve(tc,coeffs)
DivisorOnTropicalCurve{min, false}(An abstract tropical curve, [0, 1, 1, 1])
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

base_curve(dtc::DivisorOnTropicalCurve{M, EMB}) where {M, EMB} = dtc.base_curve

###
# 3.Basic properties
#
@doc Markdown.doc"""
    coefficients(dtc::DivisorOnTropicalCurve{M, EMB})            

Construct a divisor with coefficients `coeffs` on an abstract tropical curve `tc`.
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

julia> tc = TropicalCurve{min}(IM)
An abstract tropical curve

julia> coeffs = [0, 1, 1, 1];

julia> dtc = DivisorOnTropicalCurve(tc,coeffs)
DivisorOnTropicalCurve{min, false}(An abstract tropical curve, [0, 1, 1, 1])

julia> coefficients(dtc)
4-element Vector{Int64}:
 0
 1
 1
 1
```
"""
coefficients(dtc::DivisorOnTropicalCurve{M, EMB}) where {M, EMB} = dtc.coefficients

degree(dtc::DivisorOnTropicalCurve{M, EMB}) where {M, EMB} = sum(coefficients(dtc))

is_effective(dtc::DivisorOnTropicalCurve{M, EMB}) where {M, EMB} = all(e -> e>=0, coefficients(dtc))

function chip_firing_move(dtc::DivisorOnTropicalCurve{M, EMB}, position::Int) where {M, EMB}
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

function outdegree(tc::TropicalCurve{M,EMB}, W::Set{Int}, v::Int) where {M, EMB}
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

function v_reduced(dtc::DivisorOnTropicalCurve{M, EMB}, vertex::Int) where {M, EMB}
    tc = base_curve(dtc)
    G = graph(base_curve(dtc))
    n = Polymake.ncols(G)
    newcoeff = Vector{Int}(coefficients(dtc))
    S = Set{Int}(1:n)
    W = setdiff(S,vertex)
    @assert all(j -> newcoeff[j]>=0, W) "Divisor not effective outside the vertex number $vertex"
    while !(isempty(W)) 
	for w in W 
	    if newcoeff[w] < outdegree(tc,W,w)
	        W = setdiff(W,w)	    	
  	    else 
		coeffdelta = zeros(Int, n)
		delta = DivisorOnTropicalCurve(tc,coeffdelta)
	   	for j in W 
	            delta = chip_firing_move(delta,j)
	        end
	        newcoeff = newcoeff + coefficients(delta)
	    end
       end 
    end
    reduced = DivisorOnTropicalCurve(tc, newcoeff)
    return reduced   
end

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
    outdegree, 
    v_reduced, 
    is_linearly_equivalent
