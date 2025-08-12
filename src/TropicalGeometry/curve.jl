###############################################################################
#
#  Tropical curves
#  ===============
#  concrete subtype of TropicalVarietySupertype in variety_supertype.jl
#
###############################################################################

@attributes mutable struct TropicalCurve{minOrMax,isEmbedded} <: TropicalVarietySupertype{minOrMax,isEmbedded}
    polyhedralComplex::Union{PolyhedralComplex, Graph}
    multiplicities::Vector{ZZRingElem}

    # embedded tropical curves contain a PolyhedralComplex
    function TropicalCurve{minOrMax,true}(Sigma::PolyhedralComplex, multiplicities::Vector{ZZRingElem}) where {minOrMax<:Union{typeof(min),typeof(max)}}
        @req dim(Sigma)==1 "input not one-dimensional"
        return new{minOrMax,true}(Sigma,multiplicities)
    end

    # abstract tropical curves contain a Graph
    function TropicalCurve{minOrMax,false}(Sigma::Graph, multiplicities::Vector{ZZRingElem}) where {minOrMax<:Union{typeof(min),typeof(max)}}
        return new{minOrMax,false}(Sigma,multiplicities)
    end
end



###############################################################################
#
#  Printing
#
###############################################################################

function Base.show(io::IO, th::TropicalCurve{typeof(min),true})
    print(io, "Embedded min tropical curve")
end
function Base.show(io::IO, th::TropicalCurve{typeof(max),true})
    print(io, "Embedded max tropical curve")
end
function Base.show(io::IO, th::TropicalCurve{typeof(min),false})
    print(io, "Abstract min tropical curve")
end
function Base.show(io::IO, th::TropicalCurve{typeof(max),false})
    print(io, "Abstract max tropical curve")
end


###############################################################################
#
#  Constructors
#
###############################################################################

@doc raw"""
    tropical_curve(Sigma::PolyhedralComplex, multiplicities::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)

Return the embedded tropical curve consisting of the polyhedral complex `Sigma` and multiplicities `multiplicities`.

# Examples
```jldoctest
julia> verticesAndRays = [0 0; 1 0; 0 1; -1 -1];

julia> incidenceMatrix = incidence_matrix([[1,2],[1,3],[1,4]]);

julia> rayIndices = [2,3,4];

julia> Sigma = polyhedral_complex(incidenceMatrix, verticesAndRays, rayIndices)
Polyhedral complex in ambient dimension 2

julia> multiplicities = ones(ZZRingElem, n_maximal_polyhedra(Sigma))
3-element Vector{ZZRingElem}:
 1
 1
 1

julia> tropical_curve(Sigma,multiplicities)
Embedded min tropical curve

```
"""
function tropical_curve(Sigma::PolyhedralComplex, multiplicities::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)
    return TropicalCurve{typeof(minOrMax),true}(Sigma,multiplicities)
end
function tropical_curve(Sigma::PolyhedralComplex, minOrMax::Union{typeof(min),typeof(max)}=min)
    multiplicities = ones(ZZRingElem, n_maximal_polyhedra(Sigma))
    return tropical_curve(Sigma,multiplicities,minOrMax)
end


@doc raw"""
    tropical_curve(Sigma::Graph, multiplicities::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)

Return the abstract tropical curve consisting of the graph `Sigma` and multiplicities `multiplicities`.

# Examples
```jldoctest; filter = r"Edge\(.*\)"
julia> Sigma = graph_from_adjacency_matrix(Undirected,[0 1 1; 1 0 1; 1 1 0]);

julia> multiplicities = ones(ZZRingElem, n_edges(Sigma))
3-element Vector{ZZRingElem}:
 1
 1
 1

julia> tropical_curve(Sigma,multiplicities)
Abstract min tropical curve

```
"""
function tropical_curve(Sigma::Graph, multiplicities::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)
   return TropicalCurve{typeof(minOrMax), false}(Sigma,multiplicities)
end
function tropical_curve(Sigma::Graph,minOrMax::Union{typeof(min),typeof(max)}=min)
    multiplicities = ones(ZZRingElem, n_edges(Sigma))
    return tropical_curve(Sigma,multiplicities,minOrMax)
end


function tropical_curve(TropV::TropicalVarietySupertype)
    @req dim(TropV)<=1 "tropical variety dimension too high"
    return tropical_curve(polyhedral_complex(TropV), multiplicities(TropV), convention(TropV))
end



################################################################################
#
#  Properties
#
################################################################################

@doc raw"""
    graph(TropC::TropicalCurve{minOrMax,false})

Return the graph of an abstract tropical curve `TropC`.  Same as `polyhedral_complex(tc)`.
"""
function graph(TropC::TropicalCurve{minOrMax,false}) where minOrMax
    return TropC.polyhedralComplex
end



################################################################################
#
#  Outdated code (to be updated)
#
################################################################################

# @doc raw"""
#     n_vertices(tc::TropicalCurve)

# Return the number of vertices of an abstract tropical curve `tc`.
# """
# function n_vertices(tc::TropicalCurve)
#     G = graph(tc)
#     return n_vertices(G)
# end


# struct DivisorOnTropicalCurve{M, EMB}
#     base_curve::TropicalCurve{M, EMB}
#     coefficients::Vector{Int}

#     DivisorOnTropicalCurve{M, EMB}(tc::TropicalCurve{M, EMB}, coeffs::Vector{Int}) where {M,EMB} = new{M, EMB}(tc, coeffs)
# end


# @doc raw"""
#     divisor_on_tropical_curve(tc::TropicalCurve, coeffs::Vector{Int})

# Construct a divisor with coefficients `coeffs` on an abstract tropical curve `tc`.

# # Examples
# ```jldoctest
# julia> Sigma = graph_from_adjacency_matrix(Undirected,[0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0]);

# julia> tc = tropical_curve(Sigma)
# Abstract min tropical curve

# julia> dtc = divisor_on_tropical_curve(tc,[0, 1, 1, 1])
# DivisorOnTropicalCurve{min, false}(Abstract min tropical curve, [0, 1, 1, 1])
# ```
# """
# function divisor_on_tropical_curve(tc::TropicalCurve{minOrMax, false}, coeffs::Vector{Int}) where minOrMax
#     @req n_vertices(tc)==length(coeffs) "Wrong number coefficients"
#     return DivisorOnTropicalCurve{minOrMax, false}(tc, coeffs)
# end

# base_curve(dtc::DivisorOnTropicalCurve) = dtc.base_curve


# ###
# # 3.Basic properties
# ###
# @doc raw"""
#     coefficients(dtc::DivisorOnTropicalCurve)

# Return the coefficients of `dtc`.
# """
# coefficients(dtc::DivisorOnTropicalCurve) = dtc.coefficients


# @doc raw"""
#    degree(dtc::DivisorOnTropicalCurve)

# Return the degree of `dtc`.
# """
# degree(dtc::DivisorOnTropicalCurve) = sum(coefficients(dtc))


# @doc raw"""
#     is_effective(dtc::DivisorOnTropicalCurve)

# Return `true` if all coefficients of `dtc` are positive, `false` otherwise.
# """
# is_effective(dtc::DivisorOnTropicalCurve) = all(>=(0), coefficients(dtc))


# @doc raw"""
#    chip_firing_move(dtc::DivisorOnTropicalCurve, position::Int)

# Given a divisor `dtc` and vertex labeled `position`, compute the linearly equivalent divisor obtained by a chip firing move from the given vertex `position`.

# # Examples
# ```jldoctest
# julia> IM = incidence_matrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

# julia> tc = TropicalCurve(IM)
# ERROR: MethodError: no method matching TropicalCurve(::Polymake.LibPolymake.IncidenceMatrixAllocated{Polymake.LibPolymake.NonSymmetric})
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> coeffs = [0, 1, 1, 1];

# julia> dtc = divisor_on_tropical_curve(tc,coeffs)
# ERROR: UndefVarError: `tc` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> chip_firing_move(dtc,1)
# ERROR: UndefVarError: `dtc` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1
# ```
# """
# function chip_firing_move(dtc::DivisorOnTropicalCurve, position::Int)
#     G = graph(base_curve(dtc))
#     newcoeffs = Vector{Int}(coefficients(dtc))
#     for i in 1:n_edges(G)
#         row = Polymake.row(incidence_matrix(G), i-1)
#         if position in row
#             newcoeffs[position] -= 1
#             for i in row
#                 if i != position
#                     newcoeffs[i] += 1
#                 end
#             end
#         end
#     end
#     return divisor_on_tropical_curve(base_curve(dtc), newcoeffs)
# end

# ### The function computes the outdegree of a vertex v  with respect to a given subset W  of vertices.
# ### This is the number of vertices not in W adjacent to v. 1,
# function outdegree(tc::TropicalCurve, W::Set{Int}, v::Int)
#     G = graph(tc)
#     m = n_edges(G) #number of edges of tc
#     @req v in W "Vertex number $v not in $W"
#     deg = 0 #outdeg
#     for i in 1:m
#         row = Polymake.row(incidence_matrix(G),i)
#         if v in row
#             for j in row
#                 if !(j in W)
#                     deg = deg +1
#                 end
#             end
#         end
#     end
#     return deg
# end

# @doc raw"""
#    v_reduced(dtc::DivisorOnTropicalCurve, vertex::Int)

# Given a divisor `dtc` and vertex labeled `vertex`, compute the unique divisor reduced with repspect to `vertex`
# as defined in [BN07](@cite).
# The divisor `dtc` must have positive coefficients apart from `vertex`.

# # Examples
# ```jldoctest
# julia> IM = incidence_matrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

# julia> tc = TropicalCurve(IM)
# ERROR: MethodError: no method matching TropicalCurve(::Polymake.LibPolymake.IncidenceMatrixAllocated{Polymake.LibPolymake.NonSymmetric})
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> coeffs = [0, 1, 1, 1];

# julia> dtc = divisor_on_tropical_curve(tc,coeffs)
# ERROR: UndefVarError: `tc` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> v_reduced(dtc,1)
# ERROR: UndefVarError: `dtc` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1
# ```
# """
# function v_reduced(dtc::DivisorOnTropicalCurve, vertex::Int)
#     tc = base_curve(dtc)
#     G = graph(base_curve(dtc))
#     n = n_vertices(G)
#     newcoeff = Vector{Int}(coefficients(dtc))
#     S = Set{Int}(1:n)
#     W = setdiff(S,vertex)
#     @req all(j -> newcoeff[j]>=0, W) "Divisor not effective outside the vertex number $vertex"
#     while !(isempty(W))
# 	      w0 = vertex
# 	      if all(w -> newcoeff[w] >= outdegree(tc,W,w),W)
# 	          coeffdelta = zeros(Int, n)
# 	          delta = divisor_on_tropical_curve(tc,coeffdelta)
# 	          for j in W
# 	              delta = chip_firing_move(delta,j)
# 	          end
# 	          newcoeff = newcoeff + coefficients(delta)
# 	      else
# 	          for w in W
# 	              if newcoeff[w] < outdegree(tc,W,w)
# 		                w0 = w
# 		                break
# 	              end
# 	          end
# 	          W = setdiff(W,w0)
#         end
#     end
#     reduced = divisor_on_tropical_curve(tc, newcoeff)
#     return reduced
# end

# @doc raw"""
#    is_linearly_equivalent(dtc1::DivisorOnTropicalCurve, dtc2::DivisorOnTropicalCurve)

# Given two effective divisors `dtc1` and `dtc2` on the same tropical curve, check whether they are linearly equivalent.

# # Examples
# ```jldoctest
# julia> IM = incidence_matrix([[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]);

# julia> tc = TropicalCurve(IM)
# ERROR: MethodError: no method matching TropicalCurve(::Polymake.LibPolymake.IncidenceMatrixAllocated{Polymake.LibPolymake.NonSymmetric})
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> coeffs1 = [0, 1, 1, 1];

# julia> dtc1 = divisor_on_tropical_curve(tc,coeffs1)
# ERROR: UndefVarError: `tc` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> coeffs2 = [3,0,0,0];

# julia> dtc2 = divisor_on_tropical_curve(tc,coeffs2)
# ERROR: UndefVarError: `tc` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> is_linearly_equivalent(dtc1, dtc2)
# ERROR: UndefVarError: `dtc1` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1
# ```
# """
# function is_linearly_equivalent(dtc1::DivisorOnTropicalCurve, dtc2::DivisorOnTropicalCurve)
#     @req is_effective(dtc1) "The divisor $dtc1 is not effective"
#     @req is_effective(dtc2) "The divisor $dtc2 is not effective"
#     @req base_curve(dtc1) === base_curve(dtc2) "The input curve needs to be the same for both divisors"
#     v = 1
#     reduced1 = v_reduced(dtc1,v)
#     reduced2 = v_reduced(dtc2,v)
#     coeff1 = coefficients(reduced1)
#     coeff2 = coefficients(reduced2)
#     return coeff1 == coeff2
# end


# ###
# # 4. More properties
# # -------------------
# ###

# @doc raw"""
#     structure_tropical_jacobian(TropC::TropicalCurve)

# Compute the elementary divisors $n_i$ of the Laplacian matrix of the tropical curve `TropC`.
# The tropical Jacobian is then isomorphic to $\prod (Z/(n_i)Z)$.

# # Examples
# ```jldoctest
# julia> cg = complete_graph(5);

# julia> IM1 = incidence_matrix([[src(e), dst(e)] for e in edges(cg)])
# 10×5 IncidenceMatrix
# [1, 2]
# [1, 3]
# [2, 3]
# [1, 4]
# [2, 4]
# [3, 4]
# [1, 5]
# [2, 5]
# [3, 5]
# [4, 5]

# julia> TropC1 = TropicalCurve(IM1)
# ERROR: MethodError: no method matching TropicalCurve(::Polymake.LibPolymake.IncidenceMatrixAllocated{Polymake.LibPolymake.NonSymmetric})
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> structure_tropical_jacobian(TropC1)
# ERROR: UndefVarError: `TropC1` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> cg2 = complete_graph(3);

# julia> IM2 = incidence_matrix([[src(e), dst(e)] for e in edges(cg2)])
# 3×3 IncidenceMatrix
# [1, 2]
# [1, 3]
# [2, 3]

# julia> TropC2 = TropicalCurve(IM2)
# ERROR: MethodError: no method matching TropicalCurve(::Polymake.LibPolymake.IncidenceMatrixAllocated{Polymake.LibPolymake.NonSymmetric})
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> structure_tropical_jacobian(TropC2)
# ERROR: UndefVarError: `TropC2` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> IM3 = incidence_matrix([[1,2],[2,3],[3,4],[4,5],[1,5]])
# 5×5 IncidenceMatrix
# [1, 2]
# [2, 3]
# [3, 4]
# [4, 5]
# [1, 5]

# julia> TropC3=TropicalCurve(IM3)
# ERROR: MethodError: no method matching TropicalCurve(::Polymake.LibPolymake.IncidenceMatrixAllocated{Polymake.LibPolymake.NonSymmetric})
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> G = structure_tropical_jacobian(TropC3)
# ERROR: UndefVarError: `TropC3` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1
# ```
# """
# function structure_tropical_jacobian(TropC::TropicalCurve)
#     gg=Graph{Undirected}(n_vertices(TropC))
#     IM = graph(TropC)
#     for i in 1:n_edges(IM)
#         row = Vector{Int}(Polymake.row(incidence_matrix(IM),i))
#         add_edge!(gg, row[1],row[2])
#     end
#     lap = Polymake.graph.laplacian(Oscar.pm_object(gg))
#     L = Polymake.@convert_to Matrix{Int} lap
#     LL = matrix(ZZ, L)
#     ED = elementary_divisors(LL)[1:nrows(LL)-1]
#     G = abelian_group(ED)
#     return G
# end


# # """
# # This recipe allows to use `Plots` without having `Plots` as a dependency.

# # Usage example:
# # ```julia
# # julia> using Revise, Plots, Oscar, Test;

# # julia> IM = incidence_matrix([[1,2],[1,3],[1,4]]);

# # julia> VR = [0 0; 1 0; -1 0; 0 1];

# # julia> PC = PolyhedralComplex{QQFieldElem}(IM, VR);

# # julia> TropC = TropicalCurve(PC);

# # julia> plot(TropC)
# # ```
# # """
# # @recipe function visualize(tc::TropicalCurve{M,EMB}) where {M,EMB}
# #     @req EMB "Tropical curve is abstract."
# #     PC = tc.polyhedralComplex
# #     MaxPoly= maximal_polyhedra(PC)
# #     list_vertices = Vector{Complex{Float64}}()
# #     for P in MaxPoly
# #         V = vertices(P)
# #         R = rays(P)
# #         #V = [Vector{Rational{Int}}(x) for x in V]
# #         V = [Vector{Float64}(Vector{Rational}(Vector{Polymake.Rational}(x))) for x in V]
# #         #R = [Vector{Rational{Int}}(x) for x in R]
# #         R = [Vector{Float64}(Vector{Rational}(Vector{Polymake.Rational}(x))) for x in R]
# #         if length(V)==2
# #             append!(list_vertices,V[1][1]+V[1][2]*im,V[2][1]+V[2][2]*im,Inf+0*im)
# #         else
# #             B = V[1]+R[1]
# #             append!(list_vertices,V[1][1]+V[1][2]*im,B[1]+B[2]*im,Inf+0*im)
# #         end
# #     end
# #     # Set default options for plot
# #     legend := false
# #     axis := false
# #     xlabel := ""
# #     ylabel := ""
# #     return list_vertices
# # end
