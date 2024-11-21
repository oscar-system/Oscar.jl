

@doc raw"""
    ProjectiveCurve

A reduced projective curve, defined as the vanishing locus of
a homogeneous (but not necessarily radical) ideal.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);

julia> M = matrix(R, 2, 3, [w x y; x y z])
[w   x   y]
[x   y   z]

julia> V = minors(M, 2)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 w*y - x^2
 w*z - x*y
 x*z - y^2

julia> I = ideal(R, V);

julia> TC = projective_curve(I)
Projective curve
  in projective 3-space over QQ with coordinates [w, x, y, z]
defined by ideal (w*y - x^2, w*z - x*y, x*z - y^2)

```
"""
@attributes mutable struct ProjectiveCurve{BaseRingType<:Field, RingType<:Ring} <: AbsProjectiveCurve{BaseRingType, RingType}
    X::ProjectiveAlgebraicSet{BaseRingType, RingType}

    function ProjectiveCurve(X::ProjectiveAlgebraicSet{S,T}; check::Bool=true) where {S,T}
        @check dim(X) == 1 || error("not of dimension one")
        new{S,T}(X)
    end
end

ProjectiveCurve(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true, is_radical::Bool=true) = ProjectiveCurve(algebraic_set(I;check=check,is_radical=is_radical);check=check)

projective_curve(I;kwargs...) = ProjectiveCurve(I;kwargs...)

underlying_scheme(X::ProjectiveCurve) = X.X
fat_scheme(X::ProjectiveCurve) = fat_scheme(underlying_scheme(X))


function Base.show(io::IO, ::MIME"text/plain", X::ProjectiveCurve)
  io = pretty(io)
  println(io, "Projective curve")
  println(io, Indent(), "in ", Lowercase(), ambient_space(X))
  if isdefined(X, :Xred)
    I = vanishing_ideal(X)
  else
    I = fat_ideal(X)
  end
  print(io, Dedent(), "defined by ", Lowercase(), I)
end

################################################################################

@doc raw"""
    invert_birational_map(phi::Vector{T}, C::ProjectiveCurve) where {T <: MPolyRingElem}

Return a dictionary where `image` represents the image of the birational map
given by `phi`, and `inverse` represents its inverse, where `phi` is a
birational map of the projective curve `C` to its image in the projective
space of dimension `size(phi) - 1`.
Note that the entries of `inverse` should be considered as
representatives of elements in `R/image`, where `R` is the basering.
"""
function invert_birational_map(phi::Vector{T}, C::ProjectiveCurve) where {T <: MPolyRingElem}
    s = parent(phi[1])
    I = ideal(s, phi)
    IC = defining_ideal(C)
    L = Singular.LibParaplanecurves.invertBirMap(singular_generators(I), singular_generators(IC))
    R = _fromsingular_ring(L[1])
    J = L[2][:J]
    psi = L[2][:psi]
    return Dict(:image => gens(ideal(R, J)), :inverse => gens(ideal(R, psi)))
end


################################################################################
# The union of two curves is a curve
Base.union(C::T, D::T) where T <: AbsProjectiveCurve = ProjectiveCurve(union(underlying_scheme(C),underlying_scheme(D)))



@doc raw"""
    graph_curve(G::Graph)

Return the graph curve of `G`, i.e., a union of lines whose dual graph is `G`, see [BE91](@cite).

Assumes that `G` is trivalent and 3-connected.

# Examples
```jldoctest
julia> G1 = edgegraph(simplex(3))
Undirected graph with 4 nodes and the following edges:
(2, 1)(3, 1)(3, 2)(4, 1)(4, 2)(4, 3)

julia> C1 = graph_curve(G1)
Projective curve
  in projective 2-space over QQ with coordinates [x1, x2, x3]
defined by ideal (x1^2*x2*x3 - x1*x2^2*x3 + x1*x2*x3^2)

julia> G2 = edgegraph(cube(3))
Undirected graph with 8 nodes and the following edges:
(2, 1)(3, 1)(4, 2)(4, 3)(5, 1)(6, 2)(6, 5)(7, 3)(7, 5)(8, 4)(8, 6)(8, 7)

julia> C2 = graph_curve(G2)
Projective curve
  in projective 4-space over QQ with coordinates [x1, x2, x3, x4, x5]
defined by ideal with 5 generators

```
"""
function graph_curve(G::Graph; check::Bool=true)
    C,_ = graph_curve_with_vertex_ideals(G; check=check)
    return C
end

function graph_curve_with_vertex_ideals(G::Graph; check::Bool=true)
    if check
        @req all(isequal(3),degree(G)) "G is not trivalent"
    end

    R,_ = graded_polynomial_ring(QQ,Int(nv(G)//2+1))
    rowOfVariables = matrix(R,[gens(R)])

    cycleMatrix = kernel(matrix(QQ,signed_incidence_matrix(G)); side=:right)
    cycleMatrix = matrix(R,cycleMatrix) # converting to matrix over R for vcat below

    edgesG = collect(edges(G)) # indexes all edges of G

    vertexIdeals = MPolyIdeal[]
    for v in vertices(G)
        edgesContainingV = findall(edge->(v in edge),edgesG)
        cycleMatrix_v = cycleMatrix[edgesContainingV,:]
        push!(vertexIdeals,ideal(minors(vcat(cycleMatrix_v,rowOfVariables),3)))
    end
    graphCurveIdeal = reduce(intersect,vertexIdeals)
    graphCurve = projective_curve(graphCurveIdeal)

    return graphCurve, vertexIdeals
end

function graph_curve_with_vertex_and_edge_ideals(G::Graph)
    C, vertexIdeals = graph_curve_with_vertex_ideals(G)
    edgeIdeals = [ vertexIdeals[src(e)] + vertexIdeals[dst(e)] for e in edges(G)]
    return C, vertexIdeals, radical.(edgeIdeals)  # is radical really necessary?
end


function edge_pairing_candidates(G::Graph, e::Edge)
    v1, v2 = src(e), dst(e)
    adj_v1 = filter(v -> v != v2, neighbors(G,v1))
    adj_v2 = filter(v -> v != v1, neighbors(G,v2))
    pairing1 = ( (Edge(v1,adj_v1[1]),Edge(v2,adj_v2[1])), (Edge(v1,adj_v1[2]),Edge(v2,adj_v2[2])) )
    pairing2 = ( (Edge(v1,adj_v1[1]),Edge(v2,adj_v2[2])), (Edge(v1,adj_v1[1]),Edge(v2,adj_v2[2])) )
    return [pairing1,pairing2]
end

function maximal_edge_pairing(G::Graph)
    pGfaces = faces(IncidenceMatrix, tutte_lifting(G), 2)
    EP = Dict{Edge, Tuple{Tuple{Edge, Edge},Tuple{Edge, Edge}}}()
    for edge in edges(G)
        pairing1, pairing2 =  edge_pairing_candidates(G,edge)
        # take the first pair of edges of each pairing candidate
        # the second pair of edges is not required
        # as the first pair is correct if and only if the second pair is
        edge11 = pairing1[1][1]
        edge12 = pairing1[1][2]
        edge21 = pairing2[1][1]
        edge22 = pairing2[1][2]
        # for each pair of edges, construct the set of four vertices
        pairingVertices1 = [src(edge11),dst(edge11),src(edge12),dst(edge12)]
        pairingVertices2 = [src(edge21),dst(edge21),src(edge22),dst(edge22)]
        # and check whether there is a tutte lifting facet containing all of them
        for i in 1:nrows(pGfaces)
            if all(pGfaces[i,pairingVertices1])
                EP[edge] = pairing1
                break
            elseif all(pGfaces[i,pairingVertices2])
                EP[edge] = pairing2
                break
            end
        end
    end
    return EP
end
