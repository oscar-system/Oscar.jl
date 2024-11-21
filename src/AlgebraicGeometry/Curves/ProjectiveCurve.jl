

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



function edge_matrix(G::Graph)
    M = zero_matrix(ZZ, nv(G), ne(G))
    for (i,e) in enumerate(edges(G))
        M[src(e),i] = 1
        M[dst(e),i] = -1
    end
    return M
end

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
    if check
        @req all(isequal(3),degree(G)) "G is not trivalent"
    end

    R,_ = graded_polynomial_ring(QQ,Int(nv(G)//2+1))
    rowOfVariables = matrix(R,[gens(R)])

    cycleMatrix = kernel(Oscar.edge_matrix(G); side=:right)
    cycleMatrix = matrix(R,cycleMatrix) # converting to matrix over R for vcat below

    edgesG = collect(edges(G)) # indexes all edges of G

    lineIdeals = MPolyIdeal[]
    for v in vertices(G)
        edgesContainingV = findall(edge->(v in edge),edgesG)
        cycleMatrix_v = cycleMatrix[edgesContainingV,:]
        push!(lineIdeals,ideal(minors(vcat(cycleMatrix_v,rowOfVariables),3)))
    end
    graphCurveIdeal = reduce(intersect,lineIdeals)

    return projective_curve(graphCurveIdeal)
end
