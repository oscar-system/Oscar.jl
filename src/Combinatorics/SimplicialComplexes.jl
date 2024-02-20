import Oscar: Polymake, pm_object


################################################################################
##  Constructing
################################################################################

struct SimplicialComplex
    pm_simplicialcomplex::Polymake.BigObject
end

pm_object(K::SimplicialComplex) = K.pm_simplicialcomplex


@doc raw"""
    simplical_complex(generators::Union{Vector{Vector{Int}}, Vector{Set{Int}}})

Construct an abstract simplicial complex from a set of faces.
While arbitrary non-negative integers are allowed as vertices, they will be relabeled to consecutive integers starting at 1.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]])
Abstract simplicial complex of dimension 2 on 4 vertices

julia> G = complete_bipartite_graph(2,3)
Undirected graph with 5 nodes and the following edges:
(3, 1)(3, 2)(4, 1)(4, 2)(5, 1)(5, 2)

julia> K = simplicial_complex(G)
Abstract simplicial complex of dimension 1 on 5 vertices
```

Simplicial complex comprising the empty set only:
```jldoctest
julia> empty = simplicial_complex(Vector{Set{Int}}([]))
Abstract simplicial complex of dimension -1 on 0 vertices
```

The original vertices can be recovered:
```jldoctest
julia> L = simplicial_complex([[0,2,17],[2,17,90]]);

julia> facets(L)
2-element Vector{Set{Int64}}:
 Set([2, 3, 1])
 Set([4, 2, 3])

julia> vertexindices(L)
4-element Vector{Int64}:
  0
  2
 17
 90
```
"""
function simplicial_complex(generators::Union{AbstractVector{<:AbstractVector{<:Base.Integer}}, AbstractVector{<:AbstractSet{<:Base.Integer}}})
  K = Polymake.topaz.SimplicialComplex(INPUT_FACES=generators)
  SimplicialComplex(K)
end

function simplicial_complex(generators::IncidenceMatrix)
  K = Polymake.@convert_to Array{Set} Polymake.common.rows(generators)
  simplicial_complex(K)
end

function simplicial_complex(G::Graph)
  IM = incidence_matrix(G)
  edges_as_rows = IncidenceMatrix(transpose(IM))
  simplicial_complex(edges_as_rows)
end

# more efficient UNEXPORTED+UNDOCUMENTED version, which requires consecutive vertices, and facets as generators;
# will produce errors / segfaults or worse if used improperly
function _SimplicialComplex(generators::Union{Vector{Vector{Int}}, Vector{Set{Int}}})
  K = Polymake.topaz.SimplicialComplex(FACETS=generators)
  SimplicialComplex(K)
end
    
################################################################################
##  Auxiliary
################################################################################

function vertexindices(K::SimplicialComplex) 
    bigobject = pm_object(K)
    if Polymake.exists(bigobject,"VERTEX_INDICES")
        return Vector{Int}(bigobject.VERTEX_INDICES)
    else
        return Vector{Int}(1:bigobject.N_VERTICES)
    end
end

function _convert_finitely_generated_abelian_group(A::Polymake.Array{Polymake.HomologyGroup{Polymake.Integer}},
                                                   h_index::Int)
  # we return non-reduced homology and cohomology
  B = A[h_index + 1] # index shift
  betti_number = is_zero(h_index) ? Polymake.betti_number(B) + 1 : Polymake.betti_number(B)
  vec = zeros(Int, betti_number)
  torsion_i = Polymake.torsion(B)
  for (p,k) in torsion_i
    append!(vec, fill(p,k))
  end
  return abelian_group(vec)
end

################################################################################
##  Properties
################################################################################

@doc raw"""
    n_vertices(K::SimplicialComplex)

Return the number of vertices of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> n_vertices(torus())
7
```
"""
n_vertices(K::SimplicialComplex) = pm_object(K).N_VERTICES::Int

@doc raw"""
    facets(K::SimplicialComplex)

Return the maximal (by inclusion) faces of the abstract simplicial complex `K`.
"""
function facets(K::SimplicialComplex)
    bigobject = pm_object(K)
    the_facets = Polymake.to_one_based_indexing(bigobject.FACETS)
    return Vector{Set{Int}}(the_facets)
end

@doc raw"""
    dim(K::SimplicialComplex)

Return the dimension of the abstract simplicial complex `K`.
"""
dim(K::SimplicialComplex) = pm_object(K).DIM::Int

@doc raw"""
    f_vector(K::SimplicialComplex)

Return the face vector (number of faces per dimension) of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> f_vector(torus())
3-element Vector{Int64}:
  7
 21
 14
```
"""
f_vector(K::SimplicialComplex) = Vector{Int}(pm_object(K).F_VECTOR)

@doc raw"""
    h_vector(K::SimplicialComplex)

Return the h-vector of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> h_vector(torus())
4-element Vector{Int64}:
  1
  4
 10
 -1
```
"""
h_vector(K::SimplicialComplex) = Vector{Int}(pm_object(K).H_VECTOR)

@doc raw"""
    betti_numbers(K::SimplicialComplex)

Return the reduced rational Betti numbers of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> betti_numbers(klein_bottle())
3-element Vector{Int64}:
 0
 1
 0
```
"""
betti_numbers(K::SimplicialComplex) = Vector{Int}(Polymake.topaz.betti_numbers(pm_object(K)))

@doc raw"""
    euler_characteristic(K::SimplicialComplex)

Return the reduced Euler characteristic of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> euler_characteristic(complex_projective_plane())
2
```
"""
euler_characteristic(K::SimplicialComplex) = pm_object(K).EULER_CHARACTERISTIC::Int

@doc raw"""
    homology(K::SimplicialComplex, i::Int)

Return `i`-th integral homology group of `K`.

# Examples
```jldoctest
julia> [ homology(real_projective_plane(), i) for i in [0,1,2] ]
3-element Vector{FinGenAbGroup}:
 Z
 Z/2
 Z/1
```
"""
homology(K::SimplicialComplex, i::Int) = _convert_finitely_generated_abelian_group(pm_object(K).HOMOLOGY, i)

@doc raw"""
    cohomology(K::SimplicialComplex, i::Int)

Return `i`-th integral cohomology group of `K`.

# Examples
```jldoctest
julia> K = simplicial_complex([[0,1],[1,2],[0,2]]);

julia> cohomology(K,1)
Z
```
"""
cohomology(K::SimplicialComplex, i::Int) = _convert_finitely_generated_abelian_group(pm_object(K).COHOMOLOGY, i)

@doc raw"""
    minimal_nonfaces(K::SimplicialComplex)

Return the minimal non-faces of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> minimal_nonfaces(K)
1-element Vector{Set{Int64}}:
 Set([4, 1])
```
"""
minimal_nonfaces(K::SimplicialComplex) = minimal_nonfaces(Vector{Set{Int}}, K)
function minimal_nonfaces(::Type{Vector{Set{Int}}}, K::SimplicialComplex)
    I = minimal_nonfaces(IncidenceMatrix, K)
    return Vector{Set{Int}}([Polymake.row(I,i) for i in 1:Polymake.nrows(I)])
end
function minimal_nonfaces(::Type{IncidenceMatrix}, K::SimplicialComplex)
    # the following line must stay to ensure polymake uses the correct algorithm for the non-faces
    nv = n_vertices(K)
    m = pm_object(K).MINIMAL_NON_FACES
    # fix column number (see #1440) until this is fixed in polymake
    if size(m, 2) < nv
      resize!(m, size(m, 1), nv)
    end
    return m
end

@doc raw"""
    alexander_dual(K::SimplicialComplex)

Return the Alexander dual of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> alexander_dual(K)
Abstract simplicial complex of dimension 1 on 2 vertices
```
"""
alexander_dual(K::SimplicialComplex) = SimplicialComplex(Polymake.topaz.alexander_dual(pm_object(K)))

@doc raw"""
    stanley_reisner_ideal(K::SimplicialComplex)

Return the Stanley-Reisner ideal of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> stanley_reisner_ideal(real_projective_plane())
Ideal generated by
  x1*x2*x3
  x1*x2*x4
  x1*x5*x6
  x2*x5*x6
  x1*x3*x6
  x1*x4*x5
  x3*x4*x5
  x3*x4*x6
  x2*x3*x5
  x2*x4*x6
```
"""
function stanley_reisner_ideal(K::SimplicialComplex)
    n = n_vertices(K)
    R, _ = polynomial_ring(QQ, n, cached=false)
    return stanley_reisner_ideal(R, K)
end

@doc raw"""
    stanley_reisner_ideal(R::MPolyRing, K::SimplicialComplex)

Return the Stanley-Reisner ideal of the abstract simplicial complex `K`, in the given ring `R`.

# Examples
```jldoctest
julia> R, _ = QQ["a","b","c","d","e","f"];

julia> stanley_reisner_ideal(R, real_projective_plane())
Ideal generated by
  a*b*c
  a*b*d
  a*e*f
  b*e*f
  a*c*f
  a*d*e
  c*d*e
  c*d*f
  b*c*e
  b*d*f
```
"""
function stanley_reisner_ideal(R::MPolyRing, K::SimplicialComplex)
    mnf = minimal_nonfaces(IncidenceMatrix, K)
    return ideal([ R([1], [Vector{Int}(mnf[i,:])]) for i in 1:Polymake.nrows(mnf) ])
end

@doc raw"""
    stanley_reisner_ring(K::SimplicialComplex)

Return the Stanley-Reisner ring of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> stanley_reisner_ring(K)
(Quotient of multivariate polynomial ring by ideal (x1*x4), Map: multivariate polynomial ring -> quotient of multivariate polynomial ring)
```
"""
function stanley_reisner_ring(K::SimplicialComplex)
    n = n_vertices(K)
    R, _ = polynomial_ring(QQ, n, cached=false)
    return stanley_reisner_ring(R, K)
end

@doc raw"""
    stanley_reisner_ring(R::MPolyRing, K::SimplicialComplex)

Return the Stanley-Reisner ring of the abstract simplicial complex `K`, as a quotient of a given ring `R`.

# Examples
```jldoctest
julia>  R, _ = ZZ["a","b","c","d","e","f"];

julia> stanley_reisner_ring(R, real_projective_plane())
(Quotient of multivariate polynomial ring by ideal (a*b*c, a*b*d, a*e*f, b*e*f, a*c*f, a*d*e, c*d*e, c*d*f, b*c*e, b*d*f), Map: multivariate polynomial ring -> quotient of multivariate polynomial ring)
```
"""
stanley_reisner_ring(R::MPolyRing, K::SimplicialComplex) = quo(R, stanley_reisner_ideal(R, K))

################################################################################
###  Fundamental group
################################################################################

@doc raw"""
    fundamental_group(K::SimplicialComplex)

Return the fundamental group of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> pi_1 = fundamental_group(torus());

julia> describe(pi_1)
"Z x Z"
```
"""
function fundamental_group(K::SimplicialComplex)
    ngens, relations = pm_object(K).FUNDAMENTAL_GROUP
    F = free_group(ngens)
    nrels = length(relations)
    if nrels==0
        return F
    else
        rvec = Vector{FPGroupElem}(undef, nrels)
        for (i, relation) in enumerate(relations)
            relem = one(F)
            for term in relation
                relem *= F[first(term) + 1]^Int(last(term))
            end
            rvec[i] = relem
        end
        pi_1, _ = quo(F, rvec)
        return pi_1
    end
end

################################################################################
###  Sphere recognition heuristics
################################################################################

"""
    is_sphere(K::SimplicialComplex)

Heuristically check if the abstract simplicial complex `K` is a combinatorial sphere; see [JLLT22](@cite).
Note that this is undecidable in general.
Returns true if recognized as a sphere.
Returns false if not a sphere.
Returns nothing if heuristics unsuccessful.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> is_sphere(K)
false
```
"""
is_sphere(K::SimplicialComplex) = pm_object(K).SPHERE::Union{Bool,Nothing}

"""
    is_ball(K::SimplicialComplex)

Heuristically check if the abstract simplicial complex `K` is a combinatorial ball; see [JLLT22](@cite).
Note that this is undecidable in general.
Returns true if recognized as a ball.
Returns false if not a ball.
Returns nothing if heuristics unsuccessful.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> is_ball(K)
true
```
"""
is_ball(K::SimplicialComplex) = pm_object(K).BALL::Union{Bool,Nothing}

"""
    is_manifold(K::SimplicialComplex)

Check if the abstract simplicial complex `K` is a combinatorial manifold, possibly with boundary.
Note that this is undecidable in general.
Returns true if recognized as a manifold.
Returns false if not a manifold.
Returns nothing if heuristics unsuccessful.

# Examples
```jldoctest
julia> is_manifold(torus())
true

julia> is_manifold(simplicial_complex([[1,2],[2,3]]))
true

julia> is_manifold(simplicial_complex([[1,2],[2,3],[2,4]]))
false
```
"""
is_manifold(K::SimplicialComplex) = pm_object(K).MANIFOLD::Union{Bool,Nothing}

################################################################################
###  Surface examples
################################################################################

"""
    torus()

Construct Möbius' (vertex-minimal) 7-vertex triangulation of the torus (surface).
"""
torus() = SimplicialComplex(Polymake.topaz.torus())

"""
    klein_bottle()

Construct a 9-vertex triangulation of the Klein bottle.
"""
klein_bottle() = SimplicialComplex(Polymake.topaz.klein_bottle())

"""
    real_projective_plane()

Construct the (vertex-minimal) 6-vertex triangulation of the real projective plane.
"""
real_projective_plane() = SimplicialComplex(Polymake.topaz.real_projective_plane())

################################################################################
###  Other examples
################################################################################

"""
    complex_projective_plane()

Construct the (vertex-minimal) 9-vertex triangulation of the complex projective plane.
"""
complex_projective_plane() = SimplicialComplex(Polymake.topaz.complex_projective_plane())

################################################################################
###  Subcomplexes
################################################################################

@doc raw"""
    star_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}})

Return the star of the face `sigma` in the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> star_subcomplex(K,[1])
Abstract simplicial complex of dimension 2 on 3 vertices
```
"""
star_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}}) = SimplicialComplex(Polymake.topaz.star_subcomplex(pm_object(K), Polymake.to_zero_based_indexing(sigma)))

@doc raw"""
    link_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}})

Return the link of the face `sigma` in the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> link_subcomplex(K,[2,3])
Abstract simplicial complex of dimension 0 on 2 vertices
```
"""
link_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}}) = SimplicialComplex(Polymake.topaz.link_subcomplex(pm_object(K), Polymake.to_zero_based_indexing(sigma)))

###############################################################################
### Display
###############################################################################

function Base.show(io::IO, K::SimplicialComplex)
    d = dim(K)
    n = n_vertices(K)
    print(io, "Abstract simplicial complex of dimension $(d) on $(n) vertices")
end

###############################################################################
### Isomorphism
###############################################################################

@doc raw"""
     is_isomorphic(K1::SimplicialComplex, K2::SimplicialComplex)

Checks if the given simplicial complexes are isomorphic.

# Examples
```jldoctest
julia> K1 = simplicial_complex([[1,2,3],[2,3,4]]);

julia> K2 = simplicial_complex([[1,2,3],[2,3,4]]);

julia> is_isomorphic(K1, K2)
true
```
"""
function is_isomorphic(K1::SimplicialComplex, K2::SimplicialComplex)
  return Polymake.topaz.isomorphic(pm_object(K1), pm_object(K2))::Bool
end

