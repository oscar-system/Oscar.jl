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
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
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

function simplicial_complex(generators::Union{Set{<:AbstractVector{<:Base.Integer}}, Set{<:AbstractSet{<:Base.Integer}}})
  return simplicial_complex(collect(generators))
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
    n_facets(K::SimplicialComplex)

Return the number of facets of the abstract simplicial complex `K`.
"""
n_facets(K::SimplicialComplex) = pm_object(K).N_FACETS::Int

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
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
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
julia> R, _ = QQ[:a, :b, :c, :d, :e, :f];

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
  gens = [ R([1], [Vector{Int}(mnf[i,:])]) for i in 1:Polymake.nrows(mnf) ]
  # currently no way to set as a universal groebner basis
  I = IdealGens(R, gens, default_ordering(R); isGB=true, isReduced=true)
  return ideal(I)
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
julia>  R, _ = ZZ[:a, :b, :c, :d, :e, :f];

julia> stanley_reisner_ring(R, real_projective_plane())
(Quotient of multivariate polynomial ring by ideal (a*b*c, a*b*d, a*e*f, b*e*f, a*c*f, a*d*e, c*d*e, c*d*f, b*c*e, b*d*f), Map: R -> quotient of multivariate polynomial ring)
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
Return `true` if recognized as a sphere.
Return `false` if not a sphere.
Return `nothing` if heuristics unsuccessful.

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
Return `true` if recognized as a ball.
Return `false` if not a ball.
Return `nothing` if heuristics unsuccessful.

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
Return `true` if recognized as a manifold.
Return `false` if not a manifold.
Return `nothing` if heuristics unsuccessful.

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

Check if the given simplicial complexes are isomorphic.

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

###############################################################################
### helpful functions
###############################################################################

@doc raw"""
    connected_sum(K1::SimplicialComplex, K2::SimplicialComplex, f1::Int=0, f2::Int=0)

Compute the connected sum of two abstract simplicial complexes. Parameters `f1` and `f2` specify which facet
 of the first and second complex correspondingly are glued together. Default is the
0-th facet of both. The vertices in the selected facets are identified with each
other according to their order in the facet (that is, in increasing index order).

# Examples
```jldoctest
julia> K = torus();

julia> surface_genus_2 = connected_sum(K, K)
Abstract simplicial complex of dimension 2 on 11 vertices

julia> homology(surface_genus_2, 1)
Z^4

julia> is_manifold(surface_genus_2)
true
```
"""
function connected_sum(K1::SimplicialComplex, K2::SimplicialComplex, f1::Int=0, f2::Int=0)
  return SimplicialComplex(Polymake.topaz.connected_sum(pm_object(K1), pm_object(K2), f1, f2))
end

@doc raw"""
    deletion(K::SimplicialComplex, face::Union{<:AbstractSet{Int},<:AbstractVector{Int}})

Remove the given face and all the faces containing it from an abstract simplicial complex `K`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> K = simplicial_complex([[1, 2, 3], [2, 3, 4]]);

julia> K_with_deletion = deletion(K, Set([1, 2]));

julia> facets(K_with_deletion)
2-element Vector{Set{Int64}}:
 Set([3, 1])
 Set([4, 2, 3])
```
"""
function deletion(K::SimplicialComplex, face::Union{<:AbstractSet{Int},<:AbstractVector{Int}})
  zero_based = Polymake.to_zero_based_indexing(face)
  return SimplicialComplex(Polymake.topaz.deletion(pm_object(K), zero_based))
end

@doc raw"""
    automorphism_group(K::SimplicialComplex; action=:on_vertices)

Given a simplicial complex `K` return its automorphism group as a `PermGroup`.
The group can be returned as a subgroup of the permutation group of the vertices
by passing `:on_vertices` to the `action` keyword argument or on the facets
 by passing `:on_facets`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1, 2, 3], [2, 3, 4]])
Abstract simplicial complex of dimension 2 on 4 vertices

julia> automorphism_group(K)
Permutation group of degree 4 with 2 generators
  (2,3)
  (1,4)
```
"""
function automorphism_group(K::SimplicialComplex; action=:on_vertices)
  pm_K = Oscar.pm_object(K)
  Polymake.topaz.combinatorial_symmetries(pm_K)
  if action == :on_vertices
    gens_G = Polymake.to_one_based_indexing(pm_K.GROUP.RAYS_ACTION.GENERATORS)
    n = n_vertices(K)
    return permutation_group(n, perm.(gens_G))
  elseif action == :on_facets
    gens_G = Polymake.to_one_based_indexing(pm_K.GROUP.FACETS_ACTION.GENERATORS)
    n = n_facets(K)
    return permutation_group(n, perm.(gens_G))
  else
    error("unsupported keyword passed to action")
  end
end

@doc raw"""
    on_simplicial_complex(K::SimplicialComplex, g::PermGroupElem)

Given a simplicial complex `K` return the simplicial complex corresponding
to a permutation on it's vertices given by `g`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> K = simplicial_complex([[1, 2, 3], [2, 3, 4]])
Abstract simplicial complex of dimension 2 on 4 vertices

julia> G = automorphism_group(K)
Permutation group of degree 4 with 2 generators
  (2,3)
  (1,4)

julia> g = collect(G)[2]
(1,4)

julia> facets(on_simplicial_complex(K, g))
2-element Vector{Set{Int64}}:
 Set([2, 3, 1])
 Set([4, 2, 3])
```
"""
function on_simplicial_complex(K::SimplicialComplex, g::PermGroupElem)
  @req degree(parent(g)) == n_vertices(K) "g needs to be an element of the permutation group on the vertices"
  new_facets = on_sets_sets(Set(facets(K)), g)
  simplicial_complex(collect(new_facets))
end

@doc raw"""
     simplicial_product(K1::SimplicialComplex, K2::SimplicialComplex)

Given simplicial complexes `K1` and `K2` return their simplicial product.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> K1 = simplicial_complex([[1, 2], [2, 3]])
Abstract simplicial complex of dimension 1 on 3 vertices

julia> K2 = simplicial_complex([[1, 2], [1, 3]])
Abstract simplicial complex of dimension 1 on 3 vertices

julia> facets(simplicial_product(K1, K2))
8-element Vector{Set{Int64}}:
 Set([5, 2, 1])
 Set([5, 4, 1])
 Set([2, 8, 1])
 Set([7, 8, 1])
 Set([6, 2, 3])
 Set([5, 6, 2])
 Set([2, 9, 3])
 Set([2, 9, 8])
```
"""
function simplicial_product(K1::SimplicialComplex, K2::SimplicialComplex)
  return SimplicialComplex(Polymake.topaz.simplicial_product(pm_object(K1), pm_object(K2)))
end

@doc raw"""
     link_subcomplex(K::SimplicialComplex, f::Union{<:AbstractSet{Int},<:AbstractVector{Int}}))

Given simplicial complex `K` and a face `f` of `K` return the link of `f`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1, 2], [2, 3], [3, 4]])
Abstract simplicial complex of dimension 1 on 4 vertices

julia> facets(link_subcomplex(K, [2]))
2-element Vector{Set{Int64}}:
 Set([1])
 Set([2])
```
"""
function link_subcomplex(K::SimplicialComplex, face::Union{<:AbstractSet{Int},<:AbstractVector{Int}})
  zero_based_face = Polymake.to_zero_based_indexing(face)
  return SimplicialComplex(Polymake.link_subcomplex(pm_object(K), zero_based_face))
end

@doc raw"""
   barycentric_subdivision(K::SimplicialComplex)

Given simplicial complex `K` returns its barycentric subdivision.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> K = simplicial_complex([[1, 2, 3]])
Abstract simplicial complex of dimension 2 on 3 vertices

julia> facets(barycentric_subdivision(K))
6-element Vector{Set{Int64}}:
 Set([5, 2, 1])
 Set([5, 3, 1])
 Set([6, 2, 1])
 Set([4, 6, 1])
 Set([7, 3, 1])
 Set([4, 7, 1])
```
"""
function barycentric_subdivision(K::SimplicialComplex)
  return SimplicialComplex(Polymake.topaz.barycentric_subdivision(pm_object(K))) 
end

function is_shifted(K::SimplicialComplex)
  return pm_object(K).SHIFTED
end
