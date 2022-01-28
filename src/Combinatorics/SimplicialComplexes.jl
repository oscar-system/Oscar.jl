using Markdown
import Oscar: Polymake, pm_object

export
    SimplicialComplex,
    facets, nvertices,
    vertexindices,
    f_vector, h_vector,
    dim,
    betti_numbers, euler_characteristic, homology, cohomology,
    fundamental_group,
    minimal_nonfaces, alexander_dual, stanley_reisner_ideal, stanley_reisner_ring,
    real_projective_plane, klein_bottle, torus, # requires a distinction from, e.g., an algebraic group
    complex_projective_plane,
    star_subcomplex, link_subcomplex,
    load_simplicialcomplex, save_simplicialcomplex

################################################################################
##  Constructing
################################################################################

struct SimplicialComplex
    pm_simplicialcomplex::Polymake.BigObject
end

pm_object(K::SimplicialComplex) = K.pm_simplicialcomplex


@doc Markdown.doc"""
    SimplicialComplex(generators::Union{Vector{Vector{Int}}, Vector{Set{Int}}})

Construct an abstract simplicial complex from a set of faces.
While arbitrary nonnegative integers are allowed as vertices, they will be relabeled to consecutive integers starting at 1.

# Example
```jldoctest
julia> K = SimplicialComplex([[1,2,3],[2,3,4]])
Abstract simplicial complex of dimension 2 on 4 vertices
```
# Simplicial complex comprising the empty set only
```jldoctest
julia> empty = SimplicialComplex(Vector{Set{Int}}([]))
Abstract simplicial complex of dimension -1 on 0 vertices
```
# Example with relabeling
The original vertices can be recovered.
```jldoctest
julia> L = SimplicialComplex([[0,2,17],[2,17,90]]);

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
function SimplicialComplex(generators::Union{AbstractVector{<:AbstractVector{<:Base.Integer}}, AbstractVector{<:AbstractSet{<:Base.Integer}}})
    K = Polymake.topaz.SimplicialComplex(INPUT_FACES=generators)
    SimplicialComplex(K)
end

function SimplicialComplex(generators::IncidenceMatrix)
    K = Polymake.@convert_to Array{Set} Polymake.common.rows(generators)
    SimplicialComplex(K)
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

function _vertexindices(K::Polymake.BigObject)
    if Polymake.exists(K,"VERTEX_INDICES")
        return Vector{Int}(K.VERTEX_INDICES)
    else
        return Vector{Int}(1:K.N_VERTICES)
    end
end

vertexindices(K::SimplicialComplex) = _vertexindices(pm_object(K))

# currently unused
_reindexset(M::Set{Int}, ind::Vector{Int}) = [ ind[x] for x in M ]

function _characteristic_vector(M::Set{Int}, n::Int)
    chi = zeros(Int, n)
    for x in M
        chi[x] = 1
    end
    return chi
end

function _convert_finitely_generated_abelian_group(A::Polymake.HomologyGroupAllocated{Polymake.Integer})
    vec = ones(Int, Polymake.betti_number(A))
    torsion_i = Polymake.torsion(A)
    for (p,k) in torsion_i
        append!(vec, fill(p,k))
    end
    return vec
end

################################################################################
##  Properties
################################################################################

@doc Markdown.doc"""
    nvertices(K::SimplicialComplex)

Return the number of vertices of the abstract simplicial complex `K`.

# Example
```jldoctest
julia> nvertices(torus())
7
```
"""
nvertices(K::SimplicialComplex) = pm_object(K).N_VERTICES::Int

@doc Markdown.doc"""
    facets(K::SimplicialComplex)

Return the maximal (by inclusion) faces of the abstract simplicial complex `K`.
"""
function facets(K::SimplicialComplex)
    bigobject = pm_object(K)
    the_facets = Polymake.to_one_based_indexing(bigobject.FACETS)
    return Vector{Set{Int}}(the_facets)
end

@doc Markdown.doc"""
    dim(K::SimplicialComplex)

Return the dimension of the abstract simplicial complex `K`.
"""
dim(K::SimplicialComplex) = pm_object(K).DIM::Int

@doc Markdown.doc"""
    f_vector(K::SimplicialComplex)

Return the face vector (number of faces per dimension) of the abstract simplicial complex `K`.

# Example
```jldoctest
julia> f_vector(torus())
3-element Vector{Int64}:
  7
 21
 14
```
"""
f_vector(K::SimplicialComplex) = Vector{Int}(pm_object(K).F_VECTOR)

@doc Markdown.doc"""
    h_vector(K::SimplicialComplex)

Return the h-vector of the abstract simplicial complex `K`.

# Example
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

@doc Markdown.doc"""
    betti_numbers(K::SimplicialComplex)

Return the reduced rational Betti numbers of the abstract simplicial complex `K`.

# Example
```jldoctest
julia> betti_numbers(klein_bottle())
3-element Vector{Int64}:
 0
 1
 0
```
"""
betti_numbers(K::SimplicialComplex) = Vector{Int}(Polymake.topaz.betti_numbers(pm_object(K)))

@doc Markdown.doc"""
    euler_characteristic(K::SimplicialComplex)

Return the reduced Euler characteristic of the abstract simplicial complex `K`.

# Example
```jldoctest
julia> euler_characteristic(complex_projective_plane())
2
```
"""
euler_characteristic(K::SimplicialComplex) = pm_object(K).EULER_CHARACTERISTIC::Int

@doc Markdown.doc"""
    homology(K::SimplicialComplex, i::Int)

Return `i`-th reduced integral homology group of `K`.
Recall that the 0-th homology group is trivial if and only if `K` is connected.

# Example
```jldoctest
julia> [ homology(real_projective_plane(), i) for i in [0,1,2] ]
3-element Vector{GrpAbFinGen}:
 GrpAb: Z/1
 GrpAb: Z/2
 GrpAb: Z/1
```
"""
homology(K::SimplicialComplex, i::Int) = abelian_group(_convert_finitely_generated_abelian_group(pm_object(K).HOMOLOGY[i+1])) # index shift

@doc Markdown.doc"""
    cohomology(K::SimplicialComplex, i::Int)

Return `i`-th reduced integral cohomology group of `K`.

# Example
```jldoctest
julia> K = SimplicialComplex([[0,1],[1,2],[0,2]]);

julia> cohomology(K,1)
(General) abelian group with relation matrix
[1]
```
"""
cohomology(K::SimplicialComplex, i::Int) = abelian_group(_convert_finitely_generated_abelian_group(pm_object(K).COHOMOLOGY[i+1])) # index shift

@doc Markdown.doc"""
    minimal_nonfaces(K::SimplicialComplex)

Return the minimal non-faces of the abstract simplicial complex `K`.

# Example
```jldoctest
julia> K = SimplicialComplex([[1,2,3],[2,3,4]]);

julia> minimal_nonfaces(K)
1-element Vector{Set{Int64}}:
 Set([4, 1])
```
"""
function minimal_nonfaces(K::SimplicialComplex)
    I = pm_object(K).MINIMAL_NON_FACES
    return Vector{Set{Int}}([Polymake.row(I,i) for i in 1:Polymake.nrows(I)])
end

@doc Markdown.doc"""
    alexander_dual(K::SimplicialComplex)

Return the Alexander dual of the abstract simplicial complex `K`.

# Example
```jldoctest
julia> K = SimplicialComplex([[1,2,3],[2,3,4]]);

julia> alexander_dual(K)
Abstract simplicial complex of dimension 1 on 2 vertices
```
"""
alexander_dual(K::SimplicialComplex) = SimplicialComplex(Polymake.topaz.alexander_dual(pm_object(K)))

@doc Markdown.doc"""
    stanley_reisner_ideal(K::SimplicialComplex)

Return the Stanley-Reisner ideal of the abstract simplicial complex `K`.

# Example
```jldoctest
julia> stanley_reisner_ideal(real_projective_plane())
ideal(x1*x2*x3, x1*x2*x4, x1*x5*x6, x2*x5*x6, x1*x3*x6, x1*x4*x5, x3*x4*x5, x3*x4*x6, x2*x3*x5, x2*x4*x6)
```
"""
function stanley_reisner_ideal(K::SimplicialComplex)
    n = nvertices(K)
    R, _ = PolynomialRing(QQ, n, cached=false)
    return stanley_reisner_ideal(R, K)
end

@doc Markdown.doc"""
    stanley_reisner_ideal(R::MPolyRing, K::SimplicialComplex)

Return the Stanley-Reisner ideal of the abstract simplicial complex `K`, in the given ring `R`.

# Example
```jldoctest
julia> R, _ = QQ["a","b","c","d","e","f"];

julia> stanley_reisner_ideal(R, real_projective_plane())
ideal(a*b*c, a*b*d, a*e*f, b*e*f, a*c*f, a*d*e, c*d*e, c*d*f, b*c*e, b*d*f)
```
"""
function stanley_reisner_ideal(R::MPolyRing, K::SimplicialComplex)
    n = nvertices(K)
    return ideal([ R([1], [_characteristic_vector(f,n)]) for f in minimal_nonfaces(K) ])
end

@doc Markdown.doc"""
    stanley_reisner_ring(K::SimplicialComplex)

Return the Stanley-Reisner ring of the abstract simplicial complex `K`.

# Example
```jldoctest
julia> K = SimplicialComplex([[1,2,3],[2,3,4]]);

julia> stanley_reisner_ring(K)
(Quotient of Multivariate Polynomial Ring in x1, x2, x3, x4 over Rational Field by ideal(x1*x4), Map from
Multivariate Polynomial Ring in x1, x2, x3, x4 over Rational Field to Quotient of Multivariate Polynomial Ring in x1, x2, x3, x4 over Rational Field by ideal(x1*x4) defined by a julia-function with inverse)
```
"""
function stanley_reisner_ring(K::SimplicialComplex)
    n = nvertices(K)
    R, _ = PolynomialRing(QQ, n, cached=false)
    return stanley_reisner_ring(R, K)
end

@doc Markdown.doc"""
    stanley_reisner_ring(R::MPolyRing, K::SimplicialComplex)

Return the Stanley-Reisner ring of the abstract simplicial complex `K`, as a quotient of a given ring `R`.

# Example
```jldoctest
julia>  R, _ = ZZ["a","b","c","d","e","f"];

julia> stanley_reisner_ring(R, real_projective_plane())
(Quotient of Multivariate Polynomial Ring in 6 variables a, b, c, d, ..., f over Integer Ring by ideal(a*b*c, a*b*d, a*e*f, b*e*f, a*c*f, a*d*e, c*d*e, c*d*f, b*c*e, b*d*f), Map from
Multivariate Polynomial Ring in 6 variables a, b, c, d, ..., f over Integer Ring to Quotient of Multivariate Polynomial Ring in 6 variables a, b, c, d, ..., f over Integer Ring by ideal(a*b*c, a*b*d, a*e*f, b*e*f, a*c*f, a*d*e, c*d*e, c*d*f, b*c*e, b*d*f) defined by a julia-function with inverse)
```
"""
stanley_reisner_ring(R::MPolyRing, K::SimplicialComplex) = quo(R, stanley_reisner_ideal(R, K))

################################################################################
###  Fundamental group
################################################################################

@doc Markdown.doc"""
    fundamental_group(K::SimplicialComplex)

Return the fundamental group of the abstract simplicial complex `K`.

# Example
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
                relem *= F[first(term) + 1]^last(term)
            end
            rvec[i] = relem
        end
        pi_1, _ = quo(F, rvec)
        return pi_1
    end
end

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

@doc Markdown.doc"""
    star_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}})

Return the star of the face `sigma` in the abstract simplicial complex `K`.

# Example
```jldoctest
julia> K = SimplicialComplex([[1,2,3],[2,3,4]]);

julia> star_subcomplex(K,[1])
Abstract simplicial complex of dimension 2 on 3 vertices
```
"""
star_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}}) = SimplicialComplex(Polymake.topaz.star_subcomplex(pm_object(K), Polymake.to_zero_based_indexing(sigma)))

@doc Markdown.doc"""
    link_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}})

Return the link of the face `sigma` in the abstract simplicial complex `K`.

# Example
```jldoctest
julia> K = SimplicialComplex([[1,2,3],[2,3,4]]);

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
    n = nvertices(K)
    print(io, "Abstract simplicial complex of dimension $(d) on $(n) vertices")
end

###############################################################################
## Serialization
###############################################################################

"""
    save_simplicialcomplex(K::SimplicialComplex, filename::String)

Save a SimplicialComplex to a file in JSON format.
"""
function save_simplicialcomplex(K::SimplicialComplex, filename::String)
    bigobject = pm_object(K)
    Polymake.save_bigobject(bigobject, filename)
end

"""
    load_simplicialcomplex(filename::String)

Load a SimplicialComplex stored in JSON format, given the filename as input.
"""
function load_simplicialcomplex(filename::String)
   bigobject = Polymake.load_bigobject(filename)
   typename = Polymake.type_name(bigobject)
   if typename[1:17] != "SimplicialComplex"
      throw(ArgumentError("Loaded object is not of type SimplicialComplex but rather " * typename))
   end
   return SimplicialComplex(bigobject)
end
