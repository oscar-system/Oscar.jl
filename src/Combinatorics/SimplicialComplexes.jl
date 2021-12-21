using Markdown
import Oscar: Polymake, pm_object

export
    SimplicialComplex,
    betti_numbers,
    dim,
    euler_characteristic,
    f_vector,
    h_vector,
    minimalnonfaces,
    stanley_reisner_ideal,
    stanley_reisner_ring,
    nvertices,
    load_simplicialcomplex,
    save_simplicialcomplex,
    complexprojectiveplane,
    realprojectiveplane,
    klein_bottle,
    torus # requires a distinction from, e.g., an algebraic group

################################################################################
##  Constructing
################################################################################

struct SimplicialComplex
    pm_simplicialcomplex::Polymake.BigObject
end

pm_object(K::SimplicialComplex) = K.pm_simplicialcomplex


@doc Markdown.doc"""
    SimplicialComplex(generators::Vector{Vector{Int}})

Construct an abstract simplicial complex from a set of faces.

# Example
```jldoctest
julia> K = SimplicialComplex([[1,2,3],[2,3,4]])
Abstract simplicial complex of dimension 2 on 4 vertices
```
"""
function SimplicialComplex(generators::Vector{Set{Int}})
    K = Polymake.topaz.SimplicialComplex(INPUT_FACES=generators)
    SimplicialComplex(K)
end

function SimplicialComplex(generators::Vector{Vector{Int}})
    K = Polymake.topaz.SimplicialComplex(INPUT_FACES=generators)
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

_reindexset(M::Set{Int},ind::Vector{Int}) = [ ind[x+1] for x in M ]

function _characteristicvector(M::Vector{Int},n::Int)
    chi = zeros(Int, n)
    chi[M] .= 1
    return chi
end

################################################################################
##  Properties
################################################################################

@doc Markdown.doc"""
    dim(SimplicialComplex)

Dimension of an abstract simplicial complex.
"""
dim(K::SimplicialComplex) = pm_object(K).DIM

@doc Markdown.doc"""
    nvertices(SimplicialComplex)

Number of vertices of an abstract simplicial complex.
"""
nvertices(K::SimplicialComplex) = pm_object(K).N_VERTICES

@doc Markdown.doc"""
    f_vector(SimplicialComplex)

Face vector (number of faces per dimension) of an abstract simplicial complex.
"""
f_vector(K::SimplicialComplex) = Vector{Int}(pm_object(K).F_VECTOR)

@doc Markdown.doc"""
    f_vector(SimplicialComplex)

H-vector of an abstract simplicial complex.
"""
h_vector(K::SimplicialComplex) = Vector{Int}(pm_object(K).H_VECTOR)

@doc Markdown.doc"""
    betti_numbers(SimplicialComplex)

Rational Betti numbers of an abstract simplicial complex.
"""
betti_numbers(K::SimplicialComplex) = Polymake.topaz.betti_numbers(pm_object(K))

@doc Markdown.doc"""
    euler_characteristic(SimplicialComplex)

Euler characteristic of an abstract simplicial complex.
"""
euler_characteristic(K::SimplicialComplex) = pm_object(K).EULER_CHARACTERISTIC

@doc Markdown.doc"""
    minimalnonfaces(SimplicialComplex)

Minimal non-faces of an abstract simplicial complex.

# Example
```jldoctest
julia> K = SimplicialComplex([[1,2,3],[2,3,4]]);

julia> minimalnonfaces(K)
1-element Vector{Vector{Int64}}:
 [1, 4]
```
"""
function minimalnonfaces(K::SimplicialComplex)
    bigobject = pm_object(K)
    mnf = Vector{Set{Int}}(bigobject.MINIMAL_NON_FACES)
    ind = _vertexindices(bigobject)
    return [ _reindexset(nonface,ind) for nonface in mnf ]
end

@doc Markdown.doc"""
    stanley_reisner_ideal(SimplicialComplex)

Stanley-Reisner ideal of an abstract simplicial complex.

# Example
```jldoctest
julia> stanley_reisner_ideal(realprojectiveplane())
ideal(x1*x2*x3, x1*x2*x4, x1*x5*x6, x2*x5*x6, x1*x3*x6, x1*x4*x5, x3*x4*x5, x3*x4*x6, x2*x3*x5, x2*x4*x6)
```
"""
function stanley_reisner_ideal(K::SimplicialComplex)
    n = nvertices(K)
    R, () = PolynomialRing(ZZ, n)
    return ideal([ R([ZZ(1)], [_characteristicvector(f,n)]) for f in minimalnonfaces(K) ])
end

@doc Markdown.doc"""
    stanley_reisner_ring(SimplicialComplex)

Stanley-Reisner ring of an abstract simplicial complex.

# Example
```jldoctest
julia> K = SimplicialComplex([[1,2,3],[2,3,4]]);

julia> stanley_reisner_ring(K)
(Quotient of Multivariate Polynomial Ring in x1, x2, x3, x4 over Integer Ring by ideal(x1*x4), Map from
Multivariate Polynomial Ring in x1, x2, x3, x4 over Integer Ring to Quotient of Multivariate Polynomial Ring in x1, x2, x3, x4 over Integer Ring by ideal(x1*x4) defined by a julia-function with inverse)
"""
function stanley_reisner_ring(K::SimplicialComplex)
    I = stanley_reisner_ideal(K)
    return quo(base_ring(I),I)
end


################################################################################
##  Standard examples
################################################################################

@doc Markdown.doc"""
    torus()

Császár's 7-vertex triangulation of the torus (surface).
"""
torus() = SimplicialComplex(Polymake.topaz.torus())

@doc Markdown.doc"""
    klein_bottle()

9-vertex triangulation of the Klein bottle.
"""
klein_bottle() = SimplicialComplex(Polymake.topaz.klein_bottle())

@doc Markdown.doc"""
    realprojectiveplane()

6-vertex triangulation of the real projective plane.
"""
realprojectiveplane() = SimplicialComplex(Polymake.topaz.real_projective_plane())

@doc Markdown.doc"""
    complexprojectiveplane()

9-vertex triangulation of the complex projective plane.
"""
complexprojectiveplane() = SimplicialComplex(Polymake.topaz.complex_projective_plane())

###############################################################################
### Display
###############################################################################

function Base.show(io::IO, K::SimplicialComplex)
    d = dim(K)
    n = nvertices(K)
    print(io, "Abstract simplicial complex of dimension $(d) on $(n) vertices")
end

###############################################################################
### Serialization
###############################################################################

"""
    save_simplicialcomplex(SimplicialComplex, String)

Save a SimplicialComplex to a file in JSON format.
"""
function save_simplicialcomplex(K::SimplicialComplex, filename::String)
    bigobject = pm_object(K)
    Polymake.save_bigobject(bigobject, filename)
end

"""
    load_simplicialcomplex(String)

Load a SimplicialComplex stored in JSON format, given the filename as input.
"""
function load_simplicialcomplex(filename::String)
   bigobject = Polymake.load_bigobject(filename)
   typename = Polymake.type_name(bigobject)
   if typename[1:4] != "SimplicialComplex"
      throw(ArgumentError("Loaded object is not of type SimplicialComplex but rather " * typename))
   end
   return SimplicialComplex(bigobject)
end
