using Markdown
import Oscar: Polymake, pm_object

export
    SimplicialComplex,
    bettinumbers,
    dim,
    eulercharacteristic,
    f_vector,
    h_vector,
    minimalnonfaces,
    nvertices,
    load_simplicialcomplex,
    save_simplicialcomplex,
    complexprojectiveplane,
    realprojectiveplane,
    kleinbottle,
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
    chi = Vector{Int}(fill(0,n))
    for x in M
        chi[x] = 1
    end
    return chi
end

################################################################################
##  Properties
################################################################################

dim(K::SimplicialComplex) = pm_object(K).DIM
nvertices(K::SimplicialComplex) = pm_object(K).N_VERTICES
f_vector(K::SimplicialComplex) = Vector{Int}(pm_object(K).F_VECTOR)
h_vector(K::SimplicialComplex) = Vector{Int}(pm_object(K).H_VECTOR)
bettinumbers(K::SimplicialComplex) = Polymake.topaz.betti_numbers(pm_object(K))
eulercharacteristic(K::SimplicialComplex) = pm_object(K).EULER_CHARACTERISTIC

@doc Markdown.doc"""
    minimalnonfaces(K::SimplicialComplex)

Compute the minimal non-faces of K.

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

################################################################################
##  Standard examples
################################################################################

torus() = SimplicialComplex(Polymake.topaz.torus())
kleinbottle() = SimplicialComplex(Polymake.topaz.klein_bottle())
realprojectiveplane() = SimplicialComplex(Polymake.topaz.real_projective_plane())
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
