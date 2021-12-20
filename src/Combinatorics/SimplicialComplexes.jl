using Markdown
import Oscar: Polymake, pm_object

export
    SimplicialComplex,
    bettinumbers,
    dim,
    euler_characteristic,
    f_vector,
    h_vector,
    nvertices,
    load_simplicialcomplex,
    save_simplicialcomplex,
    torus, # requires a distinction from, e.g., an algebraic group
    kleinbottle

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
function SimplicialComplex(generators::Vector{Vector{Int}})
    K = Polymake.topaz.SimplicialComplex(INPUT_FACES=generators)
    SimplicialComplex(K)
end


################################################################################
##  Properties
################################################################################

function dim(K::SimplicialComplex)
    return pm_object(K).DIM
end

function bettinumbers(K::SimplicialComplex)
    return Polymake.topaz.betti_numbers(pm_object(K))
end

function euler_characteristic(K::SimplicialComplex)
    return pm_object(K).EULER_CHARACTERISTIC
end

function f_vector(K::SimplicialComplex)
    return Vector{Int}(pm_object(K).F_VECTOR)
end

function h_vector(K::SimplicialComplex)
    return Vector{Int}(pm_object(K).H_VECTOR)
end

function nvertices(K::SimplicialComplex)
    return pm_object(K).N_VERTICES
end

################################################################################
##  Standard examples
################################################################################

function torus()
    return SimplicialComplex(Polymake.topaz.torus())
end

function kleinbottle()
    return SimplicialComplex(Polymake.topaz.klein_bottle())
end

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
