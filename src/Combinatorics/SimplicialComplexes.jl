using Markdown
import Oscar: Polymake, pm_object

export
    SimplicialComplex,
    dim,
    euler_characteristic,
    f_vector,
    h_vector,
    nvertices

################################################################################
##  Constructing
################################################################################

struct SimplicialComplex
    pm_complex::Polymake.BigObject
end

function pm_object(C::SimplicialComplex)
    return K.pm_complex
end

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
    return K.pm_complex.DIM
end

function euler_characteristic(K::SimplicialComplex)
    return K.pm_complex.EULER_CHARACTERISTIC
end

function f_vector(K::SimplicialComplex)
    return Vector{Int}(K.pm_complex.F_VECTOR)
end

function h_vector(K::SimplicialComplex)
    return Vector{Int}(K.pm_complex.H_VECTOR)
end

function nvertices(K::SimplicialComplex)
    return K.pm_complex.N_VERTICES
end

################################################################################
##  Standard examples
################################################################################

function torus()
    return SimplicialComplex(Polymake.topaz.torus())
end

function klein_bottle()
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
