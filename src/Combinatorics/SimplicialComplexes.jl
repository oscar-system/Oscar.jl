using Markdown
import Oscar: Polymake, pm_object

export
    SimplicialComplex

################################################################################
################################################################################
##  Constructing
################################################################################
################################################################################
struct SimplicialComplex
    pm_complex::Polymake.BigObject
end

function pm_object(C::SimplicialComplex)
    return K.pm_complex
end

@doc Markdown.doc"""
    SimplicialComplex(generators::Vector{Vector{Int64}})

Construct an abstract simplicial complex from a set of faces.

# Example
```jldoctest
julia> K = SimplicialComplex([[1,2,3],[2,3,4]]);
```
"""
function SimplicialComplex(generators::Vector{Vector{Int64}})
    K = Polymake.topaz.SimplicialComplex(INPUT_FACES=generators)
    SimplicialComplex(K)
end
