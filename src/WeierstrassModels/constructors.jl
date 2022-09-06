##############################################
# 1: The Julia type for GlobalWeierstrassModel
##############################################

@attributes mutable struct GlobalWeierstrassModel
    poly_f::MPolyElem{fmpq}
    poly_g::MPolyElem{fmpq}
    function GlobalWeierstrassModel(poly_f::MPolyElem{fmpq}, poly_g::MPolyElem{fmpq})
        if parent(poly_f) != parent(poly_g)
            throw(ArgumentError("The two polynomials f, g must be elements of the same ring, i.e. a ring encoded by the same julia variable"))
        end
        return new(poly_f, poly_g)
    end
end
export GlobalWeierstrassModel


#######################################
# 2: Generic constructor
#######################################

@doc Markdown.doc"""
    GlobalWeierstrassModel(polys::Vector{MPolyElem{fmpq}})

A global Weierstrass model is a hypersurface cut out by an
equation of the form $P_w = y^2 - x^3 - f x z^4 - g z^6$.
Consequently, it is specified by two polynomials $f$ and $g$.
For convenience, it can be constructed as
`GlobalWeierstrassModel([f,g])` and via `GlobalWeierstrassModel(f,g)`.
"""
function GlobalWeierstrassModel(polys::Vector{MPolyElem{fmpq}})
    if length(polys) != 2
        throw(ArgumentError("Exactly two polynomials f and g must be specified to define a Weierstrass model"))
    end
    return GlobalWeierstrassModel(polys[1],polys[2])
end
export GlobalWeierstrassModel


@doc Markdown.doc"""
    GenericGlobalWeierstrassModelOverToricSpace(v::Oscar.AbstractNormalToricVariety)

Given a toric variety $v$, this method determines
generic sections $f \in H^0( v, \overline{K}_v^4)$
and $g \in H^0( v, \overline{K}_v^6)$. Subsequently,
it constructs the global Weierstrass model defined by
those two polynomials.
"""
function GenericGlobalWeierstrassModelOverToricSpace(v::Oscar.AbstractNormalToricVariety)
    Kbar = anticanonical_bundle(v)
    f = sum([rand(Int)*b for b in basis_of_global_sections(Kbar^4)])
    g = sum([rand(Int)*b for b in basis_of_global_sections(Kbar^6)])
    w = GlobalWeierstrassModel(f,g)
    set_attribute!(w, :toric_base_space, v)
    return w
end
export GenericGlobalWeierstrassModelOverToricSpace


@doc Markdown.doc"""
    GenericGlobalWeierstrassModelOverProjectiveSpace(n::Int)

This method constructs the $n$-dimensional projective
space, determines generic sections $f$ and $g$ of the
4-th and 6-th power of the anticanonical bundle and then
constructs the global Weierstrass model defined by
those two polynomials.

```jldoctest
julia> w = GenericGlobalWeierstrassModelOverProjectiveSpace(3)
A global Weierstrass model
```
"""
function GenericGlobalWeierstrassModelOverProjectiveSpace(n::Int)
    if n <= 0
        throw(ArgumentError("The integer specifies the dimension of the projective space and must be at least 1"))
    end
    return GenericGlobalWeierstrassModelOverToricSpace(projective_space(NormalToricVariety,n))
end
export GenericGlobalWeierstrassModelOverProjectiveSpace


#######################################
# 3: Display
#######################################

function Base.show(io::IO, cy::GlobalWeierstrassModel)
    join(io, "A global Weierstrass model")
end
