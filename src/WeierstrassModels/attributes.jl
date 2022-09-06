@doc Markdown.doc"""
    poly_f(w::GlobalWeierstrassModel)

Return the polynomial $f$ used for the
construction of the global Weierstrass model.

```jldoctest
julia> w = GenericGlobalWeierstrassModelOverProjectiveSpace(3)
A global Weierstrass model

julia> poly_f(w);
```
"""
@attr MPolyElem{fmpq} poly_f(w::GlobalWeierstrassModel) = w.poly_f
export poly_f

@doc Markdown.doc"""
    poly_g(w::GlobalWeierstrassModel)

Return the polynomial $g$ used for the
construction of the global Weierstrass model.

```jldoctest
julia> w = GenericGlobalWeierstrassModelOverProjectiveSpace(3)
A global Weierstrass model

julia> poly_f(w);
```
"""
@attr MPolyElem{fmpq} poly_g(w::GlobalWeierstrassModel) = w.poly_g
export poly_g


@doc Markdown.doc"""
    toric_base_space(w::GlobalWeierstrassModel)

If the global Weierstrass model in question was constructed
over a toric base space, this method returns this toric space.
Otherwise an error is raised.

```jldoctest
julia> w = GenericGlobalWeierstrassModelOverProjectiveSpace(3)
A global Weierstrass model

julia> toric_base_space(w)
A normal, non-affine, smooth, projective, gorenstein, fano, 3-dimensional toric variety without torusfactor
```
"""
function toric_base_space(w::GlobalWeierstrassModel)::Oscar.AbstractNormalToricVariety
    if has_attribute(w, :toric_base_space)
       return get_attribute(w, :toric_base_space)
    else
        throw(ArgumentError("The Weierstrass model in question was not defined over a toric base space"))
    end
end
export toric_base_space
