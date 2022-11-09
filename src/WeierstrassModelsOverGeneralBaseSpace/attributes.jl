#######################################
# 1: The Tate sections
#######################################

@doc Markdown.doc"""
    weierstrass_section_f(w::WeierstrassModelOverGeneralBaseSpace)

Return the Weierstrass section ``f``.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (f, g) = QQ["f", "g"];

julia> w = GlobalWeierstrassModel(f, g, auxiliary_base_ring)
A global Weierstrass model over a general base space

julia> weierstrass_section_f(w)
f
```
"""
@attr MPolyElem{fmpq} weierstrass_section_f(w::WeierstrassModelOverGeneralBaseSpace) = w.f
export weierstrass_section_f


@doc Markdown.doc"""
    weierstrass_section_g(w::WeierstrassModelOverGeneralBaseSpace)

Return the Weierstrass section ``g``.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (f, g) = QQ["f", "g"];

julia> w = GlobalWeierstrassModel(f, g, auxiliary_base_ring)
A global Weierstrass model over a general base space

julia> weierstrass_section_g(w)
g
```
"""
@attr MPolyElem{fmpq} weierstrass_section_g(w::WeierstrassModelOverGeneralBaseSpace) = w.g
export weierstrass_section_g


#######################################
# 2: The Weierstrass polynomial
#######################################

@doc Markdown.doc"""
    weierstrass_polynomial(w::WeierstrassModelOverGeneralBaseSpace)

Return the Weierstrass polynomial.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (f, g) = QQ["f", "g"];

julia> w = GlobalWeierstrassModel(f, g, auxiliary_base_ring)
A global Weierstrass model over a general base space

julia> weierstrass_polynomial(w)
f*x*z^4 + g*z^6 + x^3 - y^2
```
"""
@attr MPolyElem{fmpq} weierstrass_polynomial(w::WeierstrassModelOverGeneralBaseSpace) = w.pw
export weierstrass_polynomial


#######################################
# 3: Discriminant
#######################################

@doc Markdown.doc"""
    discriminant(w::WeierstrassModelOverGeneralBaseSpace)

Return the discriminant ``\Delta = 4 f^3 + 27 g^2``.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (f, g) = QQ["f", "g"];

julia> w = GlobalWeierstrassModel(f, g, auxiliary_base_ring)
A global Weierstrass model over a general base space

julia> discriminant(w);
```
"""
@attr MPolyElem{fmpq} Oscar.:discriminant(w::WeierstrassModelOverGeneralBaseSpace) = 4 * weierstrass_section_f(w)^3 + 27 * weierstrass_section_g(w)^2
export discriminant


#######################################
# 4: Auxiliary toric spaces
#######################################

@doc Markdown.doc"""
    auxiliary_base_space(w::WeierstrassModelOverGeneralBaseSpace)

Return the toric base space of the global Tate model.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (f, g) = QQ["f", "g"];

julia> w = GlobalWeierstrassModel(f, g, auxiliary_base_ring)
A global Weierstrass model over a general base space

julia> auxiliary_base_space(w)
A normal, affine, 2-dimensional toric variety
```
"""
@attr Oscar.AbstractNormalToricVariety auxiliary_base_space(w::WeierstrassModelOverGeneralBaseSpace) = w.auxiliary_base_space
export auxiliary_base_space


@doc Markdown.doc"""
    auxiliary_ambient_space(w::WeierstrassModelOverGeneralBaseSpace)

Return the toric ambient space of the global Tate model.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (f, g) = QQ["f", "g"];

julia> w = GlobalWeierstrassModel(f, g, auxiliary_base_ring)
A global Weierstrass model over a general base space

julia> auxiliary_ambient_space(w)
A normal toric variety
```
"""
@attr Oscar.AbstractNormalToricVariety auxiliary_ambient_space(w::WeierstrassModelOverGeneralBaseSpace) = w.auxiliary_ambient_space
export auxiliary_ambient_space
