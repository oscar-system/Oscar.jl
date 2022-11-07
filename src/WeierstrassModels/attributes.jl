#######################################
# 1: Weierstrass sections
#######################################

@doc Markdown.doc"""
    poly_f(w::GlobalWeierstrassModel)

Return the polynomial ``f`` used for the
construction of the global Weierstrass model.

```jldoctest
julia> using Oscar

julia> test_space = hirzebruch_surface(2) * projective_space(NormalToricVariety,1)
A normal toric variety

julia> test_space1 = blowup_on_ith_minimal_torus_orbit(test_space,1,"e1")
A normal toric variety

julia> test_space2 = blowup_on_ith_minimal_torus_orbit(test_space1,1,"e2")
A normal toric variety

julia> base = blowup_on_ith_minimal_torus_orbit(test_space2,1,"e3")
A normal toric variety

julia> w = GenericGlobalWeierstrassModel(base)
A global Weierstrass model

julia> poly_f(w);
```
"""
@attr MPolyElem{fmpq} poly_f(w::GlobalWeierstrassModel) = w.poly_f
export poly_f


@doc Markdown.doc"""
    poly_g(w::GlobalWeierstrassModel)

Return the polynomial ``g`` used for the
construction of the global Weierstrass model.

```jldoctest
julia> using Oscar

julia> test_space = hirzebruch_surface(2) * projective_space(NormalToricVariety,1)
A normal toric variety

julia> test_space1 = blowup_on_ith_minimal_torus_orbit(test_space,1,"e1")
A normal toric variety

julia> test_space2 = blowup_on_ith_minimal_torus_orbit(test_space1,1,"e2")
A normal toric variety

julia> base = blowup_on_ith_minimal_torus_orbit(test_space2,1,"e3")
A normal toric variety

julia> w = GenericGlobalWeierstrassModel(base)
A global Weierstrass model

julia> poly_g(w);
```
"""
@attr MPolyElem{fmpq} poly_g(w::GlobalWeierstrassModel) = w.poly_g
export poly_g


#######################################
# 2: Weierstrass polynomial
#######################################

@doc Markdown.doc"""
    pw(w::GlobalWeierstrassModel)

Return the Weierstrass polynomial of the global Weierstrass model.

```jldoctest
julia> using Oscar

julia> test_space = hirzebruch_surface(2) * projective_space(NormalToricVariety,1)
A normal toric variety

julia> test_space1 = blowup_on_ith_minimal_torus_orbit(test_space,1,"e1")
A normal toric variety

julia> test_space2 = blowup_on_ith_minimal_torus_orbit(test_space1,1,"e2")
A normal toric variety

julia> base = blowup_on_ith_minimal_torus_orbit(test_space2,1,"e3")
A normal toric variety

julia> w = GenericGlobalWeierstrassModel(base)
A global Weierstrass model

julia> pw(w);
```
"""
@attr MPolyElem{fmpq} pw(w::GlobalWeierstrassModel) = w.pw
export pw


#######################################
# 3: Toric spaces
#######################################

@doc Markdown.doc"""
    toric_base_space(w::GlobalWeierstrassModel)

Return the toric base space of the global Weierstrass model.

```jldoctest
julia> using Oscar

julia> test_space = hirzebruch_surface(2) * projective_space(NormalToricVariety,1)
A normal toric variety

julia> test_space1 = blowup_on_ith_minimal_torus_orbit(test_space,1,"e1")
A normal toric variety

julia> test_space2 = blowup_on_ith_minimal_torus_orbit(test_space1,1,"e2")
A normal toric variety

julia> base = blowup_on_ith_minimal_torus_orbit(test_space2,1,"e3")
A normal toric variety

julia> w = GenericGlobalWeierstrassModel(base)
A global Weierstrass model

julia> is_smooth(toric_base_space(w))
true
```
"""
@attr Oscar.AbstractNormalToricVariety toric_base_space(w::GlobalWeierstrassModel) = w.toric_base_space
export toric_base_space


@doc Markdown.doc"""
    toric_ambient_space(w::GlobalWeierstrassModel)

Return the toric base space of the global Weierstrass model.

```jldoctest
julia> using Oscar

julia> test_space = hirzebruch_surface(2) * projective_space(NormalToricVariety,1)
A normal toric variety

julia> test_space1 = blowup_on_ith_minimal_torus_orbit(test_space,1,"e1")
A normal toric variety

julia> test_space2 = blowup_on_ith_minimal_torus_orbit(test_space1,1,"e2")
A normal toric variety

julia> base = blowup_on_ith_minimal_torus_orbit(test_space2,1,"e3")
A normal toric variety

julia> w = GenericGlobalWeierstrassModel(base)
A global Weierstrass model

julia> is_smooth(toric_ambient_space(w))
false
```
"""
@attr Oscar.AbstractNormalToricVariety toric_ambient_space(w::GlobalWeierstrassModel) = w.toric_ambient_space
export toric_base_space
