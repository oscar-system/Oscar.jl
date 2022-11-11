#######################################
# 1: Weierstrass sections
#######################################

@doc Markdown.doc"""
    weierstrass_section_f(w::GlobalWeierstrassModel)

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

julia> w = GlobalWeierstrassModel(base)
A global Weierstrass model over a concrete base

julia> weierstrass_section_f(w);
```
"""
@attr MPolyElem{fmpq} weierstrass_section_f(w::GlobalWeierstrassModel) = w.poly_f
export weierstrass_section_f


@doc Markdown.doc"""
    weierstrass_section_g(w::GlobalWeierstrassModel)

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

julia> w = GlobalWeierstrassModel(base)
A global Weierstrass model over a concrete base

julia> weierstrass_section_g(w);
```
"""
@attr MPolyElem{fmpq} weierstrass_section_g(w::GlobalWeierstrassModel) = w.poly_g
export weierstrass_section_g


#######################################
# 2: Weierstrass polynomial
#######################################

@doc Markdown.doc"""
    weierstrass_polynomial(w::GlobalWeierstrassModel)

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

julia> w = GlobalWeierstrassModel(base)
A global Weierstrass model over a concrete base

julia> weierstrass_polynomial(w);
```
"""
@attr MPolyElem{fmpq} weierstrass_polynomial(w::GlobalWeierstrassModel) = w.pw
export weierstrass_polynomial


#######################################
# 3: Discriminant
#######################################

@doc Markdown.doc"""
    discriminant(w::GlobalWeierstrassModel)

Return the discriminant ``\Delta = 4 f^3 + 27 g^2``.

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

julia> w = GlobalWeierstrassModel(base)
A global Weierstrass model over a concrete base

julia> discriminant(w);
```
"""
@attr MPolyElem{fmpq} Oscar.:discriminant(w::GlobalWeierstrassModel) = 4 * w.poly_f^3 + 27 * w.poly_g^2
export discriminant


#####################################################
# 4: The CY hypersurface
#####################################################

@doc Markdown.doc"""
    cy_hypersurface(w::GlobalWeierstrassModel)

Return the Calabi-Yau hypersurface in the toric ambient space
which defines the global Weierstrass model.

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

julia> w = GlobalWeierstrassModel(base)
A global Weierstrass model over a concrete base

julia> cy_hypersurface(w)
A closed subvariety of a normal toric variety
```
"""
@attr Oscar.ClosedSubvarietyOfToricVariety cy_hypersurface(w::GlobalWeierstrassModel) = w.Y4
export cy_hypersurface


#######################################
# 5: Toric spaces
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

julia> w = GlobalWeierstrassModel(base)
A global Weierstrass model over a concrete base

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

julia> w = GlobalWeierstrassModel(base)
A global Weierstrass model over a concrete base

julia> is_smooth(toric_ambient_space(w))
false
```
"""
@attr Oscar.AbstractNormalToricVariety toric_ambient_space(w::GlobalWeierstrassModel) = w.toric_ambient_space
export toric_base_space
