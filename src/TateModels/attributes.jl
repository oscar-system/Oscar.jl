#####################################################
# 1: The Tate sections
#####################################################

@doc Markdown.doc"""
    tate_section_a1(t::GlobalTateModel)

Return the Tate section ``a_1``.

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

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> tate_section_a1(t);
```
"""
@attr MPolyElem{fmpq} tate_section_a1(t::GlobalTateModel) = t.a1
export tate_section_a1


@doc Markdown.doc"""
    tate_section_a2(t::GlobalTateModel)

Return the Tate section ``a_2``.

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

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> tate_section_a2(t);
```
"""
@attr MPolyElem{fmpq} tate_section_a2(t::GlobalTateModel) = t.a2
export tate_section_a2


@doc Markdown.doc"""
    tate_section_a3(t::GlobalTateModel)

Return the Tate section ``a_3``.

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

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> tate_section_a3(t);
```
"""
@attr MPolyElem{fmpq} tate_section_a3(t::GlobalTateModel) = t.a3
export tate_section_a3


@doc Markdown.doc"""
    tate_section_a4(t::GlobalTateModel)

Return the Tate section ``a_4``.

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

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> tate_section_a4(t);
```
"""
@attr MPolyElem{fmpq} tate_section_a4(t::GlobalTateModel) = t.a4
export tate_section_a4


@doc Markdown.doc"""
    tate_section_a6(t::GlobalTateModel)

Return the Tate section ``a_6``.

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

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> tate_section_a6(t);
```
"""
@attr MPolyElem{fmpq} tate_section_a6(t::GlobalTateModel) = t.a6
export tate_section_a6


#####################################################
# 2: The Tate polynomial
#####################################################

@doc Markdown.doc"""
    tate_polynomial(t::GlobalTateModel)

Return the Tate polynomial of the global Tate model.

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

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> tate_polynomial(t);
```
"""
@attr MPolyElem{fmpq} tate_polynomial(t::GlobalTateModel) = t.pt
export tate_polynomial


#####################################################
# 3: Toric spaces for Tate models over concrete bases
#####################################################

@doc Markdown.doc"""
    toric_base_space(t::GlobalTateModel)

Return the toric base space of the global Tate model.

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

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> toric_base_space(t)
A normal toric variety without torusfactor
```
"""
@attr Oscar.AbstractNormalToricVariety toric_base_space(t::GlobalTateModel) = t.toric_base_space
export toric_base_space


@doc Markdown.doc"""
    toric_ambient_space(t::GlobalTateModel)

Return the toric ambient space of the global Tate model.

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

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> toric_ambient_space(t)
A normal, simplicial toric variety

julia> is_smooth(toric_ambient_space(t))
false
```
"""
@attr Oscar.AbstractNormalToricVariety toric_ambient_space(t::GlobalTateModel) = t.toric_ambient_space
export toric_ambient_space


#####################################################
# 4: The CY hypersurface
#####################################################

@doc Markdown.doc"""
    cy_hypersurface(t::GlobalTateModel)

Return the Calabi-Yau hypersurface in the toric ambient space
which defines the global Tate model.

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

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> cy_hypersurface(t)
A closed subvariety of a normal toric variety
```
"""
@attr Oscar.ClosedSubvarietyOfToricVariety cy_hypersurface(t::GlobalTateModel) = t.Y4
export cy_hypersurface


#####################################################
# 5: Turn global Tate model into a Weierstrass model
#####################################################

@doc Markdown.doc"""
    global_weierstrass_model(t::GlobalTateModel)

Return the global Weierstrass model which is equivalent to the given Tate model.

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

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> global_weierstrass_model(t)
A global Weierstrass model over a concrete base
```
"""
@attr GlobalWeierstrassModel function global_weierstrass_model(t::GlobalTateModel)
    b2 = 4 * tate_section_a2(t) + tate_section_a1(t)^2
    b4 = 2 * tate_section_a4(t) + tate_section_a1(t) * tate_section_a3(t)
    b6 = 4 * tate_section_a6(t) + tate_section_a3(t)^2
    f = - 1//48 * (b2^2 - 24 * b4)
    g = 1//864 * (b2^3 - 36 * b2 * b4 + 216 * b6)
    S = cox_ring(toric_ambient_space(t))
    x = gens(S)[length(gens(S))-2]
    y = gens(S)[length(gens(S))-1]
    z = gens(S)[length(gens(S))]
    pw = x^3 - y^2 + f*x*z^4 + g*z^6
    Y4 = Oscar.ClosedSubvarietyOfToricVariety(toric_ambient_space(t), [pw])
    model = GlobalWeierstrassModel(f, g, pw, toric_base_space(t), toric_ambient_space(t), Y4)
    set_attribute!(model, :base_fully_specified, base_fully_specified(t))
    return model
end
export global_weierstrass_model


#####################################################
# 6: Discriminant
#####################################################

@doc Markdown.doc"""
    discriminant(t::GlobalTateModel)

Return the discriminant of the global Tate model.

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

julia> t = GlobalTateModel(base)
A global Tate model over a concrete base

julia> discriminant(t);
```
"""
@attr MPolyElem{fmpq} Oscar.:discriminant(t::GlobalTateModel) = discriminant(global_weierstrass_model(t))
export discriminant
