#######################################
# 1: The Tate sections
#######################################

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

julia> t = GenericGlobalTateModel(base)
A global Tate model

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

julia> t = GenericGlobalTateModel(base)
A global Tate model

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

julia> t = GenericGlobalTateModel(base)
A global Tate model

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

julia> t = GenericGlobalTateModel(base)
A global Tate model

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

julia> t = GenericGlobalTateModel(base)
A global Tate model

julia> tate_section_a6(t);
```
"""
@attr MPolyElem{fmpq} tate_section_a6(t::GlobalTateModel) = t.a6
export tate_section_a6


#######################################
# 2: The Tate polynomial
#######################################

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

julia> t = GenericGlobalTateModel(base)
A global Tate model

julia> tate_polynomial(t);
```
"""
@attr MPolyElem{fmpq} tate_polynomial(t::GlobalTateModel) = t.pt
export tate_polynomial


#######################################
# 3: Toric spaces
#######################################

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

julia> t = GenericGlobalTateModel(base)
A global Tate model

julia> toric_base_space(t)
A normal, 3-dimensional toric variety without torusfactor
```
"""
@attr Oscar.AbstractNormalToricVariety toric_base_space(t::GlobalTateModel) = t.base
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

julia> t = GenericGlobalTateModel(base)
A global Tate model

julia> toric_ambient_space(t)
A normal toric variety

julia> is_smooth(toric_ambient_space(t))
false
```
"""
@attr Oscar.AbstractNormalToricVariety toric_ambient_space(t::GlobalTateModel) = t.toric_ambient_space
export toric_ambient_space
