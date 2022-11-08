#######################################
# 1: The Tate sections
#######################################

@doc Markdown.doc"""
    a1(t::GlobalTateModel)

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

julia> a1(t);
```
"""
@attr MPolyElem{fmpq} a1(t::GlobalTateModel) = t.a1
export a1


@doc Markdown.doc"""
    a2(t::GlobalTateModel)

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

julia> a2(t);
```
"""
@attr MPolyElem{fmpq} a2(t::GlobalTateModel) = t.a2
export a2


@doc Markdown.doc"""
    a3(t::GlobalTateModel)

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

julia> a3(t);
```
"""
@attr MPolyElem{fmpq} a3(t::GlobalTateModel) = t.a3
export a3


@doc Markdown.doc"""
    a4(t::GlobalTateModel)

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

julia> a4(t);
```
"""
@attr MPolyElem{fmpq} a4(t::GlobalTateModel) = t.a4
export a4


@doc Markdown.doc"""
    a6(t::GlobalTateModel)

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

julia> a6(t);
```
"""
@attr MPolyElem{fmpq} a6(t::GlobalTateModel) = t.a6
export a6


#######################################
# 2: The Tate polynomial
#######################################

@doc Markdown.doc"""
    pt(t::GlobalTateModel)

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

julia> pt(t);
```
"""
@attr MPolyElem{fmpq} pt(t::GlobalTateModel) = t.pt
export pt


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
