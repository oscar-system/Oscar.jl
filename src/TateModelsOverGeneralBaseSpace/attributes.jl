#######################################
# 1: The Tate sections
#######################################

@doc Markdown.doc"""
    tate_section_a1(t::TateModelOverGeneralBaseSpace)

Return the Tate section ``a_1``.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a10,a21,a32,a43,a65,w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = GlobalTateModel(ais, auxiliary_base_ring)
A global Tate model over a general base space

julia> tate_section_a1(t)
a10
```
"""
@attr MPolyElem{fmpq} tate_section_a1(t::TateModelOverGeneralBaseSpace) = t.a1
export tate_section_a1


@doc Markdown.doc"""
    tate_section_a2(t::TateModelOverGeneralBaseSpace)

Return the Tate section ``a_2``.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a10,a21,a32,a43,a65,w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = GlobalTateModel(ais, auxiliary_base_ring)
A global Tate model over a general base space

julia> tate_section_a2(t)
a21*w
```
"""
@attr MPolyElem{fmpq} tate_section_a2(t::TateModelOverGeneralBaseSpace) = t.a2
export tate_section_a2


@doc Markdown.doc"""
    tate_section_a3(t::TateModelOverGeneralBaseSpace)

Return the Tate section ``a_3``.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a10,a21,a32,a43,a65,w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = GlobalTateModel(ais, auxiliary_base_ring)
A global Tate model over a general base space

julia> tate_section_a3(t)
a32*w^2
```
"""
@attr MPolyElem{fmpq} tate_section_a3(t::TateModelOverGeneralBaseSpace) = t.a3
export tate_section_a3


@doc Markdown.doc"""
    a4(t::TateModelOverGeneralBaseSpace)

Return the Tate section ``a_4``.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a10,a21,a32,a43,a65,w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = GlobalTateModel(ais, auxiliary_base_ring)
A global Tate model over a general base space

julia> tate_section_a4(t)
a43*w^3
```
"""
@attr MPolyElem{fmpq} tate_section_a4(t::TateModelOverGeneralBaseSpace) = t.a4
export tate_section_a4


@doc Markdown.doc"""
    tate_section_a6(t::TateModelOverGeneralBaseSpace)

Return the Tate section ``a_6``.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a10,a21,a32,a43,a65,w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = GlobalTateModel(ais, auxiliary_base_ring)
A global Tate model over a general base space

julia> tate_section_a6(t)
a65*w^5
```
"""
@attr MPolyElem{fmpq} tate_section_a6(t::TateModelOverGeneralBaseSpace) = t.a6
export tate_section_a6


#######################################
# 2: The Tate polynomial
#######################################

@doc Markdown.doc"""
    tate_polynomial(t::TateModelOverGeneralBaseSpace)

Return the Tate polynomial of the global Tate model.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a10,a21,a32,a43,a65,w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = GlobalTateModel(ais, auxiliary_base_ring)
A global Tate model over a general base space

julia> tate_polynomial(t)
a10*x*y*z + a21*w*x^2*z^2 + a32*w^2*y*z^3 + a43*w^3*x*z^4 + a65*w^5*z^6 + x^3 - y^2
```
"""
@attr MPolyElem{fmpq} tate_polynomial(t::TateModelOverGeneralBaseSpace) = t.pt
export tate_polynomial


#######################################
# 3: Auxiliary toric spaces
#######################################

@doc Markdown.doc"""
    auxiliary_base_space(t::TateModelOverGeneralBaseSpace)

Return the toric base space of the global Tate model.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a10,a21,a32,a43,a65,w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = GlobalTateModel(ais, auxiliary_base_ring)
A global Tate model over a general base space

julia> auxiliary_base_space(t)
A normal, affine, 6-dimensional toric variety
```
"""
@attr Oscar.AbstractNormalToricVariety auxiliary_base_space(t::TateModelOverGeneralBaseSpace) = t.auxiliary_base_space
export auxiliary_base_space


@doc Markdown.doc"""
    auxiliary_ambient_space(t::TateModelOverGeneralBaseSpace)

Return the toric ambient space of the global Tate model.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a10,a21,a32,a43,a65,w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = GlobalTateModel(ais, auxiliary_base_ring)
A global Tate model over a general base space

julia> auxiliary_ambient_space(t)
A normal toric variety

julia> dim(auxiliary_ambient_space(t))
8
```
"""
@attr Oscar.AbstractNormalToricVariety auxiliary_ambient_space(t::TateModelOverGeneralBaseSpace) = t.auxiliary_ambient_space
export auxiliary_ambient_space
