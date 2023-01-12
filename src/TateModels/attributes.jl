#####################################################
# 1: The Tate sections
#####################################################

@doc Markdown.doc"""
    tate_section_a1(t::GlobalTateModel)

Return the Tate section ``a_1``.

```jldoctest
julia> using Oscar

julia> t = GlobalTateModel(TestBase())
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

julia> t = GlobalTateModel(TestBase())
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

julia> t = GlobalTateModel(TestBase())
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

julia> t = GlobalTateModel(TestBase())
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

julia> t = GlobalTateModel(TestBase())
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

julia> t = GlobalTateModel(TestBase())
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

julia> t = GlobalTateModel(TestBase())
A global Tate model over a concrete base

julia> toric_base_space(t)
A normal toric variety without torusfactor
```
"""
@attr Oscar.AbstractNormalToricVariety function toric_base_space(t::GlobalTateModel)
    base_fully_specified(t) || @info("Base space was not fully specified. Returning AUXILIARY base space.")
    return t.toric_base_space
end
export toric_base_space


@doc Markdown.doc"""
    toric_ambient_space(t::GlobalTateModel)

Return the toric ambient space of the global Tate model.

```jldoctest
julia> using Oscar

julia> t = GlobalTateModel(TestBase())
A global Tate model over a concrete base

julia> toric_ambient_space(t)
A normal, simplicial toric variety

julia> is_smooth(toric_ambient_space(t))
false
```
"""
@attr Oscar.AbstractNormalToricVariety function toric_ambient_space(t::GlobalTateModel)
    base_fully_specified(t) || @info("Base space was not fully specified. Returning AUXILIARY ambient space.")
    return t.toric_ambient_space
end
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

julia> t = GlobalTateModel(TestBase())
A global Tate model over a concrete base

julia> cy_hypersurface(t)
A closed subvariety of a normal toric variety
```
"""
@attr Oscar.ClosedSubvarietyOfToricVariety function cy_hypersurface(t::GlobalTateModel)
    base_fully_specified(t) || @info("Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.")
    return t.Y4
end
export cy_hypersurface


#####################################################
# 5: Turn global Tate model into a Weierstrass model
#####################################################

@doc Markdown.doc"""
    global_weierstrass_model(t::GlobalTateModel)

Return the global Weierstrass model which is equivalent to the given Tate model.

```jldoctest
julia> using Oscar

julia> t = GlobalTateModel(TestBase())
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
    ring_map = hom(parent(f), S, [gens(S)[i] for i in 1:length(gens(parent(f)))])
    pw = x^3 - y^2 + ring_map(f)*x*z^4 + ring_map(g)*z^6
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

julia> t = GlobalTateModel(TestBase())
A global Tate model over a concrete base

julia> discriminant(t);
```
"""
@attr MPolyElem{fmpq} Oscar.:discriminant(t::GlobalTateModel) = discriminant(global_weierstrass_model(t))
export discriminant


@doc Markdown.doc"""
    singular_loci(t::GlobalTateModel)

Return the singular loci of the global Tate model, along with the order of
vanishing of ``(f, g, \Delta)``` at each locus and the refined Tate fiber type.

For the time being, we either explicitly or implicitly focus on toric varieties
as base spaces. Explicitly, in case the user provides such a variety as base space,
and implicitly, in case we work over a non-fully specified base. This has the
advantage that we can "filter out" trivial singular loci.

Specifically, recall that every closed subvariety of a simplicial toric variety is
of the form ``V(I)``, where ``I`` is a homogeneous ideal of the Cox ring. Let ``B``
be the irrelevant ideal of this toric variety. Then, by proposition 5.2.6. of
[CLS11](@cite), ``V(I)`` is trivial/empty iff ``B^l \subseteq I`` for a suitable ``l \geq 0``.
This can be checked by checking if the saturation ``I:B^\infty`` is the ideal generated by ``1``.

By treating a non-fully specified base space implicitly as a toric space, we can extend this
result straightforwardly to this situation also. This is the reason for constructing this
auxiliary base space.

Let us demonstrate the functionality by computing the singular loci of a Type ``III`` Tate model
[KMSS11](@cite). In this case, we  will consider a global Tate model over a non-fully specified base.
The Tate sections are factored as follows:
- ``a_1 = a_{11} w^1``,
- ``a_2 = a_{21} w^1``,
- ``a_3 = a_{31} w^1``,
- ``a_4 = a_{41} w^1``,
- ``a_6 = a_{62} w^2``.
For this factorization, we expect a singularity of Kodaira type ``III`` over the divisor
``W = {w = 0}``, as desired. So this should be one irreducible component of the discriminant. Moreover,
we should find that the discriminant vanishes to order 3 on ``W = {w = 0}``, while the Weierstrass
sections ``f`` and ``g`` vanish to orders 1 and 2, respectively. Let us verify this.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a11, a21, a31, a41, a62, w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a11 * w;

julia> a2 = a21 * w;

julia> a3 = a31 * w;

julia> a4 = a41 * w;

julia> a6 = a62 * w^2;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = GlobalTateModel(ais, auxiliary_base_ring, 3)
A global Tate model over a not fully specified base

julia> length(singular_loci(t))
[ Info: Base space was not fully specified. Returning AUXILIARY ambient space.
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
2

julia> singular_loci(t)[2]
(ideal(w), (1, 2, 3), "III")
```
"""
@attr Vector{Tuple{MPolyIdeal{MPolyElem_dec{fmpq, fmpq_mpoly}}, Tuple{Int64, Int64, Int64}, String}} singular_loci(t::GlobalTateModel) = singular_loci(global_weierstrass_model(t))
export singular_loci
