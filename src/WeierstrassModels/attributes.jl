#######################################
# 1: Weierstrass sections
#######################################

@doc Markdown.doc"""
    weierstrass_section_f(w::GlobalWeierstrassModel)

Return the polynomial ``f`` used for the
construction of the global Weierstrass model.

```jldoctest
julia> using Oscar

julia> w = GlobalWeierstrassModel(TestBase())
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

julia> w = GlobalWeierstrassModel(TestBase())
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

julia> w = GlobalWeierstrassModel(TestBase())
A global Weierstrass model over a concrete base

julia> weierstrass_polynomial(w);
```
"""
@attr MPolyElem{fmpq} weierstrass_polynomial(w::GlobalWeierstrassModel) = w.pw
export weierstrass_polynomial


#######################################
# 3: Toric spaces
#######################################

@doc Markdown.doc"""
    toric_base_space(w::GlobalWeierstrassModel)

Return the toric base space of the global Weierstrass model.

```jldoctest
julia> using Oscar

julia> w = GlobalWeierstrassModel(TestBase())
A global Weierstrass model over a concrete base

julia> is_smooth(toric_base_space(w))
true
```
"""
@attr Oscar.AbstractNormalToricVariety function toric_base_space(w::GlobalWeierstrassModel)
    base_fully_specified(w) || @info("Base space was not fully specified. Returning AUXILIARY base space.")
    return w.toric_base_space
end
export toric_base_space


@doc Markdown.doc"""
    toric_ambient_space(w::GlobalWeierstrassModel)

Return the toric base space of the global Weierstrass model.

```jldoctest
julia> using Oscar

julia> w = GlobalWeierstrassModel(TestBase())
A global Weierstrass model over a concrete base

julia> is_smooth(toric_ambient_space(w))
false
```
"""
@attr Oscar.AbstractNormalToricVariety function toric_ambient_space(w::GlobalWeierstrassModel)
    base_fully_specified(w) || @info("Base space was not fully specified. Returning AUXILIARY ambient space.")
    return w.toric_ambient_space
end
export toric_base_space


#####################################################
# 4: The CY hypersurface
#####################################################

@doc Markdown.doc"""
    cy_hypersurface(w::GlobalWeierstrassModel)

Return the Calabi-Yau hypersurface in the toric ambient space
which defines the global Weierstrass model.

```jldoctest
julia> using Oscar

julia> w = GlobalWeierstrassModel(TestBase())
A global Weierstrass model over a concrete base

julia> cy_hypersurface(w)
A closed subvariety of a normal toric variety
```
"""
@attr Oscar.ClosedSubvarietyOfToricVariety function cy_hypersurface(w::GlobalWeierstrassModel)
    base_fully_specified(w) || @info("Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.")
    return w.Y4
end
export cy_hypersurface


#####################################################
# 5: Turn global Weierstrass model into Tate model
#####################################################

# TODO: To come
# TODO: To come


#######################################
# 6: Discriminant
#######################################

@doc Markdown.doc"""
    discriminant(w::GlobalWeierstrassModel)

Return the discriminant ``\Delta = 4 f^3 + 27 g^2``.

```jldoctest
julia> using Oscar

julia> w = GlobalWeierstrassModel(TestBase())
A global Weierstrass model over a concrete base

julia> discriminant(w);
```
"""
@attr MPolyElem{fmpq} Oscar.:discriminant(w::GlobalWeierstrassModel) = 4 * w.poly_f^3 + 27 * w.poly_g^2
export discriminant


@doc Markdown.doc"""
    singular_loci(w::GlobalWeierstrassModel)

Return the primary decomposition of the discriminant of the global Tate model.

For the time being, we either explicitly or implicitly focuse on toric varieties
as base spaces. Explicitly, in case the user provides such a variety as base space.
And implicitly, in case we work over a not-fully specified base. This has the
advantage that we can "filter out" trivial singular loci.

Specifically, recall that every closed subvariety of a simplicial toric variety is
of the form ``V(I)``, where ``I`` is a homogeneous ideal of the Cox ring. Let ``B``
be the irrelevant ideal of this toric variety. Then, by proposition 5.2.6. of
[CLS11](@cite), ``V(I)`` is trivial/empty iff ``B^l \subseteq I`` for a suitable ``l \geq 0``.
This can be checked by checking if the saturation ``I:B^\infty`` is the ideal generated by ``1``.

By treating a not-fully specified base space implicitly as toric space, we can extend this
result straightforwardly to this situation also. This is the reason for constructing this
auxiliary base space.

Let us demonstrate the functionality by computing the example around equ. (4.50) in
[Wei18](@cite). In this example, one considers a global Tate model over a fully-specified based.
The Tate sections as factored as follows:
- ``a_1``,
- ``a_2 = a_{21} w^1``,
- ``a_3 = a_{31} w^1``,
- ``a_4 = a_{41} w^1``,
- ``a_6 = a_{62} w^2``.
For this factorization we expect a singularity of Kodaira type ``I_2`` over the divisor
``W = {w = 0}``. So this should be one irreducible component of the discriminant. Even more,
we should find that the discriminant vanishes to order 2 on ``W = {w = 0}``.

Let us verify this by turning this global Tate model into a Weierstrass model and then
studying the primary decomposition of the discriminant of this Weierstrass model.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a10, a21, a31, a41, a62, w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a31 * w;

julia> a4 = a41 * w;

julia> a6 = a62 * w^2;

julia> t = GlobalTateModel([a1, a2, a3, a4, a6], auxiliary_base_ring, 3)
A global Tate model over a not fully specified base

julia> weier = global_weierstrass_model(t)
[ Info: Base space was not fully specified. Returning AUXILIARY ambient space.
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
A global Weierstrass model over a not fully specified base

julia> length(singular_loci(weier))
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
2

julia> singular_loci(weier)[2]
(ideal(w^2), ideal(w))
```
"""
@attr Vector{Tuple{MPolyIdeal{MPolyElem_dec{fmpq, fmpq_mpoly}}, MPolyIdeal{MPolyElem_dec{fmpq, fmpq_mpoly}}}} function singular_loci(w::GlobalWeierstrassModel)
    B = irrelevant_ideal(toric_base_space(w))
    primes = primary_decomposition(ideal([discriminant(w)]))
    nontrivial_primes = Tuple{MPolyIdeal{MPolyElem_dec{fmpq, fmpq_mpoly}}, MPolyIdeal{MPolyElem_dec{fmpq, fmpq_mpoly}}}[]
    for k in 1:length(primes)
        if !(is_one(primes[k][2]) || is_one(saturation(primes[k][2], B)))
            push!(nontrivial_primes, primes[k])
        end
    end
    return nontrivial_primes
end
export singular_loci
