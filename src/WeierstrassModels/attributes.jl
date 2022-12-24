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

# Currently no plan to include


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

Return the singular loci of the global Weierstrass model, along with the order of
vanishing of ``(f, g, \Delta)`` at each locus.

For the time being, we either explicitly or implicitly focus on toric varieties
as base spaces. Explicitly, in case the user provides such a variety as base space,
and implicitly, in case we work over a non-fully specified base. This has the
advantage that we can "filter out" trivial singular loci.

Specifically, recall that every closed subvariety of a simplicial toric variety is
of the form ``V(I)``, where ``I`` is a homogeneous ideal of the Cox ring. Let ``B``
be the irrelevant ideal of this toric variety. Then, by proposition 5.2.6. of
[CLS11](@cite), ``V(I)`` is trivial/empty iff ``B^l \subseteq I`` for a suitable ``l \geq 0``.
This can be checked by checking if the saturation ``I:B^\infty`` is the ideal generated by ``1``.

By treating a not-fully specified base space implicitly as toric space, we can extend this
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
sections ``f`` and ``g`` vanish to orders 1 and 2, respectively.

Let us verify this by turning this global Tate model into a Weierstrass model and then
computing the singular loci of the discriminant of this Weierstrass model.

```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a11, a21, a31, a41, a62, w) = QQ["a11", "a21", "a31", "a41", "a62", "w"];

julia> a1 = a11 * w;

julia> a2 = a21 * w;

julia> a3 = a31 * w;

julia> a4 = a41 * w;

julia> a6 = a62 * w^2;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = GlobalTateModel(ais, auxiliary_base_ring, 3)
A global Tate model over a not fully specified base

julia> weier = global_weierstrass_model(t)
[ Info: Base space was not fully specified. Returning AUXILIARY ambient space.
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
A global Weierstrass model over a not fully specified base

julia> length(singular_loci(weier))
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
2

julia> singular_loci(weier)[2]
(ideal(w), (1, 2, 3))
```
"""
@attr Vector{Tuple{MPolyIdeal{MPolyElem_dec{fmpq, fmpq_mpoly}}, Tuple{Int64, Int64, Int64}}} function singular_loci(w::GlobalWeierstrassModel)
    B = irrelevant_ideal(toric_base_space(w))

    d_primes = primary_decomposition(ideal([discriminant(w)]))
    nontrivial_d_primes = Tuple{MPolyIdeal{MPolyElem_dec{fmpq, fmpq_mpoly}}, MPolyIdeal{MPolyElem_dec{fmpq, fmpq_mpoly}}}[]
    for k in 1:length(d_primes)
        if !(is_one(d_primes[k][2]) || is_one(saturation(d_primes[k][2], B)))
            push!(nontrivial_d_primes, d_primes[k])
        end
    end

    f_primes = primary_decomposition(ideal([w.poly_f]))
    g_primes = primary_decomposition(ideal([w.poly_g]))

    kodaira_types = Tuple{MPolyIdeal{MPolyElem_dec{fmpq, fmpq_mpoly}}, Tuple{Int64, Int64, Int64}}[]
    for d_prime in nontrivial_d_primes
        f_index = findfirst(fp -> fp[2] == d_prime[2], f_primes)
        g_index = findfirst(gp -> gp[2] == d_prime[2], g_primes)

        f_order = !isnothing(f_index) ? saturation_with_index(f_primes[f_index][1], d_prime[2])[2] : 0
        g_order = !isnothing(g_index) ? saturation_with_index(g_primes[g_index][1], d_prime[2])[2] : 0
        d_order = saturation_with_index(d_prime[1], d_prime[2])[2]

        push!(kodaira_types, (d_prime[2], (f_order, g_order, d_order)))
    end
    
    return kodaira_types
end
export singular_loci
