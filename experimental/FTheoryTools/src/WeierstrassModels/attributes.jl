#######################################
# 1: Weierstrass sections
#######################################

@doc raw"""
    weierstrass_section_f(w::GlobalWeierstrassModel)

Return the polynomial ``f`` used for the
construction of the global Weierstrass model.

```jldoctest
julia> auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"];

julia> w = global_weierstrass_model(f, g, auxiliary_base_ring, 3)
Global Weierstrass model over a not fully specified base

julia> weierstrass_section_f(w);
```
"""
@attr MPolyRingElem{QQFieldElem} weierstrass_section_f(w::GlobalWeierstrassModel) = w.poly_f


@doc raw"""
    weierstrass_section_g(w::GlobalWeierstrassModel)

Return the polynomial ``g`` used for the
construction of the global Weierstrass model.

```jldoctest
julia> auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"];

julia> w = global_weierstrass_model(f, g, auxiliary_base_ring, 3)
Global Weierstrass model over a not fully specified base

julia> weierstrass_section_g(w);
```
"""
@attr MPolyRingElem{QQFieldElem} weierstrass_section_g(w::GlobalWeierstrassModel) = w.poly_g


#######################################
# 2: Weierstrass polynomial
#######################################

@doc raw"""
    weierstrass_polynomial(w::GlobalWeierstrassModel)

Return the Weierstrass polynomial of the global Weierstrass model.

```jldoctest
julia> auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"];

julia> w = global_weierstrass_model(f, g, auxiliary_base_ring, 3)
Global Weierstrass model over a not fully specified base

julia> weierstrass_polynomial(w);
```
"""
@attr MPolyRingElem{QQFieldElem} weierstrass_polynomial(w::GlobalWeierstrassModel) = w.pw


#######################################
# 3: Toric spaces
#######################################

@doc raw"""
    toric_base_space(w::GlobalWeierstrassModel)

Return the toric base space of the global Weierstrass model.

```jldoctest
julia> auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"];

julia> w = global_weierstrass_model(f, g, auxiliary_base_ring, 3)
Global Weierstrass model over a not fully specified base

julia> dim(toric_base_space(w))
3
```
"""
@attr AbstractNormalToricVariety function toric_base_space(w::GlobalWeierstrassModel)
    base_fully_specified(w) || @vprint :GlobalWeierstrassModel 1 "Base space was not fully specified. Returning AUXILIARY base space.\n"
    return w.toric_base_space
end


@doc raw"""
    toric_ambient_space(w::GlobalWeierstrassModel)

Return the toric base space of the global Weierstrass model.

```jldoctest
julia> auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"];

julia> w = global_weierstrass_model(f, g, auxiliary_base_ring, 3)
Global Weierstrass model over a not fully specified base

julia> dim(toric_ambient_space(w))
5
```
"""
@attr AbstractNormalToricVariety function toric_ambient_space(w::GlobalWeierstrassModel)
    base_fully_specified(w) || @vprint :GlobalWeierstrassModel 1 "Base space was not fully specified. Returning AUXILIARY ambient space.\n"
    return w.toric_ambient_space
end


#####################################################
# 4: The CY hypersurface
#####################################################

@doc raw"""
    calabi_yau_hypersurface(w::GlobalWeierstrassModel)

Return the Calabi-Yau hypersurface in the toric ambient space
which defines the global Weierstrass model.

```jldoctest
julia> auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"];

julia> w = global_weierstrass_model(f, g, auxiliary_base_ring, 3)
Global Weierstrass model over a not fully specified base

julia> calabi_yau_hypersurface(w)
Closed subvariety of a normal toric variety
```
"""
@attr ClosedSubvarietyOfToricVariety function calabi_yau_hypersurface(w::GlobalWeierstrassModel)
    base_fully_specified(w) || @vprint :GlobalWeierstrassModel 1 "Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.\n"
    return w.calabi_yau_hypersurface
end


#####################################################
# 5: Turn global Weierstrass model into Tate model
#####################################################

# Currently no plan to include


#######################################
# 6: Discriminant
#######################################

@doc raw"""
    discriminant(w::GlobalWeierstrassModel)

Return the discriminant ``\Delta = 4 f^3 + 27 g^2``.

```jldoctest
julia> auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"];

julia> w = global_weierstrass_model(f, g, auxiliary_base_ring, 3)
Global Weierstrass model over a not fully specified base

julia> discriminant(w);
```
"""
@attr MPolyRingElem{QQFieldElem} discriminant(w::GlobalWeierstrassModel) = 4 * w.poly_f^3 + 27 * w.poly_g^2


@doc raw"""
    singular_loci(w::GlobalWeierstrassModel)

Return the singular loci of the global Weierstrass model, along with the order of
vanishing of ``(f, g, \Delta)`` at each locus and the refined Tate fiber type.

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

```jldoctest
julia> auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"];

julia> w = global_weierstrass_model(f, g, auxiliary_base_ring, 3)
Global Weierstrass model over a not fully specified base

julia> discriminant(w);

julia> length(singular_loci(w))
1
```
"""
@attr Vector{Tuple{MPolyIdeal{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}, Tuple{Int64, Int64, Int64}, String}} function singular_loci(w::GlobalWeierstrassModel)
    B = irrelevant_ideal(toric_base_space(w))

    d_primes = primary_decomposition(ideal([discriminant(w)]))
    nontrivial_d_primes = Tuple{MPolyIdeal{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}, MPolyIdeal{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}[]
    for k in 1:length(d_primes)
        if _is_nontrivial(d_primes[k][2], B)
            push!(nontrivial_d_primes, d_primes[k])
        end
    end

    f_primes = primary_decomposition(ideal([weierstrass_section_f(w)]))
    g_primes = primary_decomposition(ideal([weierstrass_section_g(w)]))

    kodaira_types = Tuple{MPolyIdeal{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}, Tuple{Int64, Int64, Int64}, String}[]
    for d_prime in nontrivial_d_primes
        f_index = findfirst(fp -> fp[2] == d_prime[2], f_primes)
        g_index = findfirst(gp -> gp[2] == d_prime[2], g_primes)

        f_order = !isnothing(f_index) ? saturation_with_index(f_primes[f_index][1], d_prime[2])[2] : 0
        g_order = !isnothing(g_index) ? saturation_with_index(g_primes[g_index][1], d_prime[2])[2] : 0
        d_order = saturation_with_index(d_prime[1], d_prime[2])[2]
        ords = (f_order, g_order, d_order)

        push!(kodaira_types, (d_prime[2], ords, _kodaira_type(d_prime[2], weierstrass_section_f(w), weierstrass_section_g(w), discriminant(w), ords)))
    end
    
    return kodaira_types
end
