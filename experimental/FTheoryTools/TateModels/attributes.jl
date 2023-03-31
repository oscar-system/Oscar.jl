#####################################################
# 1: The Tate sections
#####################################################

@doc Markdown.doc"""
    tate_section_a1(t::GlobalTateModel)

Return the Tate section ``a_1``.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> tate_section_a1(t);
```
"""
@attr MPolyRingElem{QQFieldElem} tate_section_a1(t::GlobalTateModel) = t.a1


@doc Markdown.doc"""
    tate_section_a2(t::GlobalTateModel)

Return the Tate section ``a_2``.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> tate_section_a2(t);
```
"""
@attr MPolyRingElem{QQFieldElem} tate_section_a2(t::GlobalTateModel) = t.a2


@doc Markdown.doc"""
    tate_section_a3(t::GlobalTateModel)

Return the Tate section ``a_3``.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> tate_section_a3(t);
```
"""
@attr MPolyRingElem{QQFieldElem} tate_section_a3(t::GlobalTateModel) = t.a3


@doc Markdown.doc"""
    tate_section_a4(t::GlobalTateModel)

Return the Tate section ``a_4``.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> tate_section_a4(t);
```
"""
@attr MPolyRingElem{QQFieldElem} tate_section_a4(t::GlobalTateModel) = t.a4


@doc Markdown.doc"""
    tate_section_a6(t::GlobalTateModel)

Return the Tate section ``a_6``.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> tate_section_a6(t);
```
"""
@attr MPolyRingElem{QQFieldElem} tate_section_a6(t::GlobalTateModel) = t.a6


#####################################################
# 2: The Tate polynomial
#####################################################

@doc Markdown.doc"""
    tate_polynomial(t::GlobalTateModel)

Return the Tate polynomial of the global Tate model.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> tate_polynomial(t);
```
"""
@attr MPolyRingElem{QQFieldElem} tate_polynomial(t::GlobalTateModel) = t.pt


#####################################################
# 3: Toric spaces for Tate models over concrete bases
#####################################################

@doc Markdown.doc"""
    toric_base_space(t::GlobalTateModel)

Return the toric base space of the global Tate model.

```jldoctest
julia> t = global_tate_model(test_base())
Global Tate model over a concrete base

julia> toric_base_space(t)
Normal toric variety without torusfactor
```
"""
@attr AbstractNormalToricVariety function toric_base_space(t::GlobalTateModel)
    base_fully_specified(t) || @info("Base space was not fully specified. Returning AUXILIARY base space.")
    return t.toric_base_space
end


@doc Markdown.doc"""
    toric_ambient_space(t::GlobalTateModel)

Return the toric ambient space of the global Tate model.

```jldoctest
julia> t = global_tate_model(test_base())
Global Tate model over a concrete base

julia> toric_ambient_space(t)
Normal, simplicial toric variety

julia> is_smooth(toric_ambient_space(t))
false
```
"""
@attr AbstractNormalToricVariety function toric_ambient_space(t::GlobalTateModel)
    base_fully_specified(t) || @info("Base space was not fully specified. Returning AUXILIARY ambient space.")
    return t.toric_ambient_space
end


#####################################################
# 4: The CY hypersurface
#####################################################

@doc Markdown.doc"""
    cy_hypersurface(t::GlobalTateModel)

Return the Calabi-Yau hypersurface in the toric ambient space
which defines the global Tate model.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> cy_hypersurface(t)
[ Info: Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.
Closed subvariety of a normal toric variety
```
"""
@attr ClosedSubvarietyOfToricVariety function cy_hypersurface(t::GlobalTateModel)
    base_fully_specified(t) || @info("Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.")
    return t.Y4
end


#####################################################
# 5: Turn global Tate model into a Weierstrass model
#####################################################

@doc Markdown.doc"""
    global_weierstrass_model(t::GlobalTateModel)

Return the global Weierstrass model which is equivalent to the given Tate model.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> global_weierstrass_model(t)
[ Info: Base space was not fully specified. Returning AUXILIARY ambient space.
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
Global Weierstrass model over a not fully specified base
```
"""
@attr GlobalWeierstrassModel function global_weierstrass_model(t::GlobalTateModel)
    b2 = 4 * tate_section_a2(t) + tate_section_a1(t)^2
    b4 = 2 * tate_section_a4(t) + tate_section_a1(t) * tate_section_a3(t)
    b6 = 4 * tate_section_a6(t) + tate_section_a3(t)^2
    f = - 1//48 * (b2^2 - 24 * b4)
    g = 1//864 * (b2^3 - 36 * b2 * b4 + 216 * b6)
    S = cox_ring(toric_ambient_space(t))
    x, y, z = gens(S)[ngens(S)-2:ngens(S)]
    ring_map = hom(parent(f), S, gens(S)[1:ngens(parent(f))])
    pw = x^3 - y^2 + ring_map(f)*x*z^4 + ring_map(g)*z^6
    Y4 = closed_subvariety_of_toric_variety(toric_ambient_space(t), [pw])
    model = GlobalWeierstrassModel(f, g, pw, toric_base_space(t), toric_ambient_space(t), Y4)
    set_attribute!(model, :base_fully_specified, base_fully_specified(t))
    return model
end


#####################################################
# 6: Discriminant
#####################################################

@doc Markdown.doc"""
    discriminant(t::GlobalTateModel)

Return the discriminant of the global Tate model.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> discriminant(t);
[ Info: Base space was not fully specified. Returning AUXILIARY ambient space.
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
```
"""
@attr MPolyRingElem{QQFieldElem} discriminant(t::GlobalTateModel) = discriminant(global_weierstrass_model(t))


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
[KMSS11](@cite). In this case, we  will consider Global Tate model over a non-fully specified base.
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
julia> auxiliary_base_ring, (a11, a21, a31, a41, a62, w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a11 * w;

julia> a2 = a21 * w;

julia> a3 = a31 * w;

julia> a4 = a41 * w;

julia> a6 = a62 * w^2;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = global_tate_model(ais, auxiliary_base_ring, 3)
Global Tate model over a not fully specified base

julia> length(singular_loci(t))
[ Info: Base space was not fully specified. Returning AUXILIARY ambient space.
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
2

julia> singular_loci(t)[2]
(ideal(w), (1, 2, 3), "III")
```
"""
@attr Vector{Tuple{MPolyIdeal{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}, Tuple{Int64, Int64, Int64}, String}} singular_loci(t::GlobalTateModel) = singular_loci(global_weierstrass_model(t))


#####################################################
# 7: Fiber analysis
#####################################################

function analyze_fibers(model::GlobalTateModel, centers::Vector{<:Vector{<:Integer}})
    # Ideal of the defining polynomial
    hypersurface_ideal = ideal([tate_polynomial(model)])

    # Toric ambient space
    tas = toric_ambient_space(model)

    # Various important ideals
    irr = irrelevant_ideal(tas);
    sri = stanley_reisner_ideal(tas);
    lin = ideal_of_linear_relations(tas);

    # Singular loci
    sing_loc = singular_loci(model)

    # Pick out the singular loci that are more singular than an I_1
    # Then keep only the locus and not the extra info about it
    interesting_singular_loci = map(tup -> tup[1], filter(locus -> locus[2][3] > 1, sing_loc))

    # This is a kludge to map polynomials on the base into the ambient space, and should be fixed once the ambient space constructors supply such a map
    base_coords = parent(gens(interesting_singular_loci[1])[1])
    ambient_coords = parent(tate_polynomial(model))
    base_to_ambient_ring_map = hom(base_coords, ambient_coords, gens(ambient_coords)[1:end-3])

    # Resolved model
    strict_transform, exceptionals, crepant, res_irr, res_sri, res_lin, res_S, res_S_gens, res_ring_map = _blowup_global_sequence(hypersurface_ideal, centers, irr, sri, lin)
    if !crepant
        @warn "The given sequence of blowups is not crepant"
    end

    loci_fiber_intersections = Tuple{MPolyIdeal{QQMPolyRingElem}, Vector{Tuple{Tuple{Int64, Int64}, Vector{MPolyIdeal{QQMPolyRingElem}}}}}[]
    for locus in interesting_singular_loci
        # Currently have to get the ungraded ideal generators by hand using .f
        ungraded_locus = ideal(map(gen -> base_to_ambient_ring_map(gen).f, gens(locus)))

        # Potential components of the fiber over this locus
        # For now, we only consider the associated prime ideal,
        #   but we may later want to actually consider the primary ideals
        potential_components = map(pair -> pair[2], primary_decomposition(strict_transform + res_ring_map(ungraded_locus)))

        # Filter out the trivial loci among the potential components
        components = filter(component -> _is_nontrivial(component, res_irr), potential_components)

        # Check the pairwise intersections of the components
        intersections = Tuple{Tuple{Int64, Int64}, Vector{MPolyIdeal{QQMPolyRingElem}}}[]
        for i in 1:length(components) - 1
            for j in i + 1:length(components)
                intersection = filter(candidate_locus -> _is_nontrivial(candidate_locus, res_irr), map(pair -> pair[2], primary_decomposition(components[i] + components[j])))
                push!(intersections, ((i, j), intersection))
            end
        end

        push!(loci_fiber_intersections, (ungraded_locus, intersections))
    end

    return loci_fiber_intersections
end
