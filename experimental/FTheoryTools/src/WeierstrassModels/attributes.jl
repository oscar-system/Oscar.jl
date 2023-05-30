###################################################################
###################################################################
# 1: Attributes that work the same tor toric and non-toric settings
###################################################################
###################################################################


#####################################################
# 1.1 Tate sections and Tate polynomial
#####################################################

@doc raw"""
    weierstrass_section_f(w::GlobalWeierstrassModel)

Return the polynomial ``f`` used for the
construction of the global Weierstrass model.

```jldoctest
julia> w = su5_weierstrass_model_over_arbitrary_3d_base()
Global Weierstrass model over a not fully specified base

julia> weierstrass_section_f(w);
```
"""
weierstrass_section_f(w::GlobalWeierstrassModel) = w.weierstrass_f


@doc raw"""
    weierstrass_section_g(w::GlobalWeierstrassModel)

Return the polynomial ``g`` used for the
construction of the global Weierstrass model.

```jldoctest
julia> w = su5_weierstrass_model_over_arbitrary_3d_base()
Global Weierstrass model over a not fully specified base

julia> weierstrass_section_g(w);
```
"""
weierstrass_section_g(w::GlobalWeierstrassModel) = w.weierstrass_g


#####################################################
# 1.2 Weierstrass polynomial
#####################################################

@doc raw"""
    weierstrass_polynomial(w::GlobalWeierstrassModel)

Return the Weierstrass polynomial of the global Weierstrass model.

```jldoctest
julia> w = su5_weierstrass_model_over_arbitrary_3d_base()
Global Weierstrass model over a not fully specified base

julia> weierstrass_polynomial(w);
```
"""
weierstrass_polynomial(w::GlobalWeierstrassModel) = w.weierstrass_polynomial


#####################################################
# 1.3 Base, ambient space and fiber ambient space
#####################################################

@doc raw"""
    base_space(w::GlobalWeierstrassModel)

Return the base space of the global Weierstrass model.

```jldoctest
julia> w = su5_weierstrass_model_over_arbitrary_3d_base()
Global Weierstrass model over a not fully specified base

julia> base_space(w)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0]]
"""
function base_space(w::GlobalWeierstrassModel)
  base_fully_specified(w) || @vprint :GlobalWeierstrassModel 1 "Base space was not fully specified. Returning AUXILIARY base space.\n"
  return w.base_space
end


@doc raw"""
    ambient_space(w::GlobalWeierstrassModel)

Return the base space of the global Weierstrass model.

```jldoctest
julia> w = su5_weierstrass_model_over_arbitrary_3d_base()
Global Weierstrass model over a not fully specified base

julia> ambient_space(w)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0, 0, 0, 0, 0, -2, -3], [0, 0, 0, 0, 1, 0, -2, -3], [0, 0, 0, 0, 0, 1, -2, -3], [0, 1, 0, 0, 0, 0, -2, -3], [0, 0, 1, 0, 0, 0, -2, -3], [0, 0, 0, 1, 0, 0, -2, -3], [0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, -1, -3//2]]
```
"""
function ambient_space(w::GlobalWeierstrassModel)
  base_fully_specified(w) || @vprint :GlobalWeierstrassModel 1 "Base space was not fully specified. Returning AUXILIARY ambient space.\n"
  return w.ambient_space
end


@doc raw"""
    fiber_ambient_space(w::GlobalWeierstrassModel)

Return the fiber ambient space of the global Weierstrass model.

```jldoctest
julia> w = su5_weierstrass_model_over_arbitrary_3d_base()
Global Weierstrass model over a not fully specified base

julia> fiber_ambient_space(w)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[-1, 1//3], [1, -1//2], [0, 1]]
```
"""
fiber_ambient_space(w::GlobalWeierstrassModel) = w.fiber_ambient_space





###################################################################
###################################################################
# 2: Attributes that currently only works in toric settings
###################################################################
###################################################################


#####################################################
# 2.1 Calabi-Yau hypersurface
#####################################################

@doc raw"""
    calabi_yau_hypersurface(w::GlobalWeierstrassModel)

Return the Calabi-Yau hypersurface in the toric ambient space
which defines the global Weierstrass model.

```jldoctest
julia> w = su5_weierstrass_model_over_arbitrary_3d_base()
Global Weierstrass model over a not fully specified base

julia> calabi_yau_hypersurface(w)
Closed subvariety of a normal toric variety
```
"""
@attr ClosedSubvarietyOfToricVariety function calabi_yau_hypersurface(w::GlobalWeierstrassModel)
  @req typeof(base_space(w)) <: ToricCoveredScheme "Calabi-Yau hypersurface currently only supported for toric varieties/schemes as base space"
  base_fully_specified(w) || @vprint :GlobalWeierstrassModel 1 "Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.\n"
  return closed_subvariety_of_toric_variety(underlying_toric_variety(ambient_space(w)), [weierstrass_polynomial(w)])
end


#####################################################
# 2.2 Turn global Weierstrass model into Tate model
#####################################################

# Currently no plan to include


#####################################################
# 2.3 Discriminant and singular loci
#####################################################

@doc raw"""
    discriminant(w::GlobalWeierstrassModel)

Return the discriminant ``\Delta = 4 f^3 + 27 g^2``.

```jldoctest
julia> w = su5_weierstrass_model_over_arbitrary_3d_base()
Global Weierstrass model over a not fully specified base

julia> discriminant(w);
```
"""
@attr MPolyRingElem function discriminant(w::GlobalWeierstrassModel)
  @req typeof(base_space(w)) <: ToricCoveredScheme "Discriminant of global Weierstrass model is currently only supported for toric varieties/schemes as base space"
  return 4 * w.weierstrass_f^3 + 27 * w.weierstrass_g^2
end


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
julia> w = su5_weierstrass_model_over_arbitrary_3d_base()
Global Weierstrass model over a not fully specified base

julia> length(singular_loci(w))
2
```
"""
@attr Vector{<:Tuple{<:MPolyIdeal{<:MPolyRingElem}, Tuple{Int64, Int64, Int64}, String}} function singular_loci(w::GlobalWeierstrassModel)
  @req typeof(base_space(w)) <: ToricCoveredScheme "Singular loci of global Weierstrass model is currently only supported for toric varieties/schemes as base space"
  
  B = irrelevant_ideal(base_space(w))
  
  d_primes = primary_decomposition(ideal([discriminant(w)]))
  nontrivial_d_primes = Tuple{<:MPolyIdeal{<:MPolyRingElem}, <:MPolyIdeal{<:MPolyRingElem}}[]
  for k in 1:length(d_primes)
    if _is_nontrivial(d_primes[k][2], B)
      push!(nontrivial_d_primes, d_primes[k])
    end
  end
  
  f_primes = primary_decomposition(ideal([weierstrass_section_f(w)]))
  g_primes = primary_decomposition(ideal([weierstrass_section_g(w)]))
  
  kodaira_types = Tuple{<:MPolyIdeal{<:MPolyRingElem}, Tuple{Int64, Int64, Int64}, String}[]
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
