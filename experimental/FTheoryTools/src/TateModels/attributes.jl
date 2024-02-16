###################################################################
###################################################################
# 1: Attributes that work the same tor toric and non-toric settings
###################################################################
###################################################################


#####################################################
# 1.1 Tate sections and Tate polynomial
#####################################################


@doc raw"""
    tate_section_a1(t::GlobalTateModel)

Return the Tate section ``a_1``.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> tate_section_a1(t)
a1
```
"""
tate_section_a1(t::GlobalTateModel) = explicit_model_sections(t)["a1"]


@doc raw"""
    tate_section_a2(t::GlobalTateModel)

Return the Tate section ``a_2``.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> tate_section_a2(t)
w*a21
```
"""
tate_section_a2(t::GlobalTateModel) = explicit_model_sections(t)["a2"]


@doc raw"""
    tate_section_a3(t::GlobalTateModel)

Return the Tate section ``a_3``.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> tate_section_a3(t)
w^2*a32
```
"""
tate_section_a3(t::GlobalTateModel) = explicit_model_sections(t)["a3"]


@doc raw"""
    tate_section_a4(t::GlobalTateModel)

Return the Tate section ``a_4``.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> tate_section_a4(t)
w^3*a43
```
"""
tate_section_a4(t::GlobalTateModel) = explicit_model_sections(t)["a4"]


@doc raw"""
    tate_section_a6(t::GlobalTateModel)

Return the Tate section ``a_6``.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> tate_section_a6(t)
0
```
"""
tate_section_a6(t::GlobalTateModel) = explicit_model_sections(t)["a6"]


#####################################################
# 1.2 Tate polynomial
#####################################################

@doc raw"""
    tate_polynomial(t::GlobalTateModel)

Return the Tate polynomial of the global Tate model.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> tate_polynomial(t)
w^3*a43*x*z^4 - w^2*a32*y*z^3 + w*a21*x^2*z^2 - a1*x*y*z + x^3 - y^2
```
"""
tate_polynomial(t::GlobalTateModel) = t.tate_polynomial



###################################################################
###################################################################
# 2: Attributes that currently only works in toric settings
###################################################################
###################################################################


#####################################################
# 2.1 Calabi-Yau hypersurface
#####################################################

@doc raw"""
    calabi_yau_hypersurface(t::GlobalTateModel)

Return the Calabi-Yau hypersurface in the toric ambient space
which defines the global Tate model.

```jldoctest
julia> t = global_tate_model(sample_toric_variety(); completeness_check = false)
Global Tate model over a concrete base

julia> calabi_yau_hypersurface(t)
Closed subvariety of a normal toric variety
```
"""
@attr ClosedSubvarietyOfToricVariety function calabi_yau_hypersurface(t::GlobalTateModel)
  @req base_space(t) isa NormalToricVariety "Calabi-Yau hypersurface currently only supported for toric varieties as base space"
  is_base_space_fully_specified(t) || @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.\n"
  return closed_subvariety_of_toric_variety(ambient_space(t), [tate_polynomial(t)])
end


#####################################################
# 2.2 Turn a Tate into a Weierstrass model
#####################################################

@doc raw"""
    weierstrass_model(t::GlobalTateModel)

Return the Weierstrass model which is equivalent to the given Tate model.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> weierstrass_model(t)
Weierstrass model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
"""
@attr WeierstrassModel function weierstrass_model(t::GlobalTateModel)
  @req (base_space(t) isa NormalToricVariety || base_space(t) isa FamilyOfSpaces) "Conversion of global Tate model into Weierstrass model is currently only supported for toric varieties and family of spaces as base space"
  
  # Compute explicit Weierstrass sections
  b2 = 4 * tate_section_a2(t) + tate_section_a1(t)^2
  b4 = 2 * tate_section_a4(t) + tate_section_a1(t) * tate_section_a3(t)
  b6 = 4 * tate_section_a6(t) + tate_section_a3(t)^2
  f = - 1//48 * (b2^2 - 24 * b4)
  g = 1//864 * (b2^3 - 36 * b2 * b4 + 216 * b6)

  # Compute explicit_model_sections
  new_explicit_model_sections = Dict("f" => f, "g" => g)
  sections_t = collect(keys(explicit_model_sections(t)))
  for section in sections_t
    if !(section in ["a1", "a2", "a3", "a4", "a6"])
      new_explicit_model_sections[section] = explicit_model_sections(t)[section]
    end
  end

  # Compute Weierstrass polynomial
  S = cox_ring(ambient_space(t))
  x, y, z = gens(S)[ngens(S)-2:ngens(S)]
  ring_map = hom(parent(f), S, gens(S)[1:ngens(parent(f))])
  pw = x^3 - y^2 + ring_map(f)*x*z^4 + ring_map(g)*z^6

  # Compute parametrization of Weierstrass sections
  parametrization = defining_section_parametrization(t)
  param_keys = collect(keys(parametrization))
  new_defining_section_parametrization = Dict{String, MPolyRingElem}()
  if length(param_keys) > 0
    # Find ring to evaluate polynomials into
    R = parent(parametrization[param_keys[1]])

    # Identify how we parametrize a1
    if haskey(parametrization, "a1")
      param_a1 = parametrization["a1"]
    else
      param_a1 = eval_poly("a1", R)
    end

    # Identify how we parametrize a2
    if haskey(parametrization, "a2")
      param_a2 = parametrization["a2"]
    else
      param_a2 = eval_poly("a2", R)
    end

    # Identify how we parametrize a3
    if haskey(parametrization, "a3")
      param_a3 = parametrization["a3"]
    else
      param_a3 = eval_poly("a3", R)
    end

    # Identify how we parametrize a4
    if haskey(parametrization, "a4")
      param_a4 = parametrization["a4"]
    else
      param_a4 = eval_poly("a4", R)
    end

    # Identify how we parametrize a6
    if haskey(parametrization, "a6")
      param_a6 = parametrization["a6"]
    else
      param_a6 = eval_poly("a6", R)
    end

    # Compute parametrization of b2, b4, b6
    param_b2 = 4 * param_a2 + param_a1^2
    param_b4 = 2 * param_a4 + param_a1 * param_a3
    param_b6 = 4 * param_a6 + param_a3^2

    # Compute parametrization of f, g
    param_f = -1//48 * (param_b2^2 - 24 * param_b4)
    param_g = 1//864 * (param_b2^3 - 36 * param_b2 * param_b4 + 216 * param_b6)

    # Compute defining_section_parametrization
    new_defining_section_parametrization = Dict("f" => param_f, "g" => param_g)

  end
  
  # Compute Weierstrass model
  model = WeierstrassModel(new_explicit_model_sections, new_defining_section_parametrization, pw, base_space(t), ambient_space(t))

  # Copy attributes and return model
  model_attributes = t.__attrs
  for (key, value) in model_attributes
    set_attribute!(model, key, value)
  end
  return model
end


#####################################################
# 2.3 Discriminant and singular loci
#####################################################

@doc raw"""
    discriminant(t::GlobalTateModel)

Return the discriminant of the global Tate model.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> discriminant(t);
```
"""
@attr MPolyRingElem function discriminant(t::GlobalTateModel)
  @req (base_space(t) isa NormalToricVariety || base_space(t) isa FamilyOfSpaces) "Discriminant of global Tate model is currently only supported for toric varieties and family of spaces as base space"
  return discriminant(weierstrass_model(t))
end


@doc raw"""
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

julia> auxiliary_base_grading = [1 2 3 4 6 0; -1 -1 -1 -1 -2 1];

julia> a1 = a11 * w;

julia> a2 = a21 * w;

julia> a3 = a31 * w;

julia> a4 = a41 * w;

julia> a6 = a62 * w^2;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = global_tate_model(auxiliary_base_ring, auxiliary_base_grading, 3, ais)
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base

julia> length(singular_loci(t))
2

julia> singular_loci(t)[2]
(Ideal (w), (1, 2, 3), "III")
```
"""
@attr Vector{<:Tuple{<:MPolyIdeal{<:MPolyRingElem}, Tuple{Int64, Int64, Int64}, String}} function singular_loci(t::GlobalTateModel)
  @req (base_space(t) isa NormalToricVariety || base_space(t) isa FamilyOfSpaces) "Singular loci of global Tate model currently only supported for toric varieties and families of spaces as base space"
  return singular_loci(weierstrass_model(t))
end
