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
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> tate_section_a1(t)
a10
```
"""
tate_section_a1(t::GlobalTateModel) = t.tate_a1


@doc raw"""
    tate_section_a2(t::GlobalTateModel)

Return the Tate section ``a_2``.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> tate_section_a2(t)
a21*w
```
"""
tate_section_a2(t::GlobalTateModel) = t.tate_a2


@doc raw"""
    tate_section_a3(t::GlobalTateModel)

Return the Tate section ``a_3``.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> tate_section_a3(t)
a32*w^2
```
"""
tate_section_a3(t::GlobalTateModel) = t.tate_a3


@doc raw"""
    tate_section_a4(t::GlobalTateModel)

Return the Tate section ``a_4``.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> tate_section_a4(t)
a43*w^3
```
"""
tate_section_a4(t::GlobalTateModel) = t.tate_a4


@doc raw"""
    tate_section_a6(t::GlobalTateModel)

Return the Tate section ``a_6``.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> tate_section_a6(t)
0
```
"""
tate_section_a6(t::GlobalTateModel) = t.tate_a6


#####################################################
# 1.2 Tate polynomial
#####################################################

@doc raw"""
    tate_polynomial(t::GlobalTateModel)

Return the Tate polynomial of the global Tate model.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> tate_polynomial(t)
-a10*x*y*z + a21*w*x^2*z^2 - a32*w^2*y*z^3 + a43*w^3*x*z^4 + x^3 - y^2
```
"""
tate_polynomial(t::GlobalTateModel) = t.tate_polynomial


#####################################################
# 1.3 Base, ambient space and fiber ambient space
#####################################################

@doc raw"""
    base_space(t::GlobalTateModel)

Return the base space of the global Tate model.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> base_space(t)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0]]
```
"""
function base_space(t::GlobalTateModel)
  base_fully_specified(t) || @vprint :GlobalTateModel 1 "Base space was not fully specified. Returning AUXILIARY base space.\n"
  return t.base_space
end


@doc raw"""
    ambient_space(t::GlobalTateModel)

Return the ambient space of the global Tate model.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> ambient_space(t)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0, 0, 0, 0, 0, -2, -3], [0, 0, 0, 0, 1, 0, -2, -3], [0, 0, 0, 0, 0, 1, -2, -3], [0, 1, 0, 0, 0, 0, -2, -3], [0, 0, 1, 0, 0, 0, -2, -3], [0, 0, 0, 1, 0, 0, -2, -3], [0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, -1, -3//2]]
```
"""
function ambient_space(t::GlobalTateModel)
  base_fully_specified(t) || @vprint :GlobalTateModel 1 "Base space was not fully specified. Returning AUXILIARY ambient space.\n"
  return t.ambient_space
end


@doc raw"""
    fiber_ambient_space(t::GlobalTateModel)

Return the fiber ambient space of the global Tate model.

```jldoctest
julia> t = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> fiber_ambient_space(t)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[-1, 1//3], [1, -1//2], [0, 1]]
```
"""
fiber_ambient_space(t::GlobalTateModel) = t.fiber_ambient_space





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
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> calabi_yau_hypersurface(t)
Closed subvariety of a normal toric variety
```
"""
@attr ClosedSubvarietyOfToricVariety function calabi_yau_hypersurface(t::GlobalTateModel)
  @req typeof(base_space(t)) <: ToricCoveredScheme "Calabi-Yau hypersurface currently only supported for toric varieties/schemes as base space"
  base_fully_specified(t) || @vprint :GlobalTateModel 1 "Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.\n"
  return closed_subvariety_of_toric_variety(underlying_toric_variety(ambient_space(t)), [tate_polynomial(t)])
end


#####################################################
# 2.2 Turn a Tate into a Weierstrass model
#####################################################

@doc raw"""
    global_weierstrass_model(t::GlobalTateModel)

Return the global Weierstrass model which is equivalent to the given Tate model.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> global_weierstrass_model(t)
Global Weierstrass model over a not fully specified base
```
"""
@attr GlobalWeierstrassModel function global_weierstrass_model(t::GlobalTateModel)
  @req typeof(base_space(t)) <: ToricCoveredScheme "Conversion of global Tate model into global Weierstrass model is currently only supported for toric varieties/schemes as base space"
  b2 = 4 * tate_section_a2(t) + tate_section_a1(t)^2
  b4 = 2 * tate_section_a4(t) + tate_section_a1(t) * tate_section_a3(t)
  b6 = 4 * tate_section_a6(t) + tate_section_a3(t)^2
  f = - 1//48 * (b2^2 - 24 * b4)
  g = 1//864 * (b2^3 - 36 * b2 * b4 + 216 * b6)
  S = cox_ring(ambient_space(t))
  x, y, z = gens(S)[ngens(S)-2:ngens(S)]
  ring_map = hom(parent(f), S, gens(S)[1:ngens(parent(f))])
  pw = x^3 - y^2 + ring_map(f)*x*z^4 + ring_map(g)*z^6
  model = GlobalWeierstrassModel(f, g, pw, base_space(t), ambient_space(t))
  set_attribute!(model, :base_fully_specified, base_fully_specified(t))
  return model
end


#####################################################
# 2.3 Discriminant and singular loci
#####################################################

@doc raw"""
    discriminant(t::GlobalTateModel)

Return the discriminant of the global Tate model.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> discriminant(t);
```
"""
@attr MPolyRingElem function discriminant(t::GlobalTateModel)
  @req typeof(base_space(t)) <: ToricCoveredScheme "Discriminant of global Tate model is currently only supported for toric varieties/schemes as base space"
  return discriminant(global_weierstrass_model(t))
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

julia> a1 = a11 * w;

julia> a2 = a21 * w;

julia> a3 = a31 * w;

julia> a4 = a41 * w;

julia> a6 = a62 * w^2;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = global_tate_model(ais, auxiliary_base_ring, 3)
Global Tate model over a not fully specified base

julia> length(singular_loci(t))
2

julia> singular_loci(t)[2]
(ideal(w), (1, 2, 3), "III")
```
"""
@attr Vector{<:Tuple{<:MPolyIdeal{<:MPolyRingElem}, Tuple{Int64, Int64, Int64}, String}} function singular_loci(t::GlobalTateModel)
  @req typeof(base_space(t)) <: ToricCoveredScheme "Singular loci of global Tate model currently only supported for toric varieties/schemes as base space"
  return singular_loci(global_weierstrass_model(t))
end


#####################################################
# 7: Resolutions for a global Tate model
#####################################################

@doc raw"""
    resolutions(t::GlobalTateModel)

Return the list of all known resolutions for the global
Tate model.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> length(resolutions(t))
2
```
"""
function resolutions(t::GlobalTateModel)
  @req has_attribute(t, :resolutions) "No resolutions known for this model"
  return get_attribute(t, :resolutions)
end


#####################################################
# 8: Attributes for literature Tate models
#####################################################

@doc raw"""
    doi(t::GlobalTateModel)

Return the `doi` of the publication which introduced
the global Tate model in question. If no `doi` is
known, an error is raised.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> doi(t)
"10.1016/j.nuclphysb.2011.12.013"
```
"""
function doi(t::GlobalTateModel)
  @req has_attribute(t, :doi) "No doi known for this model"
  return get_attribute(t, :doi)
end

@doc raw"""
    arxiv_id(t::GlobalTateModel)

Return the `arxiv_id` of the publication which introduced
the global Tate model in question. If no `arxiv_id` is
known, an error is raised.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> arxiv_id(t)
"1109.3454"
```
"""
function arxiv_id(t::GlobalTateModel)
  @req has_attribute(t, :arxiv_id) "No arxiv identifier known for this model"
  return get_attribute(t, :arxiv_id)
end

@doc raw"""
    version(t::GlobalTateModel)

Return the `version` of the arxiv publication which
introduced the global Tate model in question. If no
`version` is known, an error is raised.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> version(t)
"2"
```
"""
function version(t::GlobalTateModel)
  @req has_attribute(t, :version) "No version known for this model"
  return get_attribute(t, :version)
end

@doc raw"""
    equ_nr(t::GlobalTateModel)

Return the `equ_nr` in which the model was introduced
in the publication on our record (either specified by
`doi` or, if uploaded to the arxiv, the `arxiv_id` and
`version`). If no `equ_nr` is known, an error is raised.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> equ_nr(t)
"3.5"
```
"""
function equ_nr(t::GlobalTateModel)
  @req has_attribute(t, :equ_nr) "No equation number known for this model"
  return get_attribute(t, :equ_nr)
end

@doc raw"""
    description(t::GlobalTateModel)

Return the `description` of the model in question.
If no `description` is known, an error is raised.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> description(t)
"SU(5)xU(1) restricted Tate model"
```
"""
function description(t::GlobalTateModel)
  @req has_attribute(t, :description) "No description known for this model"
  return get_attribute(t, :description)
end

@doc raw"""
    link(t::GlobalTateModel)

Return the `link` (formatted as string) to the online
version of the paper which introduced the model in question.
If no link is known, an error is raised.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> link(t)
"https://export.arxiv.org/pdf/1109.3454"
```
"""
function link(t::GlobalTateModel)
  @req has_attribute(t, :link) "No link known for this model"
  return get_attribute(t, :link)
end
