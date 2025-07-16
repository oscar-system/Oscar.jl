###################################################################
###################################################################
# 1: Attributes that work the same tor toric and non-toric settings
###################################################################
###################################################################

@doc raw"""
    hypersurface_equation(h::HypersurfaceModel)

Returns the hypersurface equation of the hypersurface model.

```jldoctest
julia> b = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> D1 = 2 * anticanonical_divisor_class(b)
Divisor class on a normal toric variety

julia> D2 = 3 * anticanonical_divisor_class(b)
Divisor class on a normal toric variety

julia> new_gens = string.(vcat(gens(cox_ring(b)), gens(cox_ring(fiber_ambient_space))));

julia> ambient_ring, (x1, x2, x3, x, y, z) = polynomial_ring(QQ, new_gens, cached=false)
(Multivariate polynomial ring in 6 variables over QQ, QQMPolyRingElem[x1, x2, x3, x, y, z])

julia> p = x^3 - y^2 + x1^12 * x * z^4 + x2^18 * z^6 + 13 * x3^3*x*y*z
x1^12*x*z^4 + x2^18*z^6 + 13*x3^3*x*y*z + x^3 - y^2

julia> h = hypersurface_model(b, fiber_ambient_space, [D1, D2], p; completeness_check = false)
Hypersurface model over a concrete base

julia> hypersurface_equation(h)
x1^12*x*z^4 + x2^18*z^6 + 13*x3^3*x*y*z + x^3 - y^2
```
"""
hypersurface_equation(h::HypersurfaceModel) = h.hypersurface_equation


@doc raw"""
    hypersurface_equation_parametrization(h::HypersurfaceModel)

Returns the parametrization of the hypersurface equation by the model sections.

```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> h = literature_model(arxiv_id = "1208.2695", equation = "B.5")
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base

julia> explicit_model_sections(h)
Dict{String, MPolyRingElem} with 5 entries:
  "c2" => c2
  "c1" => c1
  "c3" => c3
  "b"  => b
  "c0" => c0

julia> hypersurface_equation_parametrization(h)
b*w*v^2 - c0*u^4 - c1*u^3*v - c2*u^2*v^2 - c3*u*v^3 + w^2
```
"""
hypersurface_equation_parametrization(h::HypersurfaceModel) = h.hypersurface_equation_parametrization



#########################################################################
#########################################################################
# 2: Attributes that define the corresponding Weierstrass and Tate models
#########################################################################
#########################################################################

@doc raw"""
    weierstrass_model(h::HypersurfaceModel)

Returns the Weierstrass model corresponding to the given hypersurface model, if known.

In the example below, we construct a hypersurface model and its corresponding Weierstrass
model (see [BMT25](@cite BMT25) for background), and then establish the relationship
between the two models.

```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> D1 = 2 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> D2 = 3 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> amb_ring, (x1, x2, x3, x, y, z) = polynomial_ring(QQ, ["x1", "x2", "x3", "x", "y", "z"])
(Multivariate polynomial ring in 6 variables over QQ, QQMPolyRingElem[x1, x2, x3, x, y, z])

julia> p = x^3 + 7*x1*x2^5*x^2*z^2 + x1^3*(x2 + x3)^9*x*z^4 - y^2 - 13*x3^3*x*y*z - x1^2*x2^4*x3^3*y*z^3;

julia> h = hypersurface_model(B2, fiber_ambient_space, [D1, D2], p, completeness_check = false)
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(cox_ring(B2));

julia> weier_f = 1//48*(-(28*x1*x2^5 + 169*x3^6)^2 + 24*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> weier_g = 1//864*(216*x1^4*x2^8*x3^6 + (28*x1*x2^5 + 169*x3^6)^3 - 36*(28*x1*x2^5 + 169*x3^6)*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> w = weierstrass_model(B2, weier_f, weier_g; completeness_check = false)
Weierstrass model over a concrete base

julia> set_weierstrass_model(h, w)

julia> weierstrass_model(h) === w
true
```
"""
function weierstrass_model(h::HypersurfaceModel)
  @req has_attribute(h, :weierstrass_model) "No corresponding Weierstrass model is known"
  w = get_attribute(h, :weierstrass_model)
  if w isa WeierstrassModel
    return w
  end
  @req w isa String "Internal inconsistency encountered"
  directory = joinpath(dirname(@__DIR__), "LiteratureModels/")
  model_indices = JSON.parsefile(directory * "model_indices.json")
  if is_base_space_fully_specified(h)
    w_model = literature_model(parse(Int, model_indices[w]), base_space = base_space(h), defining_classes = defining_classes(h), completeness_check = false)
  else
    w_model = literature_model(parse(Int, model_indices[w]))
  end
  set_weierstrass_model(h, w_model)
  return w_model
end


@doc raw"""
    global_tate_model(h::HypersurfaceModel)

Returns the global Tate model corresponding to the given hypersurface model, if known.

In the example below, we construct a hypersurface model and its corresponding global
Tate model (see [BMT25](@cite BMT25) for background), and then establish the relationship
between the two models.

```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> x1, x2, x3 = gens(cox_ring(B2));

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> D1 = 2 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> D2 = 3 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> amb_ring, (x1, x2, x3, x, y, z) = polynomial_ring(QQ, ["x1", "x2", "x3", "x", "y", "z"])
(Multivariate polynomial ring in 6 variables over QQ, QQMPolyRingElem[x1, x2, x3, x, y, z])

julia> p = x^3 + 7*x1*x2^5*x^2*z^2 + x1^3*(x2 + x3)^9*x*z^4 - y^2 - 13*x3^3*x*y*z - x1^2*x2^4*x3^3*y*z^3;

julia> h = hypersurface_model(B2, fiber_ambient_space, [D1, D2], p, completeness_check = false)
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(cox_ring(B2));

julia> a1 = 13 * x3^3;

julia> a2 = 7 * x1 * x2^5;

julia> a3 = x1^2 * x2^4 * x3^3;

julia> a4 = x1^3 * (x2 + x3)^9;

julia> a6 = zero(cox_ring(B2));

julia> t = global_tate_model(B2, [a1, a2, a3, a4, a6])
Global Tate model over a concrete base

julia> set_global_tate_model(h, t)

julia> global_tate_model(h) === t
true
```
"""
function global_tate_model(h::HypersurfaceModel)
  @req has_attribute(h, :global_tate_model) "No corresponding global Tate model is known"
  t = get_attribute(h, :global_tate_model)
  if t isa GlobalTateModel
    return t
  end
  @req t isa String "Internal inconsistency encountered"
  directory = joinpath(dirname(@__DIR__), "LiteratureModels/")
  model_indices = JSON.parsefile(directory * "model_indices.json")
  if is_base_space_fully_specified(h)
    t_model = literature_model(parse(Int, model_indices[t]), base_space = base_space(h), defining_classes = defining_classes(h), completeness_check = false)
  else
    t_model = literature_model(parse(Int, model_indices[t]))
  end
  set_global_tate_model(h, t_model)
  return t_model
end


############################################################################################################
############################################################################################################
# 3: Attributes that rest on the corresponding Weierstrass model and (currently) only work in toric settings
############################################################################################################
############################################################################################################


#####################################################
# 3.1 Calabi-Yau hypersurface
#####################################################

@doc raw"""
    calabi_yau_hypersurface(h::HypersurfaceModel)

Returns the Calabi–Yau hypersurface that defines the hypersurface model
as a closed subvariety of its toric ambient space.

```jldoctest
julia> b = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> D1 = 2 * anticanonical_divisor_class(b)
Divisor class on a normal toric variety

julia> D2 = 3 * anticanonical_divisor_class(b)
Divisor class on a normal toric variety

julia> new_gens = string.(vcat(gens(cox_ring(b)), gens(cox_ring(fiber_ambient_space))));

julia> ambient_ring, (x1, x2, x3, x, y, z) = polynomial_ring(QQ, new_gens, cached=false)
(Multivariate polynomial ring in 6 variables over QQ, QQMPolyRingElem[x1, x2, x3, x, y, z])

julia> p = x^3 - y^2 + x1^12 * x * z^4 + x2^18 * z^6 + 13 * x3^3*x*y*z
x1^12*x*z^4 + x2^18*z^6 + 13*x3^3*x*y*z + x^3 - y^2

julia> h = hypersurface_model(b, fiber_ambient_space, [D1, D2], p; completeness_check = false)
Hypersurface model over a concrete base

julia> calabi_yau_hypersurface(h)
Closed subvariety of a normal toric variety
```
"""
@attr ClosedSubvarietyOfToricVariety function calabi_yau_hypersurface(h::HypersurfaceModel)
  @req base_space(h) isa NormalToricVariety "Calabi-Yau hypersurface currently only supported for toric varieties as base space"
  is_base_space_fully_specified(h) || @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.\n"
  return closed_subvariety_of_toric_variety(ambient_space(h), [hypersurface_equation(h)])
end


#####################################################
# 3.2 Discriminant and singular loci
#####################################################

@doc raw"""
    discriminant(h::HypersurfaceModel)

Returns the discriminant ``\Delta = 4f^3 + 27g^2`` of the Weierstrass model associated
with the given hypersurface model.

Raises an error if no such Weierstrass model is known.

In the example below, we construct a hypersurface model and its corresponding Weierstrass
model (see [BMT25](@cite BMT25) for background), in order to demonstrate this functionality.

```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> D1 = 2 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> D2 = 3 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> amb_ring, (x1, x2, x3, x, y, z) = polynomial_ring(QQ, ["x1", "x2", "x3", "x", "y", "z"])
(Multivariate polynomial ring in 6 variables over QQ, QQMPolyRingElem[x1, x2, x3, x, y, z])

julia> p = x^3 + 7*x1*x2^5*x^2*z^2 + x1^3*(x2 + x3)^9*x*z^4 - y^2 - 13*x3^3*x*y*z - x1^2*x2^4*x3^3*y*z^3;

julia> h = hypersurface_model(B2, fiber_ambient_space, [D1, D2], p, completeness_check = false)
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(cox_ring(B2));

julia> weier_f = 1//48*(-(28*x1*x2^5 + 169*x3^6)^2 + 24*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> weier_g = 1//864*(216*x1^4*x2^8*x3^6 + (28*x1*x2^5 + 169*x3^6)^3 - 36*(28*x1*x2^5 + 169*x3^6)*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> w = weierstrass_model(B2, weier_f, weier_g; completeness_check = false)
Weierstrass model over a concrete base

julia> set_weierstrass_model(h, w)

julia> degree(discriminant(h))
Abelian group element [36]
```
"""
@attr MPolyRingElem function discriminant(h::HypersurfaceModel)
  @req base_space(h) isa NormalToricVariety "Discriminant currently only supported for toric varieties as base space"
  @req has_attribute(h, :weierstrass_model) "No corresponding Weierstrass model is known"
  return discriminant(weierstrass_model(h))
end


@doc raw"""
    singular_loci(h::HypersurfaceModel)

Returns the singular loci of the Weierstrass model equivalent to the given hypersurface model,
along with the order of vanishing of ``(f, g, \Delta)`` at each locus and the corresponding
refined Tate fiber type. See [singular_loci(w::WeierstrassModel)](@ref) for more details.

Raises an error if no such Weierstrass model is known.

In the example below, we construct a hypersurface model and its corresponding Weierstrass
model (see [BMT25](@cite BMT25) for background), in order to demonstrate this functionality.

!!! warning
    The classification of singularities is performed using a Monte Carlo algorithm, involving randomized sampling.
    While reliable in practice, this probabilistic method may occasionally yield non-deterministic results.

```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> D1 = 2 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> D2 = 3 * anticanonical_divisor_class(B2)
Divisor class on a normal toric variety

julia> amb_ring, (x1, x2, x3, x, y, z) = polynomial_ring(QQ, ["x1", "x2", "x3", "x", "y", "z"])
(Multivariate polynomial ring in 6 variables over QQ, QQMPolyRingElem[x1, x2, x3, x, y, z])

julia> p = x^3 + 7*x1*x2^5*x^2*z^2 + x1^3*(x2 + x3)^9*x*z^4 - y^2 - 13*x3^3*x*y*z - x1^2*x2^4*x3^3*y*z^3;

julia> h = hypersurface_model(B2, fiber_ambient_space, [D1, D2], p, completeness_check = false)
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(cox_ring(B2));

julia> weier_f = 1//48*(-(28*x1*x2^5 + 169*x3^6)^2 + 24*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> weier_g = 1//864*(216*x1^4*x2^8*x3^6 + (28*x1*x2^5 + 169*x3^6)^3 - 36*(28*x1*x2^5 + 169*x3^6)*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> w = weierstrass_model(B2, weier_f, weier_g; completeness_check = false)
Weierstrass model over a concrete base

julia> set_weierstrass_model(h, w)

julia> length(singular_loci(h))
2
```
"""
@attr Vector{<:Tuple{<:MPolyIdeal{<:MPolyRingElem}, Tuple{Int64, Int64, Int64}, String}} function singular_loci(h::HypersurfaceModel)
  @req base_space(h) isa NormalToricVariety "Singular loci currently only supported for toric varieties as base space"
  @req has_attribute(h, :weierstrass_model) || has_attribute(h, :global_tate_model) "No corresponding Weierstrass model or global Tate model is known"
  return has_attribute(h, :weierstrass_model) ? singular_loci(weierstrass_model(h)) : singular_loci(global_tate_model(h))
end
