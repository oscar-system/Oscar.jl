@doc raw"""
    hypersurface_equation(h::HypersurfaceModel)

Return the hypersurface equation of the hypersurface model.

# Examples
```jldoctest
julia> b = projective_space(NormalToricVariety, 2);

julia> kb = anticanonical_divisor_class(b);

julia> fiber_ambient = weighted_projective_space(NormalToricVariety, [2, 3, 1]);

julia> set_coordinate_names(fiber_ambient, ["x", "y", "z"]);

julia> p = "x^3 - y^2 + x1^12 * x * z^4 + x2^18 * z^6 + 13 * x3^3 * x * y * z";

julia> h = hypersurface_model(b, fiber_ambient, [2 * kb, 3 * kb], p; completeness_check=false)
Hypersurface model over a concrete base

julia> hypersurface_equation(h)
x1^12*x*z^4 + x2^18*z^6 + 13*x3^3*x*y*z + x^3 - y^2
```
"""
hypersurface_equation(h::HypersurfaceModel) = h.hypersurface_equation

@doc raw"""
    hypersurface_equation_parametrization(h::HypersurfaceModel)

Return the parametrization of the hypersurface equation by the model sections.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> using Random;

julia> h = literature_model(arxiv_id = "1208.2695", equation = "B.5", rng = Random.Xoshiro(1234))
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
hypersurface_equation_parametrization(h::HypersurfaceModel) =
  h.hypersurface_equation_parametrization

@doc raw"""
    weierstrass_model(h::HypersurfaceModel)

Return the Weierstrass model corresponding to the given hypersurface model, if known.

If needed, this function constructs a literature model, the process of which may include the
creation of generic sections. The random source used in the creation of said generic sections
can be set with the optional argument `rng`.

In the example below, we construct a hypersurface model and its corresponding Weierstrass
model (see [BMT25](@cite BMT25) for background), and then establish the relationship
between the two models.

# Examples
```jldoctest
julia> b = projective_space(NormalToricVariety, 2);

julia> kb = anticanonical_divisor_class(b);

julia> fiber_ambient = weighted_projective_space(NormalToricVariety, [2, 3, 1]);

julia> set_coordinate_names(fiber_ambient, ["x", "y", "z"]);

julia> p = "x^3 + 7*x1*x2^5*x^2*z^2 + x1^3*(x2 + x3)^9*x*z^4 - y^2 - 13*x3^3*x*y*z - x1^2*x2^4*x3^3*y*z^3";

julia> h = hypersurface_model(b, fiber_ambient, [2 * kb, 3 * kb], p; completeness_check=false)
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(coordinate_ring(b));

julia> weier_f = 1//48*(-(28*x1*x2^5 + 169*x3^6)^2 + 24*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> weier_g = 1//864*(216*x1^4*x2^8*x3^6 + (28*x1*x2^5 + 169*x3^6)^3 - 36*(28*x1*x2^5 + 169*x3^6)*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> w = weierstrass_model(b, weier_f, weier_g; completeness_check = false)
Weierstrass model over a concrete base

julia> set_weierstrass_model(h, w)

julia> weierstrass_model(h) === w
true
```
"""
function weierstrass_model(h::HypersurfaceModel; rng::AbstractRNG=Random.default_rng())
  @req has_attribute(h, :weierstrass_model) "No corresponding Weierstrass model is known"
  w = get_attribute(h, :weierstrass_model)
  if w isa WeierstrassModel
    return w
  end
  @req w isa String "Internal inconsistency encountered"
  directory = joinpath(dirname(@__DIR__), "LiteratureModels/")
  model_indices = JSON.parsefile(directory * "model_indices.json")
  if is_base_space_fully_specified(h)
    w_model = literature_model(
      parse(Int, model_indices[w]);
      base_space=base_space(h),
      defining_classes=defining_classes(h),
      completeness_check=false,
      rng=rng,
    )
  else
    w_model = literature_model(parse(Int, model_indices[w]); rng=rng)
  end
  set_weierstrass_model(h, w_model)
  return w_model
end

@doc raw"""
    global_tate_model(h::HypersurfaceModel)

Return the global Tate model corresponding to the given hypersurface model, if known.

If needed, this function constructs a literature model, the process of which may include the
creation of generic sections. The random source used in the creation of said generic sections
can be set with the optional argument `rng`.

In the example below, we construct a hypersurface model and its corresponding global
Tate model (see [BMT25](@cite BMT25) for background), and then establish the relationship
between the two models.

# Examples
```jldoctest
julia> b = projective_space(NormalToricVariety, 2);

julia> kb = anticanonical_divisor_class(b);

julia> fiber_ambient = weighted_projective_space(NormalToricVariety, [2, 3, 1]);

julia> set_coordinate_names(fiber_ambient, ["x", "y", "z"]);

julia> p = "x^3 + 7*x1*x2^5*x^2*z^2 + x1^3*(x2 + x3)^9*x*z^4 - y^2 - 13*x3^3*x*y*z - x1^2*x2^4*x3^3*y*z^3";

julia> h = hypersurface_model(b, fiber_ambient, [2 * kb, 3 * kb], p; completeness_check=false)
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(coordinate_ring(b));

julia> a1 = 13 * x3^3;

julia> a2 = 7 * x1 * x2^5;

julia> a3 = x1^2 * x2^4 * x3^3;

julia> a4 = x1^3 * (x2 + x3)^9;

julia> a6 = zero(coordinate_ring(b));

julia> t = global_tate_model(b, [a1, a2, a3, a4, a6])
Global Tate model over a concrete base

julia> set_global_tate_model(h, t)

julia> global_tate_model(h) === t
true
```
"""
function global_tate_model(h::HypersurfaceModel; rng::AbstractRNG=Random.default_rng())
  @req has_attribute(h, :global_tate_model) "No corresponding global Tate model is known"
  t = get_attribute(h, :global_tate_model)
  if t isa GlobalTateModel
    return t
  end
  @req t isa String "Internal inconsistency encountered"
  directory = joinpath(dirname(@__DIR__), "LiteratureModels/")
  model_indices = JSON.parsefile(directory * "model_indices.json")
  if is_base_space_fully_specified(h)
    t_model = literature_model(
      parse(Int, model_indices[t]);
      base_space=base_space(h),
      defining_classes=defining_classes(h),
      completeness_check=false,
      rng=rng,
    )
  else
    t_model = literature_model(parse(Int, model_indices[t]); rng=rng)
  end
  set_global_tate_model(h, t_model)
  return t_model
end

@doc raw"""
    discriminant(h::HypersurfaceModel)

Return the discriminant ``\Delta = 4f^3 + 27g^2`` of the Weierstrass model associated
with the given hypersurface model.

Raises an error if no such Weierstrass model is known.

In the example below, we construct a hypersurface model and its corresponding Weierstrass
model (see [BMT25](@cite BMT25) for background), in order to demonstrate this functionality.

# Examples
```jldoctest
julia> b = projective_space(NormalToricVariety, 2);

julia> kb = anticanonical_divisor_class(b);

julia> fiber_ambient = weighted_projective_space(NormalToricVariety, [2, 3, 1]);

julia> set_coordinate_names(fiber_ambient, ["x", "y", "z"]);

julia> p = "x^3 + 7*x1*x2^5*x^2*z^2 + x1^3*(x2 + x3)^9*x*z^4 - y^2 - 13*x3^3*x*y*z - x1^2*x2^4*x3^3*y*z^3";

julia> h = hypersurface_model(b, fiber_ambient, [2 * kb, 3 * kb], p; completeness_check=false)
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(coordinate_ring(b));

julia> weier_f = 1//48*(-(28*x1*x2^5 + 169*x3^6)^2 + 24*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> weier_g = 1//864*(216*x1^4*x2^8*x3^6 + (28*x1*x2^5 + 169*x3^6)^3 - 36*(28*x1*x2^5 + 169*x3^6)*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> w = weierstrass_model(b, weier_f, weier_g; completeness_check = false)
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

Return the singular loci of the Weierstrass model equivalent to the given hypersurface model,
along with the order of vanishing of ``(f, g, \Delta)`` at each locus and the corresponding
refined Tate fiber type. See [`singular_loci(w::WeierstrassModel)`](@ref) for more details.

Raises an error if no such Weierstrass model is known.

In the example below, we construct a hypersurface model and its corresponding Weierstrass
model (see [BMT25](@cite BMT25) for background), in order to demonstrate this functionality.

!!! warning
    The classification of singularities is performed using a Monte Carlo algorithm, involving randomized sampling.
    While reliable in practice, this probabilistic method may occasionally yield non-deterministic results.
    The random source can be set with the optional argument `rng`.

# Examples
```jldoctest
julia> b = projective_space(NormalToricVariety, 2);

julia> kb = anticanonical_divisor_class(b);

julia> fiber_ambient = weighted_projective_space(NormalToricVariety, [2, 3, 1]);

julia> set_coordinate_names(fiber_ambient, ["x", "y", "z"]);

julia> p = "x^3 + 7*x1*x2^5*x^2*z^2 + x1^3*(x2 + x3)^9*x*z^4 - y^2 - 13*x3^3*x*y*z - x1^2*x2^4*x3^3*y*z^3";

julia> h = hypersurface_model(b, fiber_ambient, [2 * kb, 3 * kb], p; completeness_check=false)
Hypersurface model over a concrete base

julia> x1, x2, x3 = gens(coordinate_ring(b));

julia> weier_f = 1//48*(-(28*x1*x2^5 + 169*x3^6)^2 + 24*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> weier_g = 1//864*(216*x1^4*x2^8*x3^6 + (28*x1*x2^5 + 169*x3^6)^3 - 36*(28*x1*x2^5 + 169*x3^6)*(2*x1^3*(x2 + x3)^9 + 13*x1^2*x2^4*x3^6));

julia> w = weierstrass_model(b, weier_f, weier_g; completeness_check = false)
Weierstrass model over a concrete base

julia> set_weierstrass_model(h, w)

julia> using Random;

julia> length(singular_loci(h; rng = Random.Xoshiro(1234)))
2
```
"""
@attr Vector{<:Tuple{<:MPolyIdeal{<:MPolyRingElem},Tuple{Int64,Int64,Int64},String}} function singular_loci(
  h::HypersurfaceModel; rng::AbstractRNG=Random.default_rng()
)
  @req base_space(h) isa NormalToricVariety "Singular loci currently only supported for toric varieties as base space"
  @req has_attribute(h, :weierstrass_model) || has_attribute(h, :global_tate_model) "No corresponding Weierstrass model or global Tate model is known"
  return if has_attribute(h, :weierstrass_model)
    singular_loci(weierstrass_model(h); rng)
  else
    singular_loci(global_tate_model(h); rng)
  end
end
