#####################################################
# 1 Basic properties
#####################################################

@doc raw"""
    is_well_quantized(gf::G4Flux)

Check whether the given ``G_4``-flux candidate satisfies necessary consistency conditions
for flux quantization as formulated in [Wit97](@cite):

$G_4 + \frac{1}{2} c_2(\widehat{Y}_4) \in H^{(2,2)}(\widehat{Y}_4, \mathbb{Z})\,.$

Since verifying this integrality condition is generally very difficult, this method performs
a series of simpler checks. The flux candidate is modelled by ``g \in H^{2,2}(X_\Sigma, \mathbb{Q})``.
This method evaluates

$\int_{X_\Sigma} \left(g + \frac{1}{2} \hat{c}_2 \right) \wedge [H] \wedge [D_i] \wedge [D_j]$

for all pairs of toric divisors ``D_i``, ``D_j`` in the ambient variety, where ``[H]`` denotes
the class of the hypersurface divisor defining ``\widehat{Y}_4`` and ``\hat{c}_2 \in H^{2,2}(X_\Sigma, \mathbb{Q})``
restricts to ``c_2(\widehat{Y}_4)`` on the hypersurface.

If all these integrals evaluate to integers, this method returns `true`; otherwise, it returns `false`.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> g4_class = cohomology_class(anticanonical_divisor_class(ambient_space(qsm_model)), quick = true)^2;

julia> g4 = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> is_well_quantized(g4)
true

julia> g4
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed
```
"""
@attr Bool function is_well_quantized(g4::G4Flux)
  m = model(g4)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Elementary quantization checks for  G4-fluxes only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Elementary quantization checks for G4-flux currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Elementary quantization checks for G4-flux currently supported only for toric ambient space"

  # Compute the cohomology class corresponding to the hypersurface equation
  cy = polynomial(cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))))

  # Now check quantization condition G4 + 1/2 c2 is integral.
  c_ds = [polynomial(cohomology_class(d)) for d in torusinvariant_prime_divisors(ambient_space(m))]

  # explicitly switched off an expensive test in the following line
  twist_g4 = polynomial(cohomology_class(g4) + 1//2 * chern_class(m, 2; check = false))

  # now execute elementary checks of the quantization condition
  for i in 1:length(c_ds)
    for j in i:length(c_ds)
      numb = integrate(cohomology_class(ambient_space(m), twist_g4 * c_ds[i] * c_ds[j] * cy); check = false)
      !is_integer(numb) && return false
    end
  end
  return true
end


@doc raw"""
    passes_transversality_checks(gf::G4Flux)

Check whether the ``G_4``-flux satisfies the transversality conditions
(cf. [Wei18](@cite)). Return `true` if all conditions are met, otherwise `false`.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> divs = torusinvariant_prime_divisors(ambient_space(qsm_model));

julia> e1 = cohomology_class(divs[35]);e2 = cohomology_class(divs[32]);e4 = cohomology_class(divs[34]);

julia> u = cohomology_class(divs[33]);v = cohomology_class(divs[30]);pb_Kbar = cohomology_class(sum([divs[k] for k in 1:29]));

julia> g4_class = (-3) // kbar3(qsm_model) * (5 * e1 * e4 + pb_Kbar * (-3 * e1 - 2 * e2 - 6 * e4 + pb_Kbar - 4 * u + v));

julia> g4 = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> passes_transversality_checks(g4)
true

julia> g4
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed
```
"""
@attr Bool function passes_transversality_checks(g4::G4Flux)
  m = model(g4)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Transversality checks supported only for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Transversality checks supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Transversality checks supported only for toric ambient space"
  @req has_attribute(m, :zero_section_class) "Transversality checks require zero section class"
  
  # Compute the cohomology class corresponding to the hypersurface equation
  cy = polynomial(cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))))
   
  n = ngens(coordinate_ring(base_space(m)))
  c_ds = [polynomial(cohomology_class(d)) for d in torusinvariant_prime_divisors(ambient_space(m))[1:n]]
  zero_sec = zero_section_class(m)

  # now execute checks to verify if the transversality conditions are satisfied
  for i in 1:n
    numb = integrate(cohomology_class(ambient_space(m), polynomial(cohomology_class(g4)) * c_ds[i] * cy) * zero_sec; check = false)
    numb!=0 && return false
  end
  for i in 1:n
    for j in i:n
      numb = integrate(cohomology_class(ambient_space(m), polynomial(cohomology_class(g4)) * c_ds[i] * c_ds[j] * cy); check = false)
      numb!=0 && return false
    end
  end
  return true
end


@doc raw"""
    passes_tadpole_cancellation_check(gf::G4Flux)

Check whether the given ``G_4``-flux satisfies the D3-tadpole cancellation condition. This
amounts to verifying that

$\frac{\chi(\widehat{Y}_4)}{24} - \frac{1}{2} \int_{\widehat{Y}_4} G_4 \wedge G_4$

is a non-negative integer.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> divs = torusinvariant_prime_divisors(ambient_space(qsm_model));

julia> e1 = cohomology_class(divs[35]);e2 = cohomology_class(divs[32]);e4 = cohomology_class(divs[34]);

julia> u = cohomology_class(divs[33]);v = cohomology_class(divs[30]);pb_Kbar = cohomology_class(sum([divs[k] for k in 1:29]));

julia> g4_class = (-3) // kbar3(qsm_model) * (5 * e1 * e4 + pb_Kbar * (-3 * e1 - 2 * e2 - 6 * e4 + pb_Kbar - 4 * u + v));

julia> g4 = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> passes_tadpole_cancellation_check(g4)
true

julia> g4
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: satisfied
```
"""
@attr Bool function passes_tadpole_cancellation_check(g4::G4Flux)
  m = model(g4)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Tadpole cancellation checks for G4-fluxes only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Tadpole cancellation checks for G4-flux currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Tadpole cancellation checks for G4-flux currently supported only for toric ambient space"
  numb = d3_tadpole_constraint(g4, check = false)
  return numb >= 0 && is_integer(numb)
end


@doc raw"""
    breaks_non_abelian_gauge_group(gf::G4Flux)

Check whether the given ``G_4``-flux candidate breaks any non-abelian gauge
symmetries. Return `true` if any breaking occurs, and `false` otherwise.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> divs = torusinvariant_prime_divisors(ambient_space(qsm_model));

julia> e1 = cohomology_class(divs[35]);e2 = cohomology_class(divs[32]);e4 = cohomology_class(divs[34]);

julia> u = cohomology_class(divs[33]);v = cohomology_class(divs[30]);pb_Kbar = cohomology_class(sum([divs[k] for k in 1:29]));

julia> g4_class = (-3) // kbar3(qsm_model) * (5 * e1 * e4 + pb_Kbar * (-3 * e1 - 2 * e2 - 6 * e4 + pb_Kbar - 4 * u + v));

julia> g4 = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> breaks_non_abelian_gauge_group(g4)
false

julia> g4
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: unbroken
  - Tadpole cancellation check: not computed
```
"""
@attr Bool function breaks_non_abelian_gauge_group(g4::G4Flux)
  m = model(g4)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Checks for breaking non-abelian gauge group factors only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Checks for breaking non-abelian gauge group factors currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Checks for breaking non-abelian gauge group factors currently supported only for toric ambient space"
  
  # Compute the cohomology class corresponding to the hypersurface equation
  cy = polynomial(cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))))

  # Identify the cohomology classes of all base divisors
  n = ngens(coordinate_ring(base_space(m)))
  c_ds = [polynomial(cohomology_class(d)) for d in torusinvariant_prime_divisors(ambient_space(m))[1:n]]

  # Identify the cohomology classes of all exceptional divisors
  gS = gens(coordinate_ring(ambient_space(m)))
  exceptional_divisor_positions = exceptional_divisor_indices(m)
  exceptional_divisors = torusinvariant_prime_divisors(ambient_space(m))[exceptional_divisor_positions]
  c_ei = [polynomial(cohomology_class(d)) for d in exceptional_divisors]

  # now execute the checks if any non-abelian gauge group factor is broken
  for i in 1:n
    for j in 1:length(exceptional_divisors)
      numb = integrate(cohomology_class(ambient_space(m), polynomial(cohomology_class(g4)) * c_ds[i] * c_ei[j] * cy); check = false)
      numb!=0 && return true
    end
  end
  return false
end
