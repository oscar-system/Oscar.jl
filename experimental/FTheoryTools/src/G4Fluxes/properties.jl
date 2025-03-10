#####################################################
# 1 Basic properties
#####################################################

@doc raw"""
    is_well_quantized(gf::G4Flux)

G4-fluxes are subject to the quantization condition
[Wit97](@cite) $G_4 + \frac{1}{2} c_2(Y) \in H^{(2,2)}(Y, \mathbb{Z})$.
It is hard to verify that this condition is met. However,
we can execute a number of simple consistency checks, by
verifying that $\int_{Y}{G_4 \wedge [D_1] \wedge [D_2]} \in \mathbb{Z}$
for any two toric divisors $D_1$, $D_2$. If all of these
simple consistency checks are met, this method will return
`true` and otherwise `false`.

It is worth mentioning that currently (August 2024), we only
support this check for $G_4$-fluxes defined on Weierstrass,
global Tate and hypersurface models. If this condition is not
met, this method will return an error.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> cohomology_ring(ambient_space(qsm_model), check = false);

julia> g4_class = cohomology_class(anticanonical_divisor_class(ambient_space(qsm_model)))^2;

julia> g4 = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

julia> is_well_quantized(g4)
true

julia> g4
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed
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

G4-fluxes are subject to transversality conditions (cf. [Wei18](@cite)).
If these conditions are met, this method will return `true` and otherwise `false`

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
  - Tadpole cancellation check: not executed

julia> passes_transversality_checks(g4)
true

julia> g4
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed
```
"""
@attr Bool function passes_transversality_checks(g4::G4Flux)
  m = model(g4)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Transversality checks supported only for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Transversality checks supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Transversality checks supported only for toric ambient space"
  @req has_zero_section_class(m) "Transversality checks require zero section class"
  
  # Compute the cohomology class corresponding to the hypersurface equation
  cy = polynomial(cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))))
   
  n = ngens(cox_ring(base_space(m)))
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

G4-fluxes are subject to the D3-tadpole cancellation condition described in [Wei18](@cite).
This check verifies that $euler_characteristic(Y)/24 - 1/2 * \int_{Y}{G_4 \wedge G_4}$ is a non-negative integer.

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
  - Tadpole cancellation check: not executed

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

G4-fluxes may break the non-abelian gauge group (cf. [Wei18](@cite)).
This function verifies if this is the case for the given G4-flux.
If it does not break any non-abelian gauge factor, we return 
`true` and otherwise `false`

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
  - Tadpole cancellation check: not executed

julia> breaks_non_abelian_gauge_group(g4)
false

julia> g4
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: not broken
  - Tadpole cancellation check: not executed
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
  n = ngens(cox_ring(base_space(m)))
  c_ds = [polynomial(cohomology_class(d)) for d in torusinvariant_prime_divisors(ambient_space(m))[1:n]]

  # Identify the cohomology classes of all exceptional divisors
  gS = gens(cox_ring(ambient_space(m)))
  exceptional_divisor_positions = findall(x -> occursin(r"^e\d+$", x), string.(symbols(cox_ring(ambient_space(m)))))
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
