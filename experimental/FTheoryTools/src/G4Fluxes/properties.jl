#####################################################
# 1 Basic properties
#####################################################

@doc raw"""
    passes_elementary_quantization_checks(gf::G4Flux)

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
G4-flux candidate lacking elementary quantization checks

julia> passes_elementary_quantization_checks(g4)
true
```
"""
@attr Bool function passes_elementary_quantization_checks(g4::G4Flux)
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
    passes_verticality_checks(gf::G4Flux)

G4-fluxes are subject to verticality conditions described first in [GH12](@cite) and in more detail in [Wei18](@cite).
It is hard to verify that these condition are met. However,
we can execute a number of simple consistency checks, by
verifying that $\int_{Y}{G_4 \wedge [D_1] \wedge [zero section]} = 0$ and $\int_{Y}{G_4 \wedge [D_1] \wedge [D_2]} = 0$
for all toric base divisors $D_1$ and $D_2$. If all of these
simple consistency checks are met, this method will return
`true` and otherwise `false`

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> divs = torusinvariant_prime_divisors(ambient_space(qsm_model));

julia> e1 = cohomology_class(divs[35]);e2 = cohomology_class(divs[32]);e4 = cohomology_class(divs[34]);

julia> u = cohomology_class(divs[33]);v = cohomology_class(divs[30]);pb_Kbar = cohomology_class(sum([divs[k] for k in 1:29]));

julia> g4_class = (-3) // kbar3(qsm_model) * (5 * e1 * e4 + pb_Kbar * (-3 * e1 - 2 * e2 - 6 * e4 + pb_Kbar - 4 * u + v));

julia> g4 = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate lacking elementary quantization checks

julia> passes_verticality_checks(g4)
true
```
"""
@attr Bool function passes_verticality_checks(g4::G4Flux)
  m = model(g4)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Tadpole cancellation checks for  G4-fluxes only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Tadpole cancellation checks for G4-flux currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Tadpole cancellation checks for G4-flux currently supported only for toric ambient space"
  @req has_zero_section_class(m) "For verticality checks, a model zero section class needs to be specified"
  
  # Compute the cohomology class corresponding to the hypersurface equation
  cy = polynomial(cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))))
   
  n = length(gens(cox_ring(base_space(m))))
  c_ds = [polynomial(cohomology_class(d)) for d in torusinvariant_prime_divisors(ambient_space(m))[1:n]]
  zero_sec = zero_section_class(m)

  # now execute verticality checks
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
G4-flux candidate lacking elementary quantization checks

julia> passes_tadpole_cancellation_check(g4)
true
```
"""
@attr Bool function passes_tadpole_cancellation_check(g4::G4Flux)
  m = model(g4)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Tadpole cancellation checks for  G4-fluxes only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Tadpole cancellation checks for G4-flux currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Tadpole cancellation checks for G4-flux currently supported only for toric ambient space"
  
  # Compute the cohomology class corresponding to the hypersurface equation
  cy = polynomial(cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))))
  
  # Now check if the D3-tadpole cancellation condition holds
  numb = euler_characteristic(m; check = false)/24 - 1/2*integrate(cohomology_class(ambient_space(m), polynomial(cohomology_class(g4)) * polynomial(cohomology_class(g4)) * cy); check = false)
  if numb > 0 && is_integer(numb)
    return true
  end
  return false
end
