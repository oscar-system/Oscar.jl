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

It is worth mentioning that currently (July 2024), we only
support this check for $G_4$-fluxes defined on Weierstrass,
global Tate and hypersurface models. If this condition is not
met, this method will return an error.

```jldoctest
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
  if m isa WeierstrassModel
    cl = toric_divisor_class(ambient_space(m), degree(weierstrass_polynomial(m)))
  end
  if m isa GlobalTateModel
    cl = toric_divisor_class(ambient_space(m), degree(tate_polynomial(m)))
  end
  if m isa HypersurfaceModel
    cl = toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))
  end
  cy = polynomial(cohomology_class(cl))

  # Now check quantization condition G4 + 1/2 c2 is integral.
  c_ds = [polynomial(cohomology_class(d)) for d in torusinvariant_prime_divisors(ambient_space(m))]

  # explicitly switched off an expensive test in the following line
  twist_g4 = polynomial(cohomology_class(g4) + 1//2 * chern_class_c2(m; check = false))

  # now execute elementary checks of the quantization condition
  for i in 1:length(c_ds)
    for j in i:length(c_ds)
      numb = integrate(cohomology_class(ambient_space(m), twist_g4 * cy * c_ds[i] * c_ds[j]); check = false)
      !is_integer(numb) && return false
    end
  end
  return true
end
