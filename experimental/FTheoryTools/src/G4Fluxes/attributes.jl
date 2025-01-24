#####################################################
# 1 Basic attributes
#####################################################

@doc raw"""
    model(gf::G4Flux)

Return the F-theory model for which this $G_4$-flux candidate is defined.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> cohomology_ring(ambient_space(qsm_model), check = false);

julia> g4_class = cohomology_class(anticanonical_divisor_class(ambient_space(qsm_model)))^2;

julia> g4f = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Tadpole cancellation check: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed  

julia> model(g4f)
Hypersurface model over a concrete base
```
"""
model(gf::G4Flux) = gf.model


@doc raw"""
    cohomology_class(gf::G4Flux)

Return the cohomology class which defines the $G_4$-flux candidate.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> cohomology_ring(ambient_space(qsm_model), check = false);

julia> g4_class = cohomology_class(anticanonical_divisor_class(ambient_space(qsm_model)))^2;

julia> g4f = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Tadpole cancellation check: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> cohomology_class(g4f) == g4_class
true
```
"""
cohomology_class(gf::G4Flux) = gf.class



#####################################################
# 2 Compute the D3-tadpole constraint
#####################################################

@doc raw"""
    d3_tadpole_constraint(gf::G4Flux; check::Bool = true)

Return the d3-tapdole constraint of a G4-flux, that is compute the quantity
$- \frac{1}{2} \cdot G_4^2 + \frac{1}{24} \cdot \chi(\widehat{Y}_4) \stackrel{!}{\geq} 0$.

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
  - Tadpole cancellation check: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> d3_tadpole_constraint(g4, check = false)
12
```
"""
function d3_tadpole_constraint(gf::G4Flux; check::Bool = true)
  if has_attribute(gf, :d3_tadpole_constraint)
    return get_attribute(gf, :d3_tadpole_constraint)::QQFieldElem
  end
  m = model(gf)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Tadpole cancellation checks for G4-fluxes only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Tadpole cancellation checks for G4-flux currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Tadpole cancellation checks for G4-flux currently supported only for toric ambient space"
  if check
    @req is_complete(ambient_space(m)) "Computation of D3-tadpole constraint only supported for complete toric ambient spaces"
    @req is_simplicial(ambient_space(m)) "Computation of D3-tadpole constraint only supported for simplicial toric ambient space"
  end
  cy = polynomial(cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))))
  numb = QQ(euler_characteristic(m; check = check)//24 - 1//2*integrate(cohomology_class(ambient_space(m), polynomial(cohomology_class(gf)) * polynomial(cohomology_class(gf)) * cy); check = check))
  set_attribute!(gf, :d3_tadpole_constraint, numb)
  set_attribute!(gf, :passes_tadpole_cancellation_check, (numb >= 0 && is_integer(numb)))
  return numb::QQFieldElem
end
