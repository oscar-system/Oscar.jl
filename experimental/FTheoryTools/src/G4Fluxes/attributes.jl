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
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed  
  - Tadpole cancellation check: not executed

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
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

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
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

julia> d3_tadpole_constraint(g4, check = false)
12

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gfs = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: failed
  - Non-abelian gauge group: broken
  - Tadpole constraint: not analyzed

julia> g4_2 = random_flux_instance(gfs, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

julia> dv2 = d3_tadpole_constraint(g4_2, check = false);

julia> int_comb = integral_coefficients(g4_2);

julia> rat_comb = rational_coefficients(g4_2);

julia> g4_3 = flux_instance(gfs, int_comb, rat_comb, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

julia> d3_tadpole_constraint(gfs, check = false);

julia> dv3 = d3_tadpole_constraint(g4_3, check = false);

julia> dv2 == dv3
true
```
"""
@attr QQFieldElem function d3_tadpole_constraint(gf::G4Flux; check::Bool = true)
  m = model(gf)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Tadpole cancellation checks for G4-fluxes only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "Tadpole cancellation checks for G4-flux currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Tadpole cancellation checks for G4-flux currently supported only for toric ambient space"
  if check
    @req is_complete(ambient_space(m)) "Computation of D3-tadpole constraint only supported for complete toric ambient spaces"
    @req is_simplicial(ambient_space(m)) "Computation of D3-tadpole constraint only supported for simplicial toric ambient space"
  end
  # Opt for (potentially?) quicker algorithm when possible
  if has_attribute(gf, :int_combination) && has_attribute(gf, :rat_combination)
    gfs = g4_flux_family(gf, check = check)
    if has_attribute(gfs, :d3_tadpole_constraint)
      values_to_evaluate_at = vcat(matrix(QQ, get_attribute(gf, :int_combination)), get_attribute(gf, :rat_combination))
      values_to_evaluate_at = values_to_evaluate_at[:, 1]
      return evaluate(d3_tadpole_constraint(gfs), values_to_evaluate_at)
    end
  end
  cy = polynomial(cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))))
  numb = QQ(euler_characteristic(m; check = check)//24 - 1//2*integrate(cohomology_class(ambient_space(m), polynomial(cohomology_class(gf)) * polynomial(cohomology_class(gf)) * cy); check = check))
  set_attribute!(gf, :passes_tadpole_cancellation_check, (numb >= 0 && is_integer(numb)))
  return numb::QQFieldElem
end


#####################################################
# 3 The "position" of a flux within a G4-flux family
#####################################################

@doc raw"""
    g4_flux_family(gf::G4Flux; check::Bool = true)

Return the family of $G_4$-fluxes that possesses, such
that all fluxes in this family share the following properties
with the given $G_4$-flux: Verticality and breaking of the
non-Abelian gauge group.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gfs = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: failed
  - Non-abelian gauge group: broken
  - Tadpole constraint: not analyzed

julia> g4 = random_flux_instance(gfs, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

julia> g4_flux_family(g4, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: failed
  - Non-abelian gauge group: broken
  - Tadpole constraint: not analyzed
```
"""
@attr FamilyOfG4Fluxes function g4_flux_family(gf::G4Flux; check::Bool = true)
  nv = passes_verticality_checks(gf)
  nb = breaks_non_abelian_gauge_group(gf)
  gfs = special_flux_family(model(gf), vert = nv, not_breaking = nb, check = check)
  return gfs
end


@doc raw"""
    integral_coefficients(gf::G4Flux)

Return the integral coefficients of a given $G_4$-flux.
If these coefficients are not known, an error is raised.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gfs = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: failed
  - Non-abelian gauge group: broken
  - Tadpole constraint: not analyzed

julia> g4 = random_flux_instance(gfs, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

julia> integral_coefficients(g4);
```
"""
function integral_coefficients(gf::G4Flux)
  @req has_attribute(gf, :int_combination) "Integral coefficients not known for the given G4-flux"
  return get_attribute(gf, :int_combination)::ZZMatrix
end


@doc raw"""
    rational_coefficients(gf::G4Flux)

Return the integral coefficients of a given $G_4$-flux.
If these coefficients are not known, an error is raised.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gfs = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Verticality checks: failed
  - Non-abelian gauge group: broken
  - Tadpole constraint: not analyzed

julia> g4 = random_flux_instance(gfs, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

julia> rational_coefficients(g4);
```
"""
function rational_coefficients(gf::G4Flux)
  @req has_attribute(gf, :rat_combination) "Rational coefficients not known for the given G4-flux"
  return get_attribute(gf, :rat_combination)::QQMatrix
end
