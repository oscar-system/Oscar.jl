#####################################################
# 1 Basic attributes
#####################################################

@doc raw"""
    model(gf::G4Flux)

Return the F-theory model used to construct the ``G_4``-flux candidate.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> g4_candidate = qsm_flux(qsm_model)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: unbroken
  - Tadpole cancellation check: not computed

julia> model(g4_candidate) == qsm_model
true
```
"""
model(gf::G4Flux) = gf.model


@doc raw"""
    cohomology_class(gf::G4Flux)

Return the ambient space cohomology class which defines the ``G_4``-flux candidate.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> g4_class = cohomology_class(anticanonical_divisor_class(ambient_space(qsm_model)), quick = true)^2;

julia> g4f = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> cohomology_class(g4f);
```
"""
cohomology_class(gf::G4Flux) = gf.class



#####################################################
# 2 Compute the D3-tadpole constraint
#####################################################

@doc raw"""
    d3_tadpole_constraint(gf::G4Flux; check::Bool = true)

Return the d3-tapdole of a G4-flux, that is compute and return the quantity
``- \frac{1}{2} \cdot \int_{\widehat{Y_4}}{G_4 \wedge G_4} + \frac{1}{24} \cdot \chi(\widehat{Y}_4)``.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> g4 = qsm_flux(qsm_model)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: unbroken
  - Tadpole cancellation check: not computed

julia> d3_tadpole_constraint(g4, check = false)
12
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

Return the family of ``G_4``-fluxes sharing the following properties
with the given ``G_4``-flux: transversality and breaking of the
non-abelian gauge group.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> g4 = qsm_flux(qsm_model)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: unbroken
  - Tadpole cancellation check: not computed

julia> g4_flux_family(g4, check = false)
Family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: unbroken
```
"""
@attr FamilyOfG4Fluxes function g4_flux_family(gf::G4Flux; check::Bool = true)
  nb = breaks_non_abelian_gauge_group(gf)
  gfs = special_flux_family(model(gf), not_breaking = !nb, check = check)
  return gfs
end


function _g4_flux_coordinate_computer(gf::G4Flux)
  # Thus far, I cannot seem to create an example which would test the entires code below. So I hard-code something here...
  #=
  gff_mat_int = matrix(QQ, [[2],[4],[6]])
  gff_mat_rat = matrix(QQ, [[1],[2],[3]])
  mat_rhs = hcat(gff_mat_int, gff_mat_rat)
  lhs = matrix(QQ, [[2],[4],[6]])
  up = matrix(QQ, [[1//2],[1]])
  =#

  # 1. Identify the family that contains the flux
  gff = g4_flux_family(gf)
  gff_mat_int = matrix_integral(gff)
  gff_mat_rat = matrix_rational(gff)
  gff_offset = offset(gff)

  # 2. Find the coordiates of the flux
  flux_coordinates = get_attribute(gf, :flux_coordinates)

  # 3. We want to solve
  # 3. flux_coordinates = gff_mat_int * vector-with-Z-entries + gff_mat_rat * vector-with-Q-entries + offset
  # 3. Define suitable quantities lhs = mat_rhs * (vector-with-Z-entries, vector-with-Q-entries)
  lhs = flux_coordinates - gff_offset
  mat_rhs = hcat(gff_mat_int, gff_mat_rat)

  # 4. Find a particular solution over the rationals
  cs, up = can_solve_with_solution(mat_rhs, lhs, side = :right)
  @req cs "Inconsistency: Could not express flux in its family"

  # 5. Check degenerate case: All entries of cs are integers
  if all(is_integer, up[1:ncols(gff_mat_int)])
    int_coeffs = transpose(matrix(ZZ, [[up[k] for k in 1:ncols(gff_mat_int)]]))::ZZMatrix
    rat_coeffs = transpose(matrix(QQ, [[up[ncols(gff_mat_int) + k] for k in 1:ncols(gff_mat_rat)]]))::QQMatrix
    set_attribute!(gf, :int_combination, int_coeffs)
    set_attribute!(gf, :rat_combination, rat_coeffs)
    return get_attribute(gf, :int_combination)::ZZMatrix
  end

  # 6. Find the kernel of mat_rhs
  K = kernel(mat_rhs, side = :right)
  @req ncols(K) > 0 "Inconcistency encountered - trivial kernel"
  P = identity_matrix(QQ, ncols(gff_mat_int) + ncols(gff_mat_rat))
  for k in ncols(gff_mat_int)+1:ncols(gff_mat_int) + ncols(gff_mat_rat)
    P[k,k] = QQ(0)
  end

  # 7. Rephrase the problem: lhs = mat_rhs * (vector-with-Z-entries, vector-with-Q-entries)
  # 7. We found a vector up with rational entries above s.t. lhs = mat_rhs * up
  # 7. Of course, we can add any element in the kernel of mat_rhs, so it still holds:
  # 7. lhs = mat_rhs * (up + K * tilde(u))
  # 7. With the above projection map, we must thus seek a tilde(u) with
  # 7. P * (up + K * tilde(u)) in Z^l
  # 7. Define C = P * K, o = P * up, then we seek tilde(u) in Q^a with
  # 7. C * tilde(u) + o in Z^l
  C = P * K
  o_vec = P * up

  # 8. To solve this equation, we compute the SNF of C over ZZ: S = T * C * U
  # 8. In this expression T and U are invertible square matrices of ZZ and S a diagonal matrix with integer entries.
  denom = lcm(unique(vcat([[denominator(k) for k in C[l,:]] for l in 1:nrows(C)]...)))
  Ctilde = matrix(ZZ, denom * C)
  S, T, U = snf_with_transform(Ctilde)
  o_vec_prime = T * o_vec

  # 9. Now solve 1/denom * S * u_tilde_prime + o_vec_prime is an integer for all entries!
  u_tilde_prime = Vector{elem_type(QQ)}()
  r = rank(S)
  @req all(a -> S[a, a] != 0, 1:r) "Inconsistency encountered"
  @req all(a -> S[a, a] == 0, r+1:min(ncols(S), nrows(S))) "Inconsistency encountered"
  for k in 1:rank(S)
    push!(u_tilde_prime, (-1) * o_vec_prime[k] * denom // S[k, k])
  end
  u_tilde_prime = transpose(matrix(QQ, [u_tilde_prime]))

  # 10. Now compute the vector u
  u_vec = U * u_tilde_prime

  # 11. Now obtain a "good" particular solution
  good_up = up + K * u_vec
  
  # 12. Execute consistency checks
  @req mat_rhs * good_up == lhs "Inconsistency encountered"
  @req all(isinteger, good_up[:, 1:1:ncols(gff_mat_int)]) "Inconsistency encountered"

  # 13. Set the attributes and return the computed value
  int_coeffs = transpose(matrix(ZZ, [[good_up[k] for k in 1:ncols(gff_mat_int)]]))::ZZMatrix
  rat_coeffs = transpose(matrix(QQ, [[good_up[ncols(gff_mat_int) + k] for k in 1:ncols(gff_mat_rat)]]))::QQMatrix
  set_attribute!(gf, :int_combination, int_coeffs)
  set_attribute!(gf, :rat_combination, rat_coeffs)
  return get_attribute(gf, :int_combination)::ZZMatrix
end


@doc raw"""
    integral_coefficients(gf::G4Flux)

Return the integral coefficients of a ``G_4``-flux.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> gfs = special_flux_family(qsm_model, check = false, algorithm = "special")
Family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed

julia> g4 = random_flux_instance(gfs, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> integral_coefficients(g4);
```
"""
function integral_coefficients(gf::G4Flux)
  _g4_flux_coordinate_computer(gf)
  @req has_attribute(gf, :int_combination) "Integral coefficients could not be determined for the given G4-flux"
  return get_attribute(gf, :int_combination)::ZZMatrix
end


@doc raw"""
    rational_coefficients(gf::G4Flux)

Return the rational coefficients of a ``G_4``-flux.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> gfs = special_flux_family(qsm_model, check = false, algorithm = "special")
Family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed

julia> g4 = random_flux_instance(gfs, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> rational_coefficients(g4);
```
"""
function rational_coefficients(gf::G4Flux)
  _g4_flux_coordinate_computer(gf)
  @req has_attribute(gf, :rat_combination) "Rational coefficients could not be determined for the given G4-flux"
  return get_attribute(gf, :rat_combination)::QQMatrix
end


@doc raw"""
    offset(gf::G4Flux)

Return the offset of a ``G_4``-flux.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> gfs = special_flux_family(qsm_model, check = false)
Family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed

julia> g4 = random_flux_instance(gfs, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> offset(g4);
```
"""
function offset(gf::G4Flux)
  @req has_attribute(gf, :offset) "Offset not known for the given G4-flux"
  return get_attribute(gf, :offset)::Vector{QQFieldElem}
end
