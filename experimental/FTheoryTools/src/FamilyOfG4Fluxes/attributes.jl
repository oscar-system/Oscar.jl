#####################################################
# 1 Basic attributes
#####################################################

@doc raw"""
    model(gf::FamilyOfG4Fluxes)

Return the F-theory model for which this family of $G_4$-flux candidates is defined.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> f_gs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole constraint: not analyzed

julia> model(f_gs) == qsm_model
true
```
"""
model(gf::FamilyOfG4Fluxes) = gf.model


@doc raw"""
    matrix_integral(gf::FamilyOfG4Fluxes)

Return the matrix whose columns specify those combinations of ambient space G4-flux
candidates, of which integral linear combinations are contained in this family
of G4-fluxes.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> f_gs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole constraint: not analyzed

julia> matrix_integral(f_gs) == mat_int
true
```
"""
matrix_integral(gf::FamilyOfG4Fluxes) = gf.mat_int


@doc raw"""
    matrix_rational(gf::FamilyOfG4Fluxes)

Return the matrix whose columns specify those combinations of ambient space G4-flux
candidates, of which rational linear combinations are contained in this family
of G4-fluxes.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> f_gs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole constraint: not analyzed

julia> matrix_rational(f_gs) == mat_rat
true
```
"""
matrix_rational(gf::FamilyOfG4Fluxes) = gf.mat_rat


#####################################################
# 2 Compute the D3-tadpole constraint
#####################################################

@doc raw"""
    d3_tadpole_constraint(fgs::FamilyOfG4Fluxes; check::Bool = true)

Return the d3-tapdole constraint of a family of G4-fluxes. Recall that for a given $G_4$-flux, this constraint
is $- \frac{1}{2} \cdot G_4^2 + \frac{1}{24} \cdot \chi(\widehat{Y}_4) \stackrel{!}{\geq} 0$.

Note that the family of fluxes is specified by linear combination of cohomology classes, some with rational
and some with integral coefficients. In terms of these coefficients the d3-tadpole constraint is the demand
that a quadratic polynomial in the coefficients evaluates to a non-negative number. This method returns said
polynomial in the cofficients. In order to evaluate the D3-tadpole for a particular $G_4$-flux, one has to
evaluate this polynomial for the numeric coefficient values that correspond to the given $G_4$-flux.

We use symbols $a_i$ to indicate integral coefficients and $r_i$ to indicate rational coefficients.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> fgs = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: broken
  - Tadpole constraint: not analyzed

julia> d3_tadpole_constraint(fgs);

julia> fgs
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: broken
  - Tadpole constraint: evaluated
```
"""
@attr QQMPolyRingElem function d3_tadpole_constraint(fgs::FamilyOfG4Fluxes; check::Bool = true)

  # Entry checks
  m = model(fgs)
  @req base_space(m) isa NormalToricVariety "Computation of D3-tadpole constraint only supported for toric base and ambient spaces"
  @req dim(ambient_space(m)) == 5 "Computation of D3-tadpole constraint only supported for 5-dimensional toric ambient spaces"
  if check
    @req is_complete(ambient_space(m)) "Computation of D3-tadpole constraint only supported for complete toric ambient spaces"
    @req is_simplicial(ambient_space(m)) "Computation of D3-tadpole constraint only supported for simplicial toric ambient space"
  end

  # Are intersection numbers known?
  inter_dict = get_attribute!(m, :inter_dict) do
    Dict{NTuple{4, Int64}, ZZRingElem}()
  end::Dict{NTuple{4, Int64}, ZZRingElem}
  s_inter_dict = get_attribute!(m, :s_inter_dict) do
    Dict{String, ZZRingElem}()
  end::Dict{String, ZZRingElem}
  
  # Compute data, that is used by the default/sophisticated intersection product
  if arxiv_doi(m) == "10.48550/arXiv.1511.03209"
    S = cox_ring(ambient_space(m))
    gS = gens(cox_ring(ambient_space(m)))
    linear_relations = matrix(QQ, rays(ambient_space(m)))
    scalings = [c.coeff for c in S.d]
    mnf = Oscar._minimal_nonfaces(ambient_space(m))
    sr_ideal_pos = [Vector{Int}(Polymake.row(mnf, i)) for i in 1:Polymake.nrows(mnf)]
    data = (
      S = S,
      gS = gS,
      linear_relations = linear_relations,
      scalings = scalings,
      sr_ideal_pos = sr_ideal_pos
    )
  else
    cy = polynomial(cohomology_class(toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))))
  end

  # Find the number of integral and rational parameters
  numb_int_parameters = ncols(matrix_integral(fgs))
  numb_rat_parameters = ncols(matrix_rational(fgs))

  # Create a polynomial ring with parameters ai for the integral_parameters and ri for the rational parameters
  amb_ring, my_gens = polynomial_ring(QQ, "a#" => 1:numb_int_parameters, "r#" => 1: numb_rat_parameters)

  # Extract ambient space basis of G4-flux candidates used to express flux family in
  basis = _ambient_space_models_of_g4_fluxes(m, check = check)
  basis_indices = get_attribute(m, :ambient_space_models_of_g4_fluxes_indices)::Vector{Tuple{Int64, Int64}}

  # Use MPolyBuildCtx to compute the tadpole constraint.
  C = MPolyBuildCtx(amb_ring)
  exp_vec = fill(0, numb_int_parameters + numb_rat_parameters)
  for k1 in 1:ngens(amb_ring)
    for k2 in k1:ngens(amb_ring)

      # Extract generator k1
      if numb_int_parameters >= k1
        gen1 = matrix_integral(fgs)[:, k1]
      else
        gen1 = matrix_rational(fgs)[:, k1 - numb_int_parameters]
      end

      # Extract generator k2
      if numb_int_parameters >= k2
        gen2 = matrix_integral(fgs)[:, k2]
      else
        gen2 = matrix_rational(fgs)[:, k2 - numb_int_parameters]
      end

      # Compute the intersection number of generator k1 and generator k2
      inter_number = ZZ(0)
      for l1 in 1:length(basis)
        for l2 in 1:length(basis)

          val = gen1[l1] * gen2[l2]
          is_zero(val) && continue

          my_tuple = Tuple(sort([basis_indices[l1]..., basis_indices[l2]...]))

          if arxiv_doi(m) == "10.48550/arXiv.1511.03209"
            change = sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data)
          else
            change = get!(inter_dict, my_tuple) do
              return QQ(integrate(cohomology_class(ambient_space(m), polynomial(basis[l1]) * polynomial(basis[l2]) * cy); check = check))
            end        
          end

          inter_number += val * change
          
        end
      end

      # Update the D3-tadpole constraint polynomial
      exp_vec[k1] += 1
      exp_vec[k2] += 1
      if k1 == k2
        push_term!(C, inter_number, exp_vec)
      else
        push_term!(C, 2 * inter_number, exp_vec)
      end
      exp_vec[k1] = 0
      exp_vec[k2] = 0

    end
  end
  tadpole_constraint_polynomial = finish(C)
  tadpole_constraint_polynomial = -1//2 * tadpole_constraint_polynomial + 1//24 * euler_characteristic(m, check = check)

  # Update the computed intersection numbers
  set_attribute!(m, :inter_dict, inter_dict)
  set_attribute!(m, :s_inter_dict, s_inter_dict)

  # Finally, return the result
  return tadpole_constraint_polynomial::QQMPolyRingElem

end
