##########################################
### (1) Fluxes within family of G4-fluxes
##########################################


@doc raw"""
    flux_instance(fgs::FamilyOfG4Fluxes, int_combination::ZZMatrix, rat_combination::QQMatrix; check::Bool = true)

Create an element of a family of G4-fluxes.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> shift = [zero(QQ) for k in 1:37];

julia> fgs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, shift, check = false)
Family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> int_combination = matrix(ZZ, [[3]])
[3]

julia> rat_combination = matrix(QQ, [[5//2]])
[5//2]

julia> flux_instance(fgs, int_combination, rat_combination, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> flux_instance(fgs, Int[], [], check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> flux_instance(fgs, [3], [], check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> flux_instance(fgs, [], [5//2], check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> flux_instance(fgs, [3], [5//2], check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed
```
"""
function flux_instance(fgs::FamilyOfG4Fluxes, int_combination::ZZMatrix, rat_combination::QQMatrix; check::Bool = true)
  @req ncols(int_combination) == 1 "int_combination is expected to be a single column vector"
  @req ncols(rat_combination) == 1 "rat_combination is expected to be a single column vector"
  @req nrows(int_combination) == ncols(matrix_integral(fgs)) "Number of specified integers must match the number of integral combinations in G4-flux family"
  @req nrows(rat_combination) == ncols(matrix_rational(fgs)) "Number of specified rationals must match the number of rational combinations in G4-flux family"
  m1 = matrix_integral(fgs) * int_combination
  m2 = matrix_rational(fgs) * rat_combination
  shift = offset(fgs)
  gens = chosen_g4_flux_gens(model(fgs), check = check)
  c1 = sum(m1[k,1] * gens[k] for k in 1:length(gens))
  c2 = sum(m2[k,1] * gens[k] for k in 1:length(gens))
  c3 = sum(shift[k] * gens[k] for k in 1:length(gens))
  flux = c1 + c2 + c3
  set_attribute!(flux, :int_combination, int_combination)
  set_attribute!(flux, :rat_combination, rat_combination)
  set_attribute!(flux, :offset, shift)
  set_attribute!(flux, :g4_flux_family, fgs)
  if has_attribute(fgs, :is_well_quantized)
    if is_well_quantized(fgs)
      set_attribute!(flux, :is_well_quantized, true)
    end
  end
  if has_attribute(fgs, :passes_transversality_checks)
    if passes_transversality_checks(fgs)
      set_attribute!(flux, :passes_transversality_checks, true)
    end
  end
  if has_attribute(fgs, :breaks_non_abelian_gauge_group)
    if !breaks_non_abelian_gauge_group(fgs)
      set_attribute!(flux, :breaks_non_abelian_gauge_group, false)
    end
  end
  return flux
end


function flux_instance(fgs::FamilyOfG4Fluxes, int_coeffs::Vector{Int}, rat_coeffs::Vector{Rational{Int}}; check::Bool = true)
  if length(int_coeffs) == 0 && length(rat_coeffs) == 0
    m = model(fgs)
    r = cohomology_ring(ambient_space(m), check = check)
    return g4_flux(m, cohomology_class(ambient_space(m), zero(r), quick = true), check = check)
  end
  @req all(x -> x isa Int, int_coeffs) "Provided integral coefficient is not an integer"
  @req all(x -> x isa Rational{Int64}, rat_coeffs) "Provided integral coefficient is not an integer"

  m_int = matrix(ZZ, [int_coeffs])
  if length(int_coeffs) == 0
    m_int = zero_matrix(ZZ, ncols(matrix_integral(fgs)), 1)
  end
  if length(int_coeffs) > 0
    @req length(int_coeffs) == ncols(matrix_integral(fgs)) "Number of specified integers must match the number of integral combinations in G4-flux family"
  end
  m_rat = matrix(QQ, [rat_coeffs])
  if length(rat_coeffs) == 0
    m_rat = zero_matrix(QQ, ncols(matrix_rational(fgs)), 1)
  end
  if length(rat_coeffs) > 0
    @req length(rat_coeffs) == ncols(matrix_rational(fgs)) "Number of specified rationals must match the number of rational combinations in G4-flux family"
  end
  return flux_instance(fgs, m_int, m_rat, check = check)
end

function flux_instance(fgs::FamilyOfG4Fluxes, int_coeffs::Vector{Any}, rat_coeffs::Vector{Rational{Int}}; check::Bool = true)
  return flux_instance(fgs, Vector{Int}(int_coeffs), rat_coeffs, check = check)
end

function flux_instance(fgs::FamilyOfG4Fluxes, int_coeffs::Vector{Int}, rat_coeffs::Vector{Any}; check::Bool = true)
  return flux_instance(fgs, int_coeffs, Vector{Rational{Int}}(rat_coeffs), check = check)
end

function flux_instance(fgs::FamilyOfG4Fluxes, int_coeffs::Vector{Any}, rat_coeffs::Vector{Any}; check::Bool = true)
  return flux_instance(fgs, Vector{Int}(int_coeffs), Vector{Rational{Int}}(rat_coeffs), check = check)
end


@doc raw"""
    random_flux_instance(fgs::FamilyOfG4Fluxes; check::Bool = true)

Create a random element of a family of G4-fluxes.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> shift = [zero(QQ) for k in 1:37];

julia> fgs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, shift, check = false)
Family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> random_flux_instance(fgs, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed
```
"""
function random_flux_instance(fgs::FamilyOfG4Fluxes; check::Bool = true)
  int_combination = zero_matrix(ZZ, ncols(matrix_integral(fgs)), 1)
  rat_combination = zero_matrix(QQ, ncols(matrix_rational(fgs)), 1)
  for i in 1:ncols(matrix_integral(fgs))
    int_combination[i] = rand(-100:100)
  end
  for i in 1:ncols(matrix_rational(fgs))
    numerator = rand(-100:100)
    denominator = rand(1:100)
    rat_combination[i] = numerator // denominator
  end
  return flux_instance(fgs, int_combination, rat_combination, check = check)
end



##########################################
### (2) Special random flux on model
##########################################

@doc raw"""
    random_flux(m::AbstractFTheoryModel; not_breaking::Bool = false, check::Bool = true)

Create a random ``G_4``-flux on a given F-theory model.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> rf = random_flux(qsm_model, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed
```
"""
function random_flux(m::AbstractFTheoryModel; not_breaking::Bool = false, check::Bool = true)
  family = special_flux_family(m, not_breaking = not_breaking, check = check)
  return random_flux_instance(family, check = check)
end
