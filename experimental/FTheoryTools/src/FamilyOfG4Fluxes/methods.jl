##########################################
### (1) Fluxes within family of G4-fluxes
##########################################


@doc raw"""
    flux_instance(fgs::FamilyOfG4Fluxes, int_combination::ZZMatrix, rat_combination::QQMatrix)

Create an element of a family of G4-fluxes.

!!! note "Completeness check"
  The implemented algorithm is guaranteed to work only for toric ambient spaces
  that are smooth and **complete**. Verifying completeness can be very time 
  consuming. To skip this check, pass the optional keyword argument 
  `completeness_check=false`.

!!! note "Consistency check"
  G4-fluxes must always be properly quantized and pass the transversality checks.
  Verifying those properties can be very time consuming. To skip these consistency checks,
  pass the optional keyword argument `consistency_check=false`.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> shift = [zero(QQ) for k in 1:37];

julia> fgs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, shift, completeness_check = false)
Family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> int_combination = matrix(ZZ, [[3]])
[3]

julia> rat_combination = matrix(QQ, [[5//2]])
[5//2]

julia> flux_instance(fgs, int_combination, rat_combination, completeness_check = false, consistency_check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> flux_instance(fgs, Int[], [], completeness_check = false, consistency_check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> flux_instance(fgs, [3], [], completeness_check = false, consistency_check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> flux_instance(fgs, [], [5//2], completeness_check = false, consistency_check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> flux_instance(fgs, [3], [5//2], completeness_check = false, consistency_check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed
```
"""
function flux_instance(fgs::FamilyOfG4Fluxes, int_combination::ZZMatrix, rat_combination::QQMatrix; completeness_check::Bool = true, consistency_check::Bool = true)
  @req ncols(int_combination) == 1 "int_combination is expected to be a single column vector"
  @req ncols(rat_combination) == 1 "rat_combination is expected to be a single column vector"
  @req nrows(int_combination) == ncols(matrix_integral(fgs)) "Number of specified integers must match the number of integral combinations in G4-flux family"
  @req nrows(rat_combination) == ncols(matrix_rational(fgs)) "Number of specified rationals must match the number of rational combinations in G4-flux family"
  m1 = matrix_integral(fgs) * int_combination
  m2 = matrix_rational(fgs) * rat_combination
  shift = offset(fgs)
  gens = chosen_g4_flux_gens(model(fgs); completeness_check)
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
  if consistency_check
    @req (is_well_quantized(flux) && passes_transversality_checks(flux)) "G4-flux candidate violates quantization and/or transversality condition"
  end
  return flux
end


function flux_instance(fgs::FamilyOfG4Fluxes, int_coeffs::Vector{Int}, rat_coeffs::Vector{Rational{Int}}; completeness_check::Bool = true, consistency_check::Bool = true)
  if length(int_coeffs) == 0 && length(rat_coeffs) == 0
    m = model(fgs)
    r = cohomology_ring(ambient_space(m); completeness_check)
    tcc = cohomology_class(ambient_space(m), zero(r); completeness_check)
    return g4_flux(m, tcc; completeness_check, consistency_check)
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
  return flux_instance(fgs, m_int, m_rat; completeness_check, consistency_check)
end

function flux_instance(fgs::FamilyOfG4Fluxes, int_coeffs::Vector{Any}, rat_coeffs::Vector{Rational{Int}}; completeness_check::Bool = true, consistency_check::Bool = true)
  return flux_instance(fgs, Vector{Int}(int_coeffs), rat_coeffs; completeness_check, consistency_check)
end

function flux_instance(fgs::FamilyOfG4Fluxes, int_coeffs::Vector{Int}, rat_coeffs::Vector{Any}; completeness_check::Bool = true, consistency_check::Bool = true)
  return flux_instance(fgs, int_coeffs, Vector{Rational{Int}}(rat_coeffs); completeness_check, consistency_check)
end

function flux_instance(fgs::FamilyOfG4Fluxes, int_coeffs::Vector{Any}, rat_coeffs::Vector{Any}; completeness_check::Bool = true, consistency_check::Bool = true)
  return flux_instance(fgs, Vector{Int}(int_coeffs), Vector{Rational{Int}}(rat_coeffs); completeness_check, consistency_check)
end


@doc raw"""
    random_flux_instance(fgs::FamilyOfG4Fluxes)

Create a random element of a family of G4-fluxes.

!!! note "Completeness check"
  The implemented algorithm is guaranteed to work only for toric ambient spaces
  that are smooth and **complete**. Verifying completeness can be very time 
  consuming. To skip this check, pass the optional keyword argument 
  `completeness_check=false`.

!!! note "Consistency check"
  G4-fluxes must always be properly quantized and pass the transversality checks.
  Verifying those properties can be very time consuming. To skip these consistency checks,
  pass the optional keyword argument `consistency_check=false`.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> mat_int = zero_matrix(QQ, 37, 1);

julia> mat_int[1,1] = 1;

julia> mat_rat = zero_matrix(QQ, 37, 1);

julia> mat_rat[2,1] = 1;

julia> shift = [zero(QQ) for k in 1:37];

julia> fgs = family_of_g4_fluxes(qsm_model, mat_int, mat_rat, shift, completeness_check = false)
Family of G4 fluxes:
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed

julia> random_flux_instance(fgs, completeness_check = false, consistency_check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed
```
"""
function random_flux_instance(fgs::FamilyOfG4Fluxes; completeness_check::Bool = true, consistency_check::Bool = true)
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
  return flux_instance(fgs, int_combination, rat_combination; completeness_check, consistency_check)
end



##########################################
### (2) Special random flux on model
##########################################

@doc raw"""
    random_flux(m::AbstractFTheoryModel)

Create a random ``G_4``-flux on a given F-theory model.

!!! note "Completeness check"
  The implemented algorithm is guaranteed to work only for toric ambient spaces
  that are smooth and **complete**. Verifying completeness can be very time 
  consuming. To skip this check, pass the optional keyword argument 
  `completeness_check=false`.

!!! note "Consistency check"
  G4-fluxes must always be properly quantized and pass the transversality checks.
  Verifying those properties can be very time consuming. To skip these consistency checks,
  pass the optional keyword argument `consistency_check=false`.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> rf = random_flux(qsm_model, completeness_check = false, consistency_check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed
```
"""
function random_flux(m::AbstractFTheoryModel; not_breaking::Bool = false, completeness_check::Bool = true, consistency_check::Bool = true)
  family = special_flux_family(m; not_breaking, completeness_check)
  return random_flux_instance(family; completeness_check, consistency_check)
end
