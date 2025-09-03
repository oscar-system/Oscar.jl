################
# 1: Constructor
################

@doc raw"""
    g4_flux(model::AbstractFTheoryModel, class::CohomologyClass)

Construct a candidate ``G_4``-flux for a resolved F-theory model from a given cohomology class
on the toric ambient space.

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
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> g4_class = cohomology_class(anticanonical_divisor_class(ambient_space(qsm_model)), completeness_check = false)^2;

julia> g4f = g4_flux(qsm_model, g4_class, completeness_check = false, consistency_check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed
```
"""
function g4_flux(m::AbstractFTheoryModel, g4_class::CohomologyClass; completeness_check::Bool = true, consistency_check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "G4-fluxes only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "G4-flux currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "G4-flux currently supported only for toric ambient space"

  # Step 1: Extract exponent pairs from input class
  poly = lift(polynomial(g4_class))
  coeffs = collect(coefficients(poly))
  exps = extract_exponent_pairs(collect(exponents(poly)))
  @req length(coeffs) == length(exps) "Mismatch between coefficients and exponent pairs"
  original_dict = Dict(zip(exps, coeffs))

  # Step 2: Transform using converter dictionary
  converter_dict = converter_dict_h22_hypersurface(m; completeness_check)
  b_ring = base_ring(cohomology_ring(ambient_space(m); completeness_check))
  b_ring_gens = gens(b_ring)
  converted_poly = zero(b_ring)
  for (exp_pair, coeff) in original_dict
    if haskey(converter_dict, exp_pair)
      terms = converter_dict[exp_pair]
      sum_term = sum(k[1] * b_ring_gens[k[2][1]] * b_ring_gens[k[2][2]] for k in terms)
      converted_poly += coeff * sum_term
    end
  end

  # Step 3: Extract exponent pairs again from converted poly
  new_coeffs = collect(coefficients(converted_poly))
  new_exps = extract_exponent_pairs(collect(exponents(converted_poly)))
  new_dict = Dict(zip(new_exps, new_coeffs))

  # Step 4: Read off flux coordinates from basis indices
  basis_indices = gens_of_h22_hypersurface_indices(m; completeness_check)
  flux_coords = [get(new_dict, b, 0) for b in basis_indices]

  # Step 5: Build cohomology class
  coh_ring = cohomology_ring(ambient_space(m); completeness_check)
  converted_poly = coh_ring(converted_poly)
  converted_class = cohomology_class(ambient_space(m), converted_poly; completeness_check)

  # Step 6: Build G4Flux and assign attributes
  g4 = G4Flux(m, converted_class)
  set_attribute!(g4, :offset, zeros(Int, length(chosen_g4_flux_gens(m))))
  set_attribute!(g4, :flux_coordinates, flux_coords)

  # Step 7: Final checks
  if consistency_check
    @req (is_well_quantized(g4) && passes_transversality_checks(g4)) "G4-flux candidate violates quantization and/or transversality condition"
  end

  return g4
end


# One helper function to avoid repeating exponent extraction logic
function extract_exponent_pairs(M::Vector{Vector{Int64}})
  pairs = Tuple{Int, Int}[]
  for k in 1:length(M)
    row = copy(M[k])
    i1 = findfirst(!=(0), row)
    row[i1] -= 1
    i2 = findfirst(!=(0), row)
    push!(pairs, (i1, i2))
  end
  return pairs
end


@doc raw"""
    qsm_flux(qsm_model::AbstractFTheoryModel)

Return the ``G_4``-flux associated with one of the Quadrillion F-theory
Standard models, as described in [CHLLT19](@cite CHLLT19).

This flux has been pre-validated to pass essential consistency checks.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> qsm_flux(qsm_model)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: unbroken
  - Tadpole cancellation check: not computed
```
"""
function qsm_flux(qsm_model::AbstractFTheoryModel)
  @req arxiv_doi(qsm_model) == "10.48550/arXiv.1903.00009" "Can only compute the QSM flux for a QSM model"
  divs = torusinvariant_prime_divisors(ambient_space(qsm_model))
  gens_strings = symbols(coordinate_ring(ambient_space(qsm_model)))
  e1 = cohomology_class(divs[findfirst(x -> x == :e1, gens_strings)])
  e2 = cohomology_class(divs[findfirst(x -> x == :e2, gens_strings)])
  e4 = cohomology_class(divs[findfirst(x -> x == :e4, gens_strings)])
  u = cohomology_class(divs[findfirst(x -> x == :u, gens_strings)])
  v = cohomology_class(divs[findfirst(x -> x == :v, gens_strings)])
  pb_Kbar = cohomology_class(sum([divs[k] for k in 1:length(gens_strings)-7]))
  g4_class = (-3) // kbar3(qsm_model) * (5 * e1 * e4 + pb_Kbar * (-3 * e1 - 2 * e2 - 6 * e4 + pb_Kbar - 4 * u + v))
  my_flux = g4_flux(qsm_model, g4_class, completeness_check = false, consistency_check = false)
  set_attribute!(my_flux, :is_well_quantized, true)
  set_attribute!(my_flux, :passes_transversality_checks, true)
  set_attribute!(my_flux, :breaks_non_abelian_gauge_group, false)
  return my_flux
end


################################################
# 2: Equality and hash
################################################

function Base.:(==)(gf1::G4Flux, gf2::G4Flux)
  # G4-fluxes can only be equal if they are defined for identically the same model
  model(gf1) !== model(gf2) && return false

  # Currently, can only decide equality for Weierstrass, global Tate and hypersurface models
  m = model(gf1)
  if (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) == false
    error("Can currently only decide equality of G4-fluxes for Weierstrass, global Tate and hypersurface models")
  end

  # Compute the cohomology class corresponding to the hypersurface equation
  cl = toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))
  cy = cohomology_class(cl)

  # Now can return the result
  return cy * cohomology_class(gf1) == cy * cohomology_class(gf2)

end

function Base.hash(gf::G4Flux, h::UInt)
  b = 0x92bd6ac4f87d834e % UInt
  h = hash(model(gf), h)
  h = hash(cohomology_class(gf), h)
  return xor(h, b)
end



################################################
# 3: Arithmetics
################################################

function Base.:+(g1::G4Flux, g2::G4Flux)
  @req model(g1) === model(g2) "The G4-fluxes must be defined on the same model"
  R = parent(polynomial(cohomology_class(g1)))
  new_poly = R(lift(polynomial(cohomology_class(g1))) + lift(polynomial(cohomology_class(g2))))
  new_cohomology_class = CohomologyClass(ambient_space(model(g1)), new_poly, true)
  return g4_flux(model(g1), new_cohomology_class, completeness_check = false, consistency_check = false)
end

Base.:-(g1::G4Flux, g2::G4Flux) = g1 + (-1) * g2

Base.:-(g::G4Flux) = (-1) * g

function Base.:*(c::T, g::G4Flux) where {T <: Union{IntegerUnion, QQFieldElem, Rational{Int64}}}
  R = parent(polynomial(cohomology_class(g)))
  new_poly = R(c * lift(polynomial(cohomology_class(g))))
  new_cohomology_class = CohomologyClass(ambient_space(model(g)), new_poly, true)
  return g4_flux(model(g), new_cohomology_class, completeness_check = false, consistency_check = false)
end



################################################
# 4: Display
################################################

# Detailed printing
function Base.show(io::IO, ::MIME"text/plain", g4::G4Flux)
  io = pretty(io)
  properties_string = ["G4-flux candidate"]

  # Check for elementary quantization checks
  if has_attribute(g4, :is_well_quantized) && get_attribute(g4, :is_well_quantized) !== nothing
    if is_well_quantized(g4)
      push!(properties_string, "  - Elementary quantization checks: satisfied")
    else
      push!(properties_string, "  - Elementary quantization checks: violated")
    end
  else
    push!(properties_string, "  - Elementary quantization checks: not executed")
  end

  # Check for transversality checks
  if has_attribute(g4, :passes_transversality_checks) && get_attribute(g4, :passes_transversality_checks) !== nothing
    if passes_transversality_checks(g4)
      push!(properties_string, "  - Transversality checks: satisfied")
    else
      push!(properties_string, "  - Transversality checks: violated")
    end
  else
    push!(properties_string, "  - Transversality checks: not executed")
  end

  # Check for non-abelian gauge group breaking
  if has_attribute(g4, :breaks_non_abelian_gauge_group) && get_attribute(g4, :breaks_non_abelian_gauge_group) !== nothing
    if breaks_non_abelian_gauge_group(g4)
      push!(properties_string, "  - Non-abelian gauge group: broken")
    else
      push!(properties_string, "  - Non-abelian gauge group: unbroken")
    end
  else
    push!(properties_string, "  - Non-abelian gauge group: breaking pattern not analyzed")
  end

  # Check for tadpole cancellation checks
  if has_attribute(g4, :passes_tadpole_cancellation_check) && get_attribute(g4, :passes_tadpole_cancellation_check) !== nothing
    if passes_tadpole_cancellation_check(g4)
      push!(properties_string, "  - Tadpole cancellation check: satisfied")
    else
      push!(properties_string, "  - Tadpole cancellation check: violated")
    end
  else
    push!(properties_string, "  - Tadpole cancellation check: not computed")
  end

  # Print each line separately, to avoid extra line break at the end
  for (i, line) in enumerate(properties_string)
    if i == length(properties_string)
      print(io, line) # Last line without extra newline
    else
      println(io, line) # Print all other lines with line break
    end
  end
end

# Terse and one line printing
function Base.show(io::IO, g4::G4Flux)
  print(io, "G4-flux candidate")
end
