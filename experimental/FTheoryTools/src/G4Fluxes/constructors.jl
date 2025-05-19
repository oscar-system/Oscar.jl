################
# 1: Constructor
################

@doc raw"""
    g4_flux(model::AbstractFTheoryModel, class::CohomologyClass; convert::Bool = false)

Construct a G4-flux candidate on an F-theory model. This functionality is
currently limited to
- Weierstrass models,
- global Tate models,
- hypersurface models.
Furthermore, our functionality requires a concrete geometry. That is,
the base space as well as the ambient space must be toric varieties.
In the toric ambient space $X_\Sigma$, the elliptically fibered space $Y$
that defines the F-theory model, is given by a hypersurface (cut out by
the Weierstrass, Tate or hypersurface polynomial, respectively).

In this setting, we assume that a $G_4$-flux candidate is represented by a
cohomology class $h$ in $H^{(2,2)} (X_\Sigma)$. The actual $G_4$-flux candidate
is then obtained by restricting $h$ to $Y$.

It is worth recalling that the $G_4$-flux candidate is subject to the quantization
condition $G_4 + \frac{1}{2} c_2(Y) \in H^{/2,2)}( Y_, \mathbb{Z})$ (see [Wit97](@cite)).
This condition is very hard to verify. However, it is relatively easy to gather
evidence for this condition to be satisfied/show that it is violated. To this end, let
$D_1$, $D_2$ be two toric divisors in $X_\Sigma$, then the topological intersection number
$\left[ h|_Y \right] \cdot \left[ P \right] \cdot \left[ D_1 \right] \cdot \left[ D_2 \right]$
must be an integer. Even this rather elementary check can be computationally expensive.
Users can therefore decide to skip this check upon construction by setting the parameter
`check` to the value `false`.

Another bottleneck can be the computation of the cohomology ring, which is necessary to
work with cohomology classes on the toric ambient space, which in turn define the G4-flux,
as explained above. The reason for this is, that by employing the theory explained in
[CLS11](@cite), we can only work out the cohomology ring of simpicial and complete (i.e. compact)
toric varieties. However, checking if a toric variety is complete (i.e. compact) can take
a long time. If the geometry in question is involved and you already know that the variety
is simplicial and complete, then we recommend to trigger the computation of the cohomology
ring with `check = false`. This will avoid this time consuming test.

Let us mention that you can also supply the option `convert = true`. This will turn the
provided cohomology class into the basis chosen internally.

An example is in order.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4))
Hypersurface model over a concrete base

julia> cohomology_ring(ambient_space(qsm_model), check = false);

julia> g4_class = cohomology_class(anticanonical_divisor_class(ambient_space(qsm_model)))^2;

julia> g4f = g4_flux(qsm_model, g4_class)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> g4f2 = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> cohomology_class(g4f2)
Cohomology class on a normal toric variety given by 2*x5*e2 + 4*x5*u + 6*x5*e4 + 4*x5*e1 + 2*x5*w + 2*x6*e2 + 4*x6*u + 6*x6*e4 + 4*x6*e1 + 2*x6*w + 2*x7*x12 + 2*x7*e2 + 4*x7*u + 6*x7*e4 + 4*x7*e1 + 2*x7*w + 8*x8*x23 + 4*x10*x23 - 4*x15*x26 + 2*x15*e2 + 4*x15*u + 6*x15*e4 + 4*x15*e1 + 2*x15*w - 2*x16^2 + 2*x16*e2 + 4*x16*u + 6*x16*e4 + 4*x16*e1 + 2*x16*w - 22*x17*x24 + 4*x17*e2 + 8*x17*u + 12*x17*e4 + 8*x17*e1 + 4*x17*w - 10*x18^2 + 7//2*x18*x25 + 4*x18*e2 + 8*x18*u + 12*x18*e4 + 8*x18*e1 + 4*x18*w + 6*x19*e2 + 12*x19*u + 18*x19*e4 + 12*x19*e1 + 6*x19*w + 8*x20^2 + 88//3*x20*x21 - 7*x20*x25 + 4*x20*e2 + 8*x20*u + 12*x20*e4 + 8*x20*e1 + 4*x20*w + 11//3*x21^2 - 77//3*x21*x24 + 4*x21*e2 + 8*x21*u + 12*x21*e4 + 8*x21*e1 + 4*x21*w + 31//3*x22^2 + 5//3*x22*x23 + 97//6*x22*x25 + 4*x23^2 - 8*x24^2 - 17//3*x24*x27 + 2*x24*e2 + 4*x24*u + 6*x24*e4 + 4*x24*e1 + 2*x24*w + 7//2*x25^2 + 2*x25*e2 + 4*x25*u + 6*x25*e4 + 4*x25*e1 + 2*x25*w - 2//3*x26*x27 + 5//3*x27^2 + 2*x27*e2 + 4*x27*u + 6*x27*e4 + 4*x27*e1 + 2*x27*w + x28^2 - 7//3*x29^2 + 5*e1*w

julia> g4f3 = g4_flux(qsm_model, g4_class, check = false, convert = true)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> cohomology_class(g4f3)
Cohomology class on a normal toric variety given by 2*x5*e2 + 4*x5*u + 6*x5*e4 + 4*x5*e1 + 2*x5*w + 6*x6*e4 + 2*x6*w + 2*x7*x12 + 2*x7*e2 + 4*x7*u + 6*x7*e4 + 4*x7*e1 + 2*x7*w - 4*x15*x26 + 2*x15*e2 + 4*x15*u + 6*x15*e4 + 4*x15*e1 + 2*x15*w - 2*x16^2 + 6*x16*e4 + 2*x16*w - 6*x17*x24 + 4*x17*e2 + 8*x17*u + 12*x17*e4 + 8*x17*e1 + 4*x17*w - 10*x18^2 - 1//2*x18*x25 + 4*x18*e2 + 8*x18*u + 12*x18*e4 + 8*x18*e1 + 4*x18*w - 16*x19*x20 + 6*x19*e2 + 12*x19*u + 18*x19*e4 + 12*x19*e1 + 6*x19*w + 8*x20^2 + 56//3*x20*x21 + x20*x25 + 4*x20*e2 + 8*x20*u + 12*x20*e4 + 8*x20*e1 + 4*x20*w + 19//3*x21^2 - 13//3*x21*x24 + 12*x21*e4 + 4*x21*w - 1//3*x22^2 + 1//3*x22*x23 - 7//6*x22*x25 + 7//3*x24*x27 + 6*x24*e4 + 2*x24*w - 1//2*x25^2 + 6*x25*e4 + 2*x25*w - 2//3*x26*x27 + 5//3*x27^2 + 6*x27*e4 + 2*x27*w + x28^2 + 1//3*x29^2 + 5*e1*w
```
"""
function g4_flux(m::AbstractFTheoryModel, g4_class::CohomologyClass; check::Bool = true, convert::Bool = false)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "G4-fluxes only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "G4-flux currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "G4-flux currently supported only for toric ambient space"

  # If conversion to internally chosen basis desired, then modify the cohomology class
  converted_class = g4_class
  if convert == true
    to_be_transformed_poly = lift(polynomial(g4_class))
    M = collect(exponents(to_be_transformed_poly))
    non_zero_exponents = Vector{Tuple{Int64, Int64}}()
    for my_row in M
      i1 = findfirst(x -> x != 0, my_row)
      my_row[i1] -= 1
      i2 = findfirst(x -> x != 0, my_row)
      push!(non_zero_exponents, (i1, i2))
    end
    coeffs = collect(coefficients(to_be_transformed_poly))
    @req length(coeffs) == length(non_zero_exponents) "Inconsistency encountered"

    converter_dict = converter_dict_h22_hypersurface(m, check = check)
    b_ring = base_ring(cohomology_ring(ambient_space(m), check = check))
    b_ring_gens = gens(b_ring)
    new_converter_dict = Dict{Tuple{Int64, Int64}, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}()
    for (key, value) in converter_dict
      new_converter_dict[key] = sum(k[1] * b_ring_gens[k[2][1]] * b_ring_gens[k[2][2]] for k in value)
    end
    converted_poly = zero(b_ring)
    for l in 1:length(coeffs)
      if haskey(new_converter_dict, non_zero_exponents[l])
        converted_poly += coeffs[l] * new_converter_dict[non_zero_exponents[l]]
      end
    end
    converted_poly = cohomology_ring(ambient_space(m), check = check)(converted_poly)
    converted_class = CohomologyClass(ambient_space(m), converted_poly)
    
  end

  # Build the G4-flux candidate
  g4_candidate = G4Flux(m, converted_class)

  # Execute quantization checks if desired and return the created object
  if check && !is_well_quantized(g4_candidate) && !passes_transversality_checks(g4_candidate)
    error("Given G4-flux candidate found to violate quantization and/or transversality condition")
  end
  return g4_candidate
end


@doc raw"""
    qsm_flux(qsm_model::AbstractFTheoryModel)

For an F-theory QSM [CHLLT19](@cite), this method creates the particularly chosen G4-flux in these models.

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
  gens_strings = symbols(cox_ring(ambient_space(qsm_model)))
  e1 = cohomology_class(divs[findfirst(x -> x == :e1, gens_strings)])
  e2 = cohomology_class(divs[findfirst(x -> x == :e2, gens_strings)])
  e4 = cohomology_class(divs[findfirst(x -> x == :e4, gens_strings)])
  u = cohomology_class(divs[findfirst(x -> x == :u, gens_strings)])
  v = cohomology_class(divs[findfirst(x -> x == :v, gens_strings)])
  pb_Kbar = cohomology_class(sum([divs[k] for k in 1:length(gens_strings)-7]))
  g4_class = (-3) // kbar3(qsm_model) * (5 * e1 * e4 + pb_Kbar * (-3 * e1 - 2 * e2 - 6 * e4 + pb_Kbar - 4 * u + v))
  my_flux = g4_flux(qsm_model, g4_class, convert = true, check = false)
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
  new_poly = R(polynomial(cohomology_class(g1)).f + polynomial(cohomology_class(g2)).f)
  new_cohomology_class = CohomologyClass(ambient_space(model(g1)), new_poly)
  return G4Flux(model(g1), new_cohomology_class)
end

Base.:-(g1::G4Flux, g2::G4Flux) = g1 + (-1) * g2

Base.:-(g::G4Flux) = (-1) * g

function Base.:*(c::T, g::G4Flux) where {T <: Union{IntegerUnion, QQFieldElem, Rational{Int64}}}
  R = parent(polynomial(cohomology_class(g)))
  new_poly = R(c * polynomial(cohomology_class(g)).f)
  new_cohomology_class = CohomologyClass(ambient_space(model(g)), new_poly)
  return G4Flux(model(g), new_cohomology_class)
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
