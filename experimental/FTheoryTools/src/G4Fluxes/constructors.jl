################
# 1: Constructor
################

@doc raw"""
    g4_flux(model::AbstractFTheoryModel, class::CohomologyClass)

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
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

julia> g4f2 = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate
  - Elementary quantization checks: not executed
  - Verticality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed
```
"""
function g4_flux(m::AbstractFTheoryModel, g4_class::CohomologyClass; check::Bool = true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "G4-fluxes only supported for Weierstrass, global Tate and hypersurface models"
  @req base_space(m) isa NormalToricVariety "G4-flux currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "G4-flux currently supported only for toric ambient space"
  g4_candidate = G4Flux(m, g4_class)
  if check && !is_well_quantized(g4_candidate)
    error("Given G4-flux candidate found to violate quantization condition")
  end
  return g4_candidate
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

function Base.show(io::IO, g4::G4Flux)
  properties_string = ["G4-flux candidate"]

  # Check for elementary quantization checks
  if has_attribute(g4, :is_well_quantized)
    if is_well_quantized(g4)
      push!(properties_string, "  - Elementary quantization checks: satisfied")
    else
      push!(properties_string, "  - Elementary quantization checks: failed")
    end
  else
    push!(properties_string, "  - Elementary quantization checks: not executed")
  end

  # Check for verticality checks
  if has_attribute(g4, :is_vertical)
    if is_vertical(g4)
      push!(properties_string, "  - Verticality checks: satisfied")
    else
      push!(properties_string, "  - Verticality checks: failed")
    end
  else
    push!(properties_string, "  - Verticality checks: not executed")
  end

  # Check for non-abelian gauge group breaking
  if has_attribute(g4, :breaks_non_abelian_gauge_group)
    if breaks_non_abelian_gauge_group(g4)
      push!(properties_string, "  - Non-abelian gauge group: broken")
    else
      push!(properties_string, "  - Non-abelian gauge group: not broken")
    end
  else
    push!(properties_string, "  - Non-abelian gauge group: breaking pattern not analyzed")
  end

  # Check for tadpole cancellation checks
  if has_attribute(g4, :passes_tadpole_cancellation_check)
    if passes_tadpole_cancellation_check(g4)
      push!(properties_string, "  - Tadpole cancellation check: satisfied")
    else
      push!(properties_string, "  - Tadpole cancellation check: failed")
    end
  else
    push!(properties_string, "  - Tadpole cancellation check: not executed")
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
