###############################
# Blowup helper functionality #
###############################

const _BLOWUP_INVARIANT_ATTRIBUTES = (
  :classes_of_model_sections,
  :classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes,
  :discriminant,
  :gauge_algebra,
  :global_gauge_group_quotient,
  :singular_loci,
)

function _copy_blowup_model_data!(
  target::AbstractFTheoryModel, source::AbstractFTheoryModel
)
  target.defining_classes = copy(defining_classes(source))
  _copy_model_metadata!(target, source)
  for name in _BLOWUP_INVARIANT_ATTRIBUTES
    if has_attribute(source, name)
      set_attribute!(target, name => get_attribute(source, name))
    end
  end
  return target
end

function _remap_toric_divisor_indices!(
  target::AbstractFTheoryModel,
  source::AbstractFTheoryModel,
  new_exceptional_names::Vector{Symbol}=Symbol[],
)
  source_names = symbols(coordinate_ring(ambient_space(source)))
  target_names = symbols(coordinate_ring(ambient_space(target)))
  names = [source_names[exceptional_divisor_indices(source)]; new_exceptional_names]
  indices = [findfirst(==(name), target_names)::Int for name in names]
  set_attribute!(target, :exceptional_divisor_indices => indices)
  if has_attribute(source, :zero_section_index)
    name = source_names[zero_section_index(source)]
    set_attribute!(target, :zero_section_index => findfirst(==(name), target_names)::Int)
  end
  return target
end

function _strict_transform_toric_sequence(
  blowups::AbstractVector{<:ToricBlowupMorphism}, f::MPolyDecRingElem
)
  @req !isempty(blowups) "The blowup sequence must not be empty"
  @req !is_zero(f) "The polynomial must not be zero"

  # Perform every strict transform on exponent vectors, but build the resulting
  # polynomial only once.
  exps = [collect(first(exponents(monomial))) for monomial in monomials(f)]
  R = parent(f)
  for blowup in blowups
    @req R === coordinate_ring(codomain(blowup)) "The toric blowups must form a composable sequence"
    @req index_of_exceptional_ray(blowup) == ngens(coordinate_ring(domain(blowup))) "The exceptional coordinate must be the final Cox ring generator"
    @req !has_torusfactor(codomain(blowup)) "Only implemented when there are no torus factors"
    coordinates = minimal_supercone_coordinates_of_exceptional_ray(blowup)
    for exponent in exps
      push!(exponent, ceil(Int, sum(coordinates .* exponent)))
    end
    common_factor = minimum(last, exps)
    for exponent in exps
      exponent[end] -= common_factor
    end
    R = coordinate_ring(domain(blowup))
  end

  builder = MPolyBuildCtx(R)
  for (coefficient, exponent) in zip(coefficients(f), exps)
    push_term!(builder, coefficient, exponent)
  end
  return finish(builder)
end

function _toric_blow_up_plan(
  m::AbstractFTheoryModel,
  centers::Vector{BlowupCenterType},
  exceptionals::Vector{String},
)
  isdefined(m, m isa GlobalTateModel ? :tate_polynomial : :weierstrass_polynomial) ||
    return ToricBlowupMorphism[], nothing

  X = ambient_space(m)
  has_torusfactor(X) && return ToricBlowupMorphism[], nothing
  sections = explicit_model_sections(m)
  blowups = ToricBlowupMorphism[]
  for (center, exceptional) in zip(centers, exceptionals)
    center = [string(get(sections, generator, generator)) for generator in center]
    R = coordinate_ring(X)
    I = ideal([eval_poly(generator, R) for generator in center])
    coordinates = _cox_ring_ideal_to_minimal_supercone_coordinates(X, I)
    isnothing(coordinates) && return blowups, ideal_sheaf(X, I)

    blowup = blow_up_along_minimal_supercone_coordinates(
      X, coordinates; coordinate_name=exceptional
    )
    @req index_of_exceptional_ray(blowup) == ngens(coordinate_ring(domain(blowup))) "Inconsistency encountered. Contact the authors"
    push!(blowups, blowup)
    X = domain(blowup)
  end
  return blowups, nothing
end

function _resolve_toric_blowup_sequence(
  m::Union{GlobalTateModel,WeierstrassModel},
  blowups::AbstractVector{<:ToricBlowupMorphism},
)
  polynomial = m isa GlobalTateModel ? tate_polynomial(m) : weierstrass_polynomial(m)
  polynomial = _strict_transform_toric_sequence(blowups, polynomial)
  T = m isa GlobalTateModel ? GlobalTateModel : WeierstrassModel
  model = T(
    explicit_model_sections(m),
    model_section_parametrization(m),
    polynomial,
    base_space(m),
    domain(last(blowups)),
  )
  _copy_blowup_model_data!(model, m)

  # Star subdivisions preserve the lattice, so their composition is the identity.
  blowdown = toric_morphism(
    ambient_space(model),
    identity_matrix(ZZ, ambient_dim(polyhedral_fan(ambient_space(m)))),
    ambient_space(m);
    check=false,
  )
  new_exceptional_names = [
    symbols(coordinate_ring(domain(blowup)))[index_of_exceptional_ray(blowup)] for
    blowup in blowups
  ]
  set_attribute!(model, :partially_resolved, true)
  set_attribute!(model, :blow_down_morphism, blowdown)
  _remap_toric_divisor_indices!(model, m, new_exceptional_names)
  return model
end

@doc raw"""
    blow_up(m::AbstractFTheoryModel, ideal_gens::Vector{String}; kwargs...)

Resolve an F-theory model by blowing up a locus in the ambient space.
Set `coordinate_name` to choose the name of the exceptional homogeneous coordinate; it defaults to `"e"`.

# Examples
```jldoctest
julia> using Random;

julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false, rng = Random.Xoshiro(1234))
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> blow_up(t, ["x", "y", "x1"]; coordinate_name = "e1")
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
Here is an example for a Weierstrass model.

# Examples
```jldoctest
julia> using Random;

julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> w = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, defining_classes = Dict("b" => b), completeness_check = false, rng = Random.Xoshiro(1234))
Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)

julia> blow_up(w, ["x", "y", "x1"]; coordinate_name = "e1")
Partially resolved Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)
```
"""
function blow_up(
  m::AbstractFTheoryModel,
  ideal_gens::Vector{String};
  coordinate_name::String="e",
)
  R = coordinate_ring(ambient_space(m))
  I = ideal([eval_poly(k, R) for k in ideal_gens])
  return blow_up(m, I; coordinate_name)
end

@doc raw"""
    blow_up(m::AbstractFTheoryModel, I::MPolyIdeal; kwargs...)

Resolve an F-theory model by blowing up a locus in the ambient space.
Set `coordinate_name` to choose the name of the exceptional homogeneous coordinate; it defaults to `"e"`.

# Examples
```jldoctest
julia> using Random;

julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false, rng = Random.Xoshiro(1234))
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> x1, x2, x3, x4, x, y, z = gens(coordinate_ring(ambient_space(t)))
7-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4
 x
 y
 z

julia> blow_up(t, ideal([x, y, x1]); coordinate_name = "e1")
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
"""
function blow_up(
  m::AbstractFTheoryModel,
  I::MPolyIdeal;
  coordinate_name::String="e",
)
  return _blow_up(m, ideal_sheaf(ambient_space(m), I); coordinate_name)
end

function _cox_ring_ideal_to_minimal_supercone_coordinates(
  X::NormalToricVarietyType, defining_ideal::MPolyIdeal
)
  ideal_generators = gens(defining_ideal)
  ring_generators = gens(base_ring(defining_ideal))
  all(generator -> generator in ring_generators, ideal_generators) || return nothing

  R = coordinate_ring(X)
  coordinates = zeros(QQFieldElem, n_rays(X))
  for i in 1:n_rays(X)
    R[i] in ideal_generators && (coordinates[i] = 1)
  end
  is_zero(coordinates) && return nothing
  is_minimal_supercone_coordinate_vector(polyhedral_fan(X), coordinates) || return nothing
  return coordinates
end

function _ideal_sheaf_to_minimal_supercone_coordinates(
  X::AbsCoveredScheme, I::AbsIdealSheaf
)
  X isa NormalToricVarietyType || return nothing
  I isa ToricIdealSheafFromCoxRingIdeal || return nothing
  return _cox_ring_ideal_to_minimal_supercone_coordinates(X, ideal_in_cox_ring(I))
end

@doc raw"""
    blow_up(m::AbstractFTheoryModel, I::AbsIdealSheaf; kwargs...)

Resolve an F-theory model by blowing up a locus in the ambient space.
For this method, the blowup center is encoded by an ideal sheaf.
Set `coordinate_name` to choose the name of the exceptional homogeneous coordinate; it defaults to `"e"`.

# Examples
```jldoctest
julia> using Random;

julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false, rng = Random.Xoshiro(1234))
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> x1, x2, x3, x4, x, y, z = gens(coordinate_ring(ambient_space(t)))
7-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4
 x
 y
 z

julia> blowup_center = ideal_sheaf(ambient_space(t), ideal([x, y, x1]))
Sheaf of ideals
  on q-factorial normal toric variety
with restrictions
   1: Ideal (x_5_1, x_4_1, x_1_1)
   2: Ideal (1)
   3: Ideal (x_5_3, x_4_3, x_1_3)
   4: Ideal (x_5_4, x_4_4, x_1_4)
   5: Ideal (1)
   6: Ideal (1)
   7: Ideal (1)
   8: Ideal (1)
   9: Ideal (1)
  10: Ideal (1)
  11: Ideal (1)
  12: Ideal (1)

julia> blow_up(t, blowup_center; coordinate_name = "e1")
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
"""
function blow_up(
  m::AbstractFTheoryModel,
  I::AbsIdealSheaf;
  coordinate_name::String="e",
)
  return _blow_up(m, I; coordinate_name)
end

function _blow_up(
  m::AbstractFTheoryModel,
  I::AbsIdealSheaf;
  coordinate_name::String,
  coords=missing,
)

  # Cannot (yet) blowup if this is not a Tate or Weierstrass model
  @req m isa Union{GlobalTateModel,WeierstrassModel} "Blowups are currently only supported for Tate and Weierstrass models"
  @req (base_space(m) isa FamilyOfSpaces) == false "Base space must be a concrete space for blowups to work"

  # Compute the new ambient_space
  ismissing(coords) &&
    (coords = _ideal_sheaf_to_minimal_supercone_coordinates(ambient_space(m), I))
  if !isnothing(coords)
    # Apply toric method
    bd = blow_up_along_minimal_supercone_coordinates(
      ambient_space(m), coords; coordinate_name
    )
  else
    # Reroute to scheme theory
    bd = blow_up(I)
  end
  new_ambient_space = domain(bd)

  # FIXME: Retaining the base is wrong for centers that modify the base. Supporting
  # those requires tracking the ambient-to-base projection and transforming sections.

  # Construct the new model
  if m isa GlobalTateModel
    if isdefined(m, :tate_polynomial) && new_ambient_space isa NormalToricVariety
      f = tate_polynomial(m)
      new_tate_polynomial = strict_transform(bd, f)
      model = GlobalTateModel(
        explicit_model_sections(m),
        model_section_parametrization(m),
        new_tate_polynomial,
        base_space(m),
        new_ambient_space,
      )
    else
      if bd isa ToricBlowupMorphism
        new_tate_ideal_sheaf = ideal_sheaf(
          domain(bd), strict_transform(bd, ideal_in_coordinate_ring(tate_ideal_sheaf(m)))
        )
      else
        new_tate_ideal_sheaf = strict_transform(bd, tate_ideal_sheaf(m))
      end
      model = GlobalTateModel(
        explicit_model_sections(m),
        model_section_parametrization(m),
        new_tate_ideal_sheaf,
        base_space(m),
        new_ambient_space,
      )
    end
  else
    if isdefined(m, :weierstrass_polynomial) && new_ambient_space isa NormalToricVariety
      f = weierstrass_polynomial(m)
      new_weierstrass_polynomial = strict_transform(bd, f)
      model = WeierstrassModel(
        explicit_model_sections(m),
        model_section_parametrization(m),
        new_weierstrass_polynomial,
        base_space(m),
        new_ambient_space,
      )
    else
      if bd isa ToricBlowupMorphism
        new_weierstrass_ideal_sheaf = ideal_sheaf(
          domain(bd),
          strict_transform(bd, ideal_in_coordinate_ring(weierstrass_ideal_sheaf(m))),
        )
      else
        new_weierstrass_ideal_sheaf = strict_transform(bd, weierstrass_ideal_sheaf(m))
      end
      model = WeierstrassModel(
        explicit_model_sections(m),
        model_section_parametrization(m),
        new_weierstrass_ideal_sheaf,
        base_space(m),
        new_ambient_space,
      )
    end
  end

  _copy_blowup_model_data!(model, m)
  set_attribute!(model, :partially_resolved, true)
  set_attribute!(model, :blow_down_morphism, bd)

  if ambient_space(model) isa NormalToricVariety
    index = index_of_exceptional_ray(bd)
    @req index == ngens(coordinate_ring(ambient_space(model))) "Inconsistency encountered. Contact the authors"
    names = symbols(coordinate_ring(ambient_space(model)))
    _remap_toric_divisor_indices!(model, m, [names[index]])
  end

  # Return the model
  return model
end

function _resolve_mixed_blowup_sequence(
  m::Union{GlobalTateModel,WeierstrassModel},
  centers::Vector{BlowupCenterType},
  exceptionals::Vector{String},
  toric_blowups::AbstractVector{<:ToricBlowupMorphism},
  non_toric_center::Union{Nothing,AbsIdealSheaf},
)
  n = length(centers)
  n_toric = length(toric_blowups)
  @req n_toric < n "The mixed blowup sequence must contain a non-toric blowup"

  # Retain the toric prefix constructed while planning instead of repeating it.
  resolved_model =
    isempty(toric_blowups) ? m : _resolve_toric_blowup_sequence(m, toric_blowups)
  blowdowns = AbsSimpleBlowupMorphism[toric_blowups...]
  for k in (n_toric + 1):n
    sections = explicit_model_sections(resolved_model)
    center = [string(get(sections, generator, generator)) for generator in centers[k]]

    if ambient_space(resolved_model) isa NormalToricVariety
      if !isnothing(non_toric_center) && k == n_toric + 1
        resolved_model = _blow_up(
          resolved_model,
          non_toric_center;
          coordinate_name=exceptionals[k],
          coords=nothing,
        )
      else
        X = ambient_space(resolved_model)
        R = coordinate_ring(X)
        I = ideal([eval_poly(generator, R) for generator in center])
        resolved_model = _blow_up(
          resolved_model, ideal_sheaf(X, I); coordinate_name=exceptionals[k]
        )
      end
    else
      # Transform the non-exceptional part from the initial ambient space.
      filtered_center = [c for c in center if !(c in exceptionals)]
      X = ambient_space(m)
      R = coordinate_ring(X)
      filtered_ideal = ideal_sheaf(
        X, ideal([eval_poly(generator, R) for generator in filtered_center])
      )
      for blowdown in blowdowns
        filtered_ideal = strict_transform(blowdown, filtered_ideal)
      end

      # Transform exceptional divisors only through their subsequent blowups.
      positions = [findfirst(==(c), exceptionals)::Int for c in center if c in exceptionals]
      exceptional_ideals = [
        ideal_sheaf(exceptional_divisor(blowdowns[position])) for
        position in positions
      ]
      for (i, position) in enumerate(positions)
        for l in (position + 1):length(blowdowns)
          exceptional_ideals[i] = strict_transform(blowdowns[l], exceptional_ideals[i])
        end
      end

      isempty(exceptional_ideals) || (filtered_ideal += sum(exceptional_ideals))
      resolved_model = blow_up(
        resolved_model, filtered_ideal; coordinate_name=exceptionals[k]
      )
    end

    push!(blowdowns, get_attribute(resolved_model, :blow_down_morphism))
  end
  return resolved_model
end

@doc raw"""
    resolve(m::AbstractFTheoryModel, index::Int; use_resolved_model_artifact::Bool=true)

Resolve a model with the index-th resolution that is known.

For the large model associated with arXiv:1511.03209, the resolved model is loaded from an artifact by default. Set `use_resolved_model_artifact=false` to execute the complete resolution algorithm instead. The direct computation can take days or longer.

Toric blowups are batched whenever possible, so their ambient spaces are
constructed first and the hypersurface strict transform is computed only once.

# Examples
```jldoctest
julia> using Random;

julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false, rng = Random.Xoshiro(1234))
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> t2 = resolve(t, 1)
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> coordinate_ring(ambient_space(t2))
Multivariate polynomial ring in 12 variables over QQ graded by
  x1 -> [1 0 0 0 0 0 0]
  x2 -> [0 1 0 0 0 0 0]
  x3 -> [0 1 0 0 0 0 0]
  x4 -> [0 1 0 0 0 0 0]
  x -> [0 0 1 0 0 0 0]
  y -> [0 0 0 1 0 0 0]
  z -> [0 0 0 0 1 0 0]
  e1 -> [0 0 0 0 0 1 0]
  e4 -> [0 0 0 0 0 0 1]
  e2 -> [-1 -3 -1 1 -1 -1 0]
  e3 -> [0 4 1 -1 1 0 -1]
  s -> [2 6 -1 0 2 1 1]

julia> w2 = 2 * torusinvariant_prime_divisors(B3)[1]
Torus-invariant, non-prime divisor on a normal toric variety

julia> t3 = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w2), completeness_check = false, rng = Random.Xoshiro(1234))
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> t4 = resolve(t3, 1)
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
"""
function resolve(
  m::AbstractFTheoryModel,
  resolution_index::Int;
  use_resolved_model_artifact::Bool=true,
)

  # To be extended to hypersurface models...
  @req m isa Union{GlobalTateModel,WeierstrassModel} "Resolve currently supported only for Weierstrass and Tate models"
  @req has_attribute(m, :resolutions) "No resolutions known for this model"
  @req resolution_index > 0 "The resolution must be specified by a positive integer"
  @req resolution_index <= length(resolutions(m)) "The resolution must be specified by an integer that is not larger than the number of known resolutions"
  @req (base_space(m) isa NormalToricVariety) "Currently, resolve is only supported for models over concrete toric bases"
  @req (ambient_space(m) isa NormalToricVariety) "Currently, resolve is only supported for singular models defined in a toric space"

  is_resolution_1511_03209 =
    has_attribute(m, :arxiv_id) && arxiv_id(m) == "1511.03209" && resolution_index == 1
  if is_resolution_1511_03209 && use_resolved_model_artifact
    resolved_model = load(
      artifact"FTM-1511-03209/1511-03209-resolved.mrdi"
    )::GlobalTateModel
    delete!.(Ref(resolved_model.__attrs), (:resolutions, :weighted_resolutions))
    return resolved_model
  end

  # Gather information for resolution
  centers, exceptionals = resolutions(m)[resolution_index]
  n = length(centers)
  @req n > 0 "A resolution must contain at least one blowup"

  toric_blowups, non_toric_center = _toric_blow_up_plan(m, centers, exceptionals)
  n_toric = length(toric_blowups)
  resolved_model =
    if n_toric == n
      _resolve_toric_blowup_sequence(m, toric_blowups)
    else
      _resolve_mixed_blowup_sequence(
        m, centers, exceptionals, toric_blowups, non_toric_center
      )
    end

  if is_resolution_1511_03209
    toric_model = resolved_model
    resolved_model = _resolve_ambient_space_1511_03209(toric_model)
    _copy_blowup_model_data!(resolved_model, toric_model)
    _remap_toric_divisor_indices!(resolved_model, toric_model)
  end
  return resolved_model
end
