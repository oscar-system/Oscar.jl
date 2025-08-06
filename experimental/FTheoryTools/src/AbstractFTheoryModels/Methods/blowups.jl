##################################################################
# 1: Blowups (old helper function, to be used for family of bases)
##################################################################

function _blowup_global(id::MPolyIdeal{QQMPolyRingElem}, center::MPolyIdeal{QQMPolyRingElem}, irr::MPolyIdeal{QQMPolyRingElem}, sri::MPolyIdeal{QQMPolyRingElem}, lin::MPolyIdeal{<:MPolyRingElem}; index::Integer = 1)
  # @warn "The function _blowup_global is experimental; absence of bugs and proper results are not guaranteed"
  
  R = base_ring(id)
  center_size = ngens(center)
  
  # Various sanity checks
  @req (!is_zero(center)) "The blowup center must be non-empty"
  
  # @req is_subset(id, center) "The ideal of the blowup center must contain the ideal to be blown up"
  @req base_ring(irr) == R "The given irrelevant ideal must share the base ring of the ideal to be blown up"
  @req base_ring(sri) == R "The given Stanley–Reisner ideal must share the base ring of the ideal to be blown up"
  @req ngens(base_ring(lin)) == ngens(R) "The base ring of ideal of linear relations must have the same number of generators as the base ring of the ideal to be blown up"
  
  # Make sure the ideal of linear relations has the same base ring as the others
  lin = ideal(map(hom(base_ring(lin), R, collect(1:ngens(R))), gens(lin)))
  
  # Create new base ring for the blown up ideal and a map between the rings
  S, S_gens = polynomial_ring(QQ, [Symbol("e_", index); [Symbol("b_", index, "_", i) for i in 1:center_size]; symbols(R)], cached = false)
  (_e, new_coords...) = S_gens[1:center_size + 1]
  ring_map = hom(R, S, S_gens[center_size + 2:end])
  
  # Compute the total transform
  center_gens_S = map(ring_map, gens(center))
  total_transform = ideal(map(ring_map, gens(id))) + ideal([new_coords[i] * _e - center_gens_S[i] for i in 1:center_size])
  
  # Compute the exceptional locus and strict transform, checking for crepancy
  # Could alternatively replace _e with center_gens_S in the exceptional locus here, then take the
  # primary decomposition and remove parts whose saturation by the irrelevant ideal is the whole ring
  exceptional_ideal = total_transform + ideal([_e])
  strict_transform, exceptional_factor = saturation_with_index(total_transform, exceptional_ideal)
  crepant = (exceptional_factor == center_size - 1)
  
  # Compute the new irrelevant ideal, SRI, and ideal of linear relations
  # These may need to be changed after reintroducing e
  new_irr = ideal(map(ring_map, gens(irr))) * ideal(new_coords)
  new_sri = ideal(map(ring_map, gens(sri))) + ideal([prod(new_coords)])
  new_lin = ideal(map(ring_map, gens(lin))) + ideal([g - new_coords[end] for g in new_coords[1:end - 1]])
  
  return total_transform, strict_transform, exceptional_ideal, crepant, new_irr, new_sri, new_lin, S, S_gens, ring_map
end
_blowup_global(id::T, center::T, irr::T, sri::T, lin::T; index::Integer = 1) where {T<:MPolyIdeal{<:MPolyRingElem}} = _blowup_global(ideal(map(g -> lift(g), gens(id))), ideal(map(g -> lift(g), gens(center))), ideal(map(g -> lift(g), gens(irr))), ideal(map(g -> lift(g), gens(sri))), lin, index = index)


function _blowup_global_sequence(id::MPolyIdeal{QQMPolyRingElem}, centers::Vector{<:Vector{<:Integer}}, irr::MPolyIdeal{QQMPolyRingElem}, sri::MPolyIdeal{QQMPolyRingElem}, lin::MPolyIdeal{<:MPolyRingElem}; index::Integer = 1)
  # @warn "The function _blowup_global_sequence is experimental; absence of bugs and proper results are not guaranteed"
  
  (cur_strict_transform, cur_irr, cur_sri, cur_lin, cur_S, cur_S_gens, cur_index) = (id, irr, sri, lin, base_ring(id), gens(base_ring((id))), index)
  crepant = true
  ring_map = hom(cur_S, cur_S, cur_S_gens) # Identity map
  
  exceptionals = MPolyIdeal{<:MPolyRingElem}[]
  for center in centers
    @req all(ind -> 1 <= ind <= length(cur_S_gens), center) "The given indices for the center generators are out of bounds"
    
    (_, cur_strict_transform, cur_ex, cur_crep, cur_irr, cur_sri, cur_lin, cur_S, cur_S_gens, cur_ring_map) = _blowup_global(cur_strict_transform, ideal(map(ind -> cur_S_gens[ind], center)), cur_irr, cur_sri, cur_lin, index = cur_index)
    
    map!(cur_ring_map, exceptionals, exceptionals)
    push!(exceptionals, cur_ex)
    
    crepant = crepant && cur_crep
    
    ring_map = compose(ring_map, cur_ring_map)
    
    cur_index += 1
  end
  
  return cur_strict_transform, exceptionals, crepant, cur_irr, cur_sri, cur_lin, cur_S, cur_S_gens, ring_map
end
_blowup_global_sequence(id::T, centers::Vector{<:Vector{<:Integer}}, irr::T, sri::T, lin::T; index::Integer = 1) where {T<:MPolyIdeal{<:MPolyRingElem}} = _blowup_global_sequence(ideal(map(g -> lift(g), gens(id))), centers, ideal(map(g -> lift(g), gens(irr))), ideal(map(g -> lift(g), gens(sri))), lin, index = index)



###########################################################################
# 2: Convenience functions for blowups
# 2: FOR INTERNAL USE ONLY (as of Feb 1, 2025 and PR 4523)
# 2: They are not in use (as of Feb 1, 2025 and PR 4523)
# 2: Gauge in the future if they are truly needed!
###########################################################################

@doc raw"""
    _martins_desired_blowup(m::NormalToricVariety, I::ToricIdealSheafFromCoxRingIdeal; coordinate_name::String = "e")

Blow up the toric variety along a toric ideal sheaf.

!!! warning
    This function is type unstable. The type of the domain of the output `f` is always a subtype of `AbsCoveredScheme` (meaning that `domain(f) isa AbsCoveredScheme` is always true). 
    Sometimes, the type of the domain will be a toric variety (meaning that `domain(f) isa NormalToricVariety` is true) if the algorithm can successfully detect this.
    In the future, the detection algorithm may be improved so that this is successful more often.

!!! warning
    This is an internal method. It is NOT exported.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> x1, x2, x3, x4 = gens(cox_ring(P3))
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4

julia> II = ideal_sheaf(P3, ideal([x1*x2]))
Sheaf of ideals
  on normal toric variety
with restrictions
  1: Ideal (x_1_1*x_2_1)
  2: Ideal (x_2_2)
  3: Ideal (x_1_3)
  4: Ideal (x_1_4*x_2_4)

julia> f = Oscar._martins_desired_blowup(P3, II);
```
"""
function _martins_desired_blowup(v::NormalToricVarietyType, I::ToricIdealSheafFromCoxRingIdeal; coordinate_name::Union{String, Nothing} = nothing)
  coords = _ideal_sheaf_to_minimal_supercone_coordinates(v, I)
  if !isnothing(coords)
    return blow_up_along_minimal_supercone_coordinates(v, coords; coordinate_name=coordinate_name) # Apply toric method
  else
    return blow_up(I) # Reroute to scheme theory
  end
end


@doc raw"""
    _martins_desired_blowup(v::NormalToricVariety, I::MPolyIdeal; coordinate_name::String = "e")

Blow up the toric variety by subdividing the cone in the list
of *all* cones of the fan of `v` which corresponds to the
provided ideal `I`. Note that this cone need not be maximal.

By default, we pick "e" as the name of the homogeneous coordinate for
the exceptional prime divisor. As third optional argument one can supply
a custom variable name.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> (x1,x2,x3,x4) = gens(cox_ring(P3))
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4

julia> I = ideal([x2,x3])
Ideal generated by
  x2
  x3

julia> bP3 = domain(Oscar._martins_desired_blowup(P3, I))
Normal toric variety

julia> cox_ring(bP3)
Multivariate polynomial ring in 5 variables over QQ graded by
  x1 -> [1 0]
  x2 -> [0 1]
  x3 -> [0 1]
  x4 -> [1 0]
  e -> [1 -1]

julia> I2 = ideal([x2 * x3])
Ideal generated by
  x2*x3

julia> b2P3 = Oscar._martins_desired_blowup(P3, I2);

julia> codomain(b2P3) == P3
true
```
"""
function _martins_desired_blowup(v::NormalToricVarietyType, I::MPolyIdeal; coordinate_name::Union{String, Nothing} = nothing)
  return _martins_desired_blowup(v, ideal_sheaf(v, I))
end


##################################################################
# 3: Currently used blowup functionality
##################################################################

@doc raw"""
    blow_up(m::AbstractFTheoryModel, ideal_gens::Vector{String}; coordinate_name::String = "e")

Resolve an F-theory model by blowing up a locus in the ambient space.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> blow_up(t, ["x", "y", "x1"]; coordinate_name = "e1")
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
Here is an example for a Weierstrass model.

# Examples
```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> w = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, defining_classes = Dict("b" => b), completeness_check = false)
Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)

julia> blow_up(w, ["x", "y", "x1"]; coordinate_name = "e1")
Partially resolved Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)
```
"""
function blow_up(m::AbstractFTheoryModel, ideal_gens::Vector{String}; coordinate_name::String = "e", nr_of_current_blow_up::Int = 1, nr_blowups_in_sequence::Int = 1)
  R = cox_ring(ambient_space(m))
  I = ideal([eval_poly(k, R) for k in ideal_gens])
  return blow_up(m, I; coordinate_name = coordinate_name, nr_of_current_blow_up = nr_of_current_blow_up, nr_blowups_in_sequence = nr_blowups_in_sequence)
end


@doc raw"""
    blow_up(m::AbstractFTheoryModel, I::MPolyIdeal; coordinate_name::String = "e")

Resolve an F-theory model by blowing up a locus in the ambient space.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> x1, x2, x3, x4, x, y, z = gens(cox_ring(ambient_space(t)))
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
function blow_up(m::AbstractFTheoryModel, I::MPolyIdeal; coordinate_name::String = "e", nr_of_current_blow_up::Int = 1, nr_blowups_in_sequence::Int = 1)
  return blow_up(m, ideal_sheaf(ambient_space(m), I); coordinate_name = coordinate_name, nr_of_current_blow_up = nr_of_current_blow_up, nr_blowups_in_sequence = nr_blowups_in_sequence)
end

function _ideal_sheaf_to_minimal_supercone_coordinates(X::AbsCoveredScheme, I::AbsIdealSheaf; coordinate_name::String = "e")
  # Return this when cannot convert ideal to minimal supercone coordinates
  not_possible = nothing

  # X needs to be a smooth toric variety
  X isa NormalToricVarietyType || return not_possible

  I isa ToricIdealSheafFromCoxRingIdeal || return not_possible
  defining_ideal = ideal_in_cox_ring(I)
  all(in(gens(base_ring(defining_ideal))), gens(defining_ideal)) || return not_possible
  R = cox_ring(X)
  coords = zeros(QQ, n_rays(X))
  for i in 1:n_rays(X)
    R[i] in gens(defining_ideal) && (coords[i] = 1)
  end
  coords == zeros(QQ, n_rays(X)) && return not_possible
  is_minimal_supercone_coordinate_vector(polyhedral_fan(X), coords) || return not_possible
  return coords
end

@doc raw"""
    blow_up(m::AbstractFTheoryModel, I::AbsIdealSheaf; coordinate_name::String = "e")

Resolve an F-theory model by blowing up a locus in the ambient space.
For this method, the blowup center is encoded by an ideal sheaf.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> x1, x2, x3, x4, x, y, z = gens(cox_ring(ambient_space(t)))
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
  on normal, simplicial toric variety
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
function blow_up(m::AbstractFTheoryModel, I::AbsIdealSheaf; coordinate_name::String = "e", nr_of_current_blow_up::Int = 1, nr_blowups_in_sequence::Int = 1)
  
  # Cannot (yet) blowup if this is not a Tate or Weierstrass model
  entry_test = (m isa GlobalTateModel) || (m isa WeierstrassModel)
  @req entry_test "Blowups are currently only supported for Tate and Weierstrass models"
  @req (base_space(m) isa FamilyOfSpaces) == false "Base space must be a concrete space for blowups to work"

  # Compute the new ambient_space
  coords = _ideal_sheaf_to_minimal_supercone_coordinates(ambient_space(m), I)
  if !isnothing(coords)
    # Apply toric method
    bd = blow_up_along_minimal_supercone_coordinates(
      ambient_space(m), coords; coordinate_name=coordinate_name
    )
  else
    # Reroute to scheme theory
    bd = blow_up(I)
  end
  new_ambient_space = domain(bd)

  # Compute the new base
  # FIXME: THIS WILL IN GENERAL BE WRONG! IN PRINCIPLE, THE ABOVE ALLOWS TO BLOW UP THE BASE AND THE BASE ONLY.
  # FIXME: We should save the projection \pi from the ambient space to the base space.
  # FIXME: This is also ties in with the model sections to be saved, see below. Should the base change, so do these sections...
  new_base = base_space(m)

  # Construct the new model
  if m isa GlobalTateModel
    if isdefined(m, :tate_polynomial) && new_ambient_space isa NormalToricVariety
      f = tate_polynomial(m)
      new_tate_polynomial = strict_transform(bd, f)
      model = GlobalTateModel(explicit_model_sections(m), model_section_parametrization(m), new_tate_polynomial, base_space(m), new_ambient_space)
    else
      if bd isa ToricBlowupMorphism
        new_tate_ideal_sheaf = ideal_sheaf(domain(bd), strict_transform(bd, ideal_in_cox_ring(tate_ideal_sheaf(m))))
      else
        new_tate_ideal_sheaf = strict_transform(bd, tate_ideal_sheaf(m))
      end
      model = GlobalTateModel(explicit_model_sections(m), model_section_parametrization(m), new_tate_ideal_sheaf, base_space(m), new_ambient_space)
    end
  else
    if isdefined(m, :weierstrass_polynomial) && new_ambient_space isa NormalToricVariety
      f = weierstrass_polynomial(m)
      new_weierstrass_polynomial = strict_transform(bd, f)
      model = WeierstrassModel(explicit_model_sections(m), model_section_parametrization(m), new_weierstrass_polynomial, base_space(m), new_ambient_space)
    else
      if bd isa ToricBlowupMorphism
        new_weierstrass_ideal_sheaf = ideal_sheaf(domain(bd), strict_transform(bd, ideal_in_cox_ring(weierstrass_ideal_sheaf(m))))
      else
        new_weierstrass_ideal_sheaf = strict_transform(bd, weierstrass_ideal_sheaf(m))
      end
      model = WeierstrassModel(explicit_model_sections(m), model_section_parametrization(m), new_weierstrass_ideal_sheaf, base_space(m), new_ambient_space)
    end
  end

  # Copy/overwrite/set attributes
  model_attributes = m.__attrs
  for (key, value) in model_attributes
    set_attribute!(model, key, value)
  end

  set_attribute!(model, :partially_resolved, true)
  set_attribute!(model, :blow_down_morphism, bd)

  if ambient_space(model) isa NormalToricVariety
    index = index_of_exceptional_ray(bd)
    @req index == ngens(cox_ring(ambient_space(model))) "Inconsistency encountered. Contact the authors"
    indices = exceptional_divisor_indices(model)
    push!(indices, index)
    set_attribute!(model, :exceptional_divisor_indices, indices)

    #Update slow attributes only at the end of a blow up sequence, if possible
    if nr_of_current_blow_up == nr_blowups_in_sequence
      # Update exceptional classes and their indices
      divs = torusinvariant_prime_divisors(ambient_space(model))

      indets = [lift(g) for g in gens(cohomology_ring(ambient_space(model), check = false))]
      coeff_ring = coefficient_ring(ambient_space(model))
      new_e_classes = Vector{CohomologyClass}()
      for i in indices
        poly = sum(coeff_ring(coefficients(divs[i])[k]) * indets[k] for k in 1:length(indets))
        push!(new_e_classes, CohomologyClass(ambient_space(model), cohomology_ring(ambient_space(model), check = false)(poly), true))
      end

      set_attribute!(model, :exceptional_classes, new_e_classes)
    end
  end

  # Return the model
  return model
end


@doc raw"""
    resolve(m::AbstractFTheoryModel, index::Int)

Resolve a model with the index-th resolution that is known.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> t2 = resolve(t, 1)
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> cox_ring(ambient_space(t2))
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

julia> t3 = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w2), completeness_check = false)
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> t4 = resolve(t3, 1)
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
"""
function resolve(m::AbstractFTheoryModel, resolution_index::Int)

  # For model 1511.03209 and resolution_index = 1, a particular resolution is available from an artifact
  if has_attribute(m, :arxiv_id)
    if resolution_index == 1 && arxiv_id(m) == "1511.03209"
      model_data_path = artifact"FTM-1511-03209/1511-03209-resolved.mrdi"
      model = load(model_data_path)
      # Modify attributes, see PR #5031 for details
      set_attribute!(model, :gens_of_h22_hypersurface, get_attribute(model, :basis_of_h22_hypersurface))
      set_attribute!(model, :gens_of_h22_hypersurface_indices, get_attribute(model, :basis_of_h22_hypersurface_indices))
      delete!(model.__attrs, :basis_of_h22_hypersurface)
      delete!(model.__attrs, :basis_of_h22_hypersurface_indices)
      return model
    end
  end

  # To be extended to hypersurface models...
  entry_test = (m isa GlobalTateModel) || (m isa WeierstrassModel)
  @req entry_test "Resolve currently supported only for Weierstrass and Tate models"
  @req (base_space(m) isa NormalToricVariety) "Currently, resolve is only supported for models over concrete toric bases"
  @req (ambient_space(m) isa NormalToricVariety) "Currently, resolve is only supported for singular models defined in a toric space"
  @req has_attribute(m, :resolutions) "No resolutions known for this model"
  @req resolution_index > 0 "The resolution must be specified by a non-negative integer"
  @req resolution_index <= length(resolutions(m)) "The resolution must be specified by an integer that is not larger than the number of known resolutions"
  
  # Gather information for resolution
  centers, exceptionals = resolutions(m)[resolution_index]
  nr_blowups = length(centers)

  # Resolve the model
  resolved_model = m
  blow_up_chain = []
  for k in 1:nr_blowups

    # Replace parameters in the blow_up_center with explicit_model_sections
    blow_up_center = centers[k]
    for l in 1:length(blow_up_center)
      model_sections = explicit_model_sections(resolved_model)
      if haskey(model_sections, blow_up_center[l])
        new_locus = string(explicit_model_sections(resolved_model)[blow_up_center[l]])
        blow_up_center[l] = new_locus
      end
    end

    # Conduct the blowup
    if ambient_space(resolved_model) isa NormalToricVariety
      # Toric case is easy...
      resolved_model = blow_up(resolved_model, blow_up_center; coordinate_name = exceptionals[k], nr_of_current_blow_up = k, nr_blowups_in_sequence = nr_blowups)
    else
      # Compute proper transform of center generated by anything but exceptional divisors
      filtered_center = [c for c in blow_up_center if !(c in exceptionals)]
      initial_ambient_space = ambient_space(m)
      initial_cox_ring = cox_ring(initial_ambient_space)
      initial_filtered_ideal_sheaf = ideal_sheaf(initial_ambient_space, ideal([eval_poly(l, initial_cox_ring) for l in filtered_center]))
      bd_morphism = get_attribute(blow_up_chain[1], :blow_down_morphism)
      filtered_ideal_sheaf = strict_transform(bd_morphism, initial_filtered_ideal_sheaf)
      for l in 2:k-1
        bd_morphism = get_attribute(blow_up_chain[l], :blow_down_morphism)
        filtered_ideal_sheaf = strict_transform(bd_morphism, filtered_ideal_sheaf)
      end

      # Compute strict transform of ideal sheaves appearing in blowup center
      exceptional_center = [c for c in blow_up_center if (c in exceptionals)]
      positions = [findfirst(==(l), exceptionals) for l in exceptional_center]
      exceptional_divisors = [exceptional_divisor(get_attribute(blow_up_chain[l], :blow_down_morphism)) for l in positions]
      exceptional_ideal_sheafs = [ideal_sheaf(d) for d in exceptional_divisors]
      for l in 1:length(positions)
        if positions[l] < k-1
          for m in positions[l]+1: k-1
            internal_bd_morphism = get_attribute(blow_up_chain[m], :blow_down_morphism)
            exceptional_ideal_sheafs[l] = strict_transform(internal_bd_morphism, exceptional_ideal_sheafs[l])
          end
        end
      end

      # Compute the prepared center
      prepared_center = filtered_ideal_sheaf
      if length(exceptional_ideal_sheafs) > 0
        prepared_center = prepared_center + sum(exceptional_ideal_sheafs)
      end

      # Execute the blow-up
      resolved_model = blow_up(resolved_model, prepared_center; coordinate_name = exceptionals[k])
    end

    # Remember the result
    push!(blow_up_chain, resolved_model)

  end


  # For model 1511.03209 and resolution_index = 1, we extend beyond what is currently saved as resolution in our json file.
  # Namely, we also resolve the ambient space. This is done by the following lines.
  if has_attribute(m, :arxiv_id) && resolution_index == 1 && arxiv_id(m) == "1511.03209"

    # Additional blowup 1:
    # Additional blowup 1:
    # Ambient space has the following rays
    #x: [0,0,0,-3,1]
    #y: [0,0,0,2,-1]
    #z: [0,0,0,0,1]
    # We add the ray m1: (0,0,0,1,0). This looks like y + z = 2 * m1.
    # So naively, I think of this as blowing up y^2 = z^2 = 0 and introducing the variable m1. For the strict transform, we thus do
    # y^2 -> y^2 * m1
    # z^2 -> z^2 * m1
    # y * z -> y * z * m1
    as = ambient_space(resolved_model);
    bl = domain(blow_up(as, [0,0,0,1,0], coordinate_name = "m1", 1, 3));
    f = hypersurface_equation(resolved_model);
    my_mons = collect(monomials(f));
    pos_1 = findfirst(k -> k == "y", [string(a) for a in gens(cox_ring(as))])
    pos_2 = findfirst(k -> k == "z", [string(a) for a in gens(cox_ring(as))])
    exp_list = [collect(exponents(m))[1] for m in my_mons];
    my_exps = [[k[pos_1], k[pos_2]] for k in exp_list];
    @req all(k -> isinteger(sum(k)), my_exps) "Inconsistency encountered in computation of strict transform. Please inform the authors."
    m_power = [Int(1//2*sum(a)) for a in my_exps]
    overall_factor = minimum(m_power)
    new_coeffs = collect(coefficients(f))
    new_exps = [vcat([exp_list[k], m_power[k] - overall_factor]...) for k in 1:length(exp_list)]
    my_builder = MPolyBuildCtx(cox_ring(bl))
    for a in 1:length(new_exps)
      push_term!(my_builder, new_coeffs[a], new_exps[a])
    end
    new_tate_polynomial = finish(my_builder);
    model_bl = GlobalTateModel(explicit_model_sections(resolved_model), model_section_parametrization(resolved_model), new_tate_polynomial, base_space(resolved_model), bl);
    set_attribute!(model_bl, :partially_resolved, true)
    # Additional blowup 2:
    # Additional blowup 2:
    # Ambient space has the following rays:
    # x: [0,0,0,-3,1]
    # y: [0,0,0,2,-1]
    # z: [0,0,0,0,1]
    # m1: [0,0,0,1,0]
    # We add the ray m2: (0,0,0,-2,1). This looks like 2 * x + z = 3 * m2. For the strict transform, we thus do
    # x^3 -> x^3 * m2^2
    # z^3 -> z^3 * m2
    as = ambient_space(model_bl);
    bl = domain(blow_up(as, [0,0,0,-2,1], coordinate_name = "m2", 2, 3));
    f = hypersurface_equation(model_bl);
    my_mons = collect(monomials(f));
    pos_1 = findfirst(k -> k == "x", [string(a) for a in gens(cox_ring(as))])
    pos_2 = findfirst(k -> k == "z", [string(a) for a in gens(cox_ring(as))])
    exp_list = [collect(exponents(m))[1] for m in my_mons];
    my_exps = [[k[pos_1], k[pos_2]] for k in exp_list];
    @req all(k -> isinteger(sum(k)), my_exps) "Inconsistency encountered in computation of strict transform. Please inform the authors."
    m_power = [Int(2//3 * a[1] + 1//3 * a[2]) for a in my_exps]
    overall_factor = minimum(m_power)
    new_coeffs = collect(coefficients(f))
    new_exps = [vcat([exp_list[k], m_power[k] - overall_factor]...) for k in 1:length(exp_list)]
    my_builder = MPolyBuildCtx(cox_ring(bl))
    for a in 1:length(new_exps)
      push_term!(my_builder, new_coeffs[a], new_exps[a])
    end
    new_tate_polynomial = finish(my_builder);
    model_bl2 = GlobalTateModel(explicit_model_sections(model_bl), model_section_parametrization(model_bl), new_tate_polynomial, base_space(model_bl), bl);
    set_attribute!(model_bl2, :partially_resolved, true)
    # Additional blowup 3:
    # Additional blowup 3:
    # Ambient space has the following rays:
    # x: [0,0,0,-3,1]
    # y: [0,0,0,2,-1]
    # z: [0,0,0,0,1]
    # m1: [0,0,0,1,0]
    # m2: [0,0,0,-2,1]
    # We add the ray m3: (0,0,0,-1,1). This looks like m2 + z = 2 * m3. For the strict transform, we thus do
    # m2^2 -> m2^2 * m3
    # z^2 -> z^2 * m3
    as = ambient_space(model_bl2);
    bl = domain(blow_up(as, [0,0,0,-1,1], coordinate_name = "m3", 3, 3));
    f = hypersurface_equation(model_bl2);
    my_mons = collect(monomials(f));
    pos_1 = findfirst(k -> k == "m2", [string(a) for a in gens(cox_ring(as))])
    pos_2 = findfirst(k -> k == "z", [string(a) for a in gens(cox_ring(as))])
    exp_list = [collect(exponents(m))[1] for m in my_mons];
    my_exps = [[k[pos_1], k[pos_2]] for k in exp_list];
    @req all(k -> isinteger(sum(k)), my_exps) "Inconsistency encountered in computation of strict transform. Please inform the authors."
    m_power = [Int(1//2 * sum(a)) for a in my_exps]
    overall_factor = minimum(m_power)
    new_coeffs = collect(coefficients(f))
    new_exps = [vcat([exp_list[k], m_power[k] - overall_factor]...) for k in 1:length(exp_list)]
    my_builder = MPolyBuildCtx(cox_ring(bl))
    for a in 1:length(new_exps)
      push_term!(my_builder, new_coeffs[a], new_exps[a])
    end
    new_tate_polynomial = finish(my_builder);
    model_bl3 = GlobalTateModel(explicit_model_sections(model_bl2), model_section_parametrization(model_bl2), new_tate_polynomial, base_space(model_bl2), bl);
    set_attribute!(model_bl3, :partially_resolved, true)

    # We confirm that after these steps, we achieve what we desire.
    @req is_smooth(ambient_space(model_bl3)) "Ambient space not yet smooth. Please inform the authors!"
    @req is_homogeneous(hypersurface_equation(model_bl3)) "Strict transform is not homogeneous. Please inform the authors!"
    return model_bl3
    
  end

  return resolved_model
end
