#####################################################
# 1: Fiber analysis
#####################################################

@doc raw"""
    analyze_fibers(model::GlobalTateModel, centers::Vector{<:Vector{<:Integer}})

Determine the fiber of a (singular) global Tate model over a particular base locus.
```
"""
function analyze_fibers(model::GlobalTateModel, centers::Vector{<:Vector{<:Integer}})
  
  # This method only works if the model is defined over a toric variety over toric scheme
  @req typeof(base_space(model)) <: NormalToricVariety "Analysis of fibers currently only supported for toric scheme/variety as base space"
  
  # Ideal of the defining polynomial
  hypersurface_ideal = ideal([tate_polynomial(model)])
  
  # Toric ambient space
  tas = ambient_space(model)
  
  # Various important ideals
  irr = irrelevant_ideal(tas);
  sri = stanley_reisner_ideal(tas);
  lin = ideal_of_linear_relations(tas);
  
  # Singular loci
  sing_loc = singular_loci(model)
  
  # Pick out the singular loci that are more singular than an I_1
  # Then keep only the locus and not the extra info about it
  interesting_singular_loci = map(tup -> tup[1], filter(locus -> locus[2][3] > 1, sing_loc))
  
  # This is a kludge to map polynomials on the base into the ambient space, and should be fixed once the ambient space constructors supply such a map
  base_coords = parent(gens(interesting_singular_loci[1])[1])
  ambient_coords = parent(tate_polynomial(model))
  base_to_ambient_ring_map = hom(base_coords, ambient_coords, gens(ambient_coords)[1:end-3])
  
  # Resolved model
  strict_transform, exceptionals, crepant, res_irr, res_sri, res_lin, res_S, res_S_gens, res_ring_map = _blowup_global_sequence(hypersurface_ideal, centers, irr, sri, lin)
  if !crepant
      @warn "The given sequence of blowups is not crepant"
  end
  
  loci_fiber_intersections = Tuple{MPolyIdeal{QQMPolyRingElem}, Vector{Tuple{Tuple{Int64, Int64}, Vector{MPolyIdeal{QQMPolyRingElem}}}}}[]
  for locus in interesting_singular_loci
    # Currently have to get the ungraded ideal generators by hand using .f
    ungraded_locus = ideal(map(gen -> base_to_ambient_ring_map(gen).f, gens(locus)))
    
    # Potential components of the fiber over this locus
    # For now, we only consider the associated prime ideal,
    # but we may later want to actually consider the primary ideals
    potential_components = map(pair -> pair[2], primary_decomposition(strict_transform + res_ring_map(ungraded_locus)))
    
    # Filter out the trivial loci among the potential components
    components = filter(component -> _is_nontrivial(component, res_irr), potential_components)
    
    # Check the pairwise intersections of the components
    intersections = Tuple{Tuple{Int64, Int64}, Vector{MPolyIdeal{QQMPolyRingElem}}}[]
    for i in 1:length(components) - 1
      for j in i + 1:length(components)
        intersection = filter(candidate_locus -> _is_nontrivial(candidate_locus, res_irr), map(pair -> pair[2], primary_decomposition(components[i] + components[j])))
        push!(intersections, ((i, j), intersection))
      end
    end
    
    push!(loci_fiber_intersections, (ungraded_locus, intersections))
  end
  
  return loci_fiber_intersections

end



#####################################################
# 2: Resolve a Tate model via blowup
#####################################################

@doc raw"""
    blow_up(t::GlobalTateModel, ideal_gens::Vector{String}; coordinate_name::String = "e")

Resolve a global Tate model by blowing up a locus in the ambient space.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, model_sections = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> blow_up(t, ["x", "y", "x1"]; coordinate_name = "e1")
Partially resolved global Tate model over a concrete base
```
"""
function blow_up(t::GlobalTateModel, ideal_gens::Vector{String}; coordinate_name::String = "e")
  R = cox_ring(ambient_space(t))
  I = ideal([eval_poly(k, R) for k in ideal_gens])
  return blow_up(t, I; coordinate_name = coordinate_name)
end

function blow_up(t::GlobalTateModel, I::MPolyIdeal; coordinate_name::String = "e")
  
  # This method only works if the model is defined over a toric variety over toric scheme
  @req typeof(base_space(t)) <: NormalToricVariety "Blowups of Tate models are currently only supported for toric bases"
  @req typeof(ambient_space(t)) <: NormalToricVariety "Blowups of Tate models are currently only supported for toric ambient spaces"

  # Compute the new ambient_space
  bd = blow_up(ambient_space(t), I; coordinate_name = coordinate_name)
  new_ambient_space = domain(bd)

  # Compute the new base
  # FIXME: THIS WILL IN GENERAL BE WRONG! IN PRINCIPLE, THE ABOVE ALLOWS TO BLOW UP THE BASE AND THE BASE ONLY.
  # FIXME: We should save the projection \pi from the ambient space to the base space.
  # FIXME: This is also ties in with the model sections to be saved, see below. Should the base change, so do these sections...
  new_base = base_space(t)

  # Prepare for the computation of the strict transform of the tate polynomial
  # FIXME: This assume that I is generated by indeterminates! Very special!
  S = cox_ring(new_ambient_space)
  _e = eval_poly(coordinate_name, S)
  images = MPolyRingElem[]
  for v in gens(S)
    v == _e && continue
    if string(v) in [string(k) for k in gens(I)]
      push!(images, v * _e)
    else
      push!(images, v)
    end
  end
  ring_map = hom(base_ring(I), S, images)
  total_transform = ring_map(ideal([tate_polynomial(t)]))
  exceptional_ideal = total_transform + ideal([_e])
  strict_transform, exceptional_factor = saturation_with_index(total_transform, exceptional_ideal)
  new_pt = gens(strict_transform)[1]

  # Extract the old Tate sections, which do not change.
  ais = [tate_section_a1(t), tate_section_a2(t), tate_section_a3(t), tate_section_a4(t), tate_section_a6(t)]

  # Construct the new model
  # This is not really a Tate model any more, as the hypersurface equation is merely a strict transform of a Tate polynomial.
  # Change/Fix? We may want to provide not only output that remains true forever but also output, while the internals may change?
  model = GlobalTateModel(ais[1], ais[2], ais[3], ais[4], ais[5], new_pt, base_space(t), new_ambient_space)

  # Set attributes
  set_attribute!(model, :base_fully_specified, true)
  set_attribute!(model, :partially_resolved, true)
  if has_attribute(t, :explicit_model_sections)
    set_attribute!(model, :explicit_model_sections => get_attribute(t, :explicit_model_sections))
  end
  return model
end
