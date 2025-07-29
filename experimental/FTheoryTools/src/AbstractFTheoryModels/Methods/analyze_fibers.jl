# Check if an ideal/subvariety is nontrivial
_is_nontrivial(id::MPolyIdeal{T}, irr::MPolyIdeal{T}) where {T<:MPolyRingElem} = !is_one(id) && !is_one(saturation(id, irr))

@doc raw"""
    analyze_fibers(model::GlobalTateModel, centers::Vector{<:Vector{<:Integer}})

Determine the fiber of a (singular) global Tate model over a particular base locus.

!!! warning
    This method may run for very long time and is currently not tested as part of the regular OSCAR CI due to its excessive run times.
"""
function analyze_fibers(model::GlobalTateModel, centers::Vector{<:Vector{<:Integer}}; rng::AbstractRNG = Random.default_rng())
  
  # This method only works if the model is defined over a toric variety over toric scheme
  @req base_space(model) isa NormalToricVariety "Analysis of fibers currently only supported for toric variety as base space"
  
  # Ideal of the defining polynomial
  hypersurface_ideal = ideal([tate_polynomial(model)])
  
  # Toric ambient space
  tas = ambient_space(model)
  
  # Various important ideals
  irr = irrelevant_ideal(tas);
  sri = stanley_reisner_ideal(tas);
  lin = ideal_of_linear_relations(tas);
  
  # Singular loci
  sing_loc = singular_loci(model; rng)
  
  # Pick out the singular loci that are more singular than an I_1
  # Then keep only the locus and not the extra info about it
  interesting_singular_loci = map(first, filter(locus -> locus[2][3] > 1, sing_loc))
  
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
    # Currently have to get the ungraded ideal generators by hand using lift
    ungraded_locus = ideal(map(gen -> lift(base_to_ambient_ring_map(gen)), gens(locus)))
    
    # Potential components of the fiber over this locus
    # For now, we only consider the associated prime ideal,
    # but we may later want to actually consider the primary ideals
    potential_components = map(last, primary_decomposition(strict_transform + res_ring_map(ungraded_locus)))
    
    # Filter out the trivial loci among the potential components
    components = filter(component -> _is_nontrivial(component, res_irr), potential_components)
    
    # Check the pairwise intersections of the components
    intersections = Tuple{Tuple{Int64, Int64}, Vector{MPolyIdeal{QQMPolyRingElem}}}[]
    for i in 1:length(components) - 1
      for j in i + 1:length(components)
        intersection = filter(candidate_locus -> _is_nontrivial(candidate_locus, res_irr), map(last, primary_decomposition(components[i] + components[j])))
        push!(intersections, ((i, j), intersection))
      end
    end
    
    push!(loci_fiber_intersections, (ungraded_locus, intersections))
  end
  
  return loci_fiber_intersections

end
