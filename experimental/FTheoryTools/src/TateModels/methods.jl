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
# 2: Tune a Tate model
#####################################################

@doc raw"""
    tune(t::GlobalTateModel, special_ai_choices::Dict{String, <:Any}; completeness_check::Bool = true)

Tune a Tate model by fixing a special choice of the Tate sections.
This choice is provided by a dictionary, to be provided as second argument.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, model_sections = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> x1, x2, x3, x4 = gens(cox_ring(base_space(t)))
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4

julia> my_choice = Dict("a1" => x1^4)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 1 entry:
  "a1" => x1^4

julia> tuned_t = tune(t, my_choice)
Global Tate model over a concrete base

julia> tate_section_a1(tuned_t) == x1^4
true
```
"""
function tune(t::GlobalTateModel, special_ai_choices::Dict{String, <:Any}; completeness_check::Bool = true)
  @req base_space(t) isa NormalToricVariety "Currently, tuning is only supported for models over concrete toric bases"
  isempty(special_ai_choices) && return t
  ais = [tate_section_a1(t), tate_section_a2(t), tate_section_a3(t), tate_section_a4(t), tate_section_a6(t)]
  if haskey(special_ai_choices, "a1")
    @req parent(special_ai_choices["a1"]) == parent(tate_section_a1(t)) "Parent mismatch between given and existing Tate section a1"
    @req degree(special_ai_choices["a1"]) == degree(tate_section_a1(t)) "Parent mismatch between given and existing Tate section a1"
    ais[1] = special_ai_choices["a1"]
  end
  if haskey(special_ai_choices, "a2")
    @req parent(special_ai_choices["a2"]) == parent(tate_section_a2(t)) "Parent mismatch between given and existing Tate section a2"
    @req degree(special_ai_choices["a2"]) == degree(tate_section_a2(t)) "Parent mismatch between given and existing Tate section a2"
    ais[2] = special_ai_choices["a2"]
  end
  if haskey(special_ai_choices, "a3")
    @req parent(special_ai_choices["a3"]) == parent(tate_section_a3(t)) "Parent mismatch between given and existing Tate section a3"
    @req degree(special_ai_choices["a3"]) == degree(tate_section_a3(t)) "Parent mismatch between given and existing Tate section a3"
    ais[3] = special_ai_choices["a3"]
  end
  if haskey(special_ai_choices, "a4")
    @req parent(special_ai_choices["a4"]) == parent(tate_section_a4(t)) "Parent mismatch between given and existing Tate section a4"
    @req degree(special_ai_choices["a4"]) == degree(tate_section_a4(t)) "Parent mismatch between given and existing Tate section a4"
    ais[4] = special_ai_choices["a4"]
  end
  if haskey(special_ai_choices, "a6")
    @req parent(special_ai_choices["a6"]) == parent(tate_section_a6(t)) "Parent mismatch between given and existing Tate section a6"
    @req degree(special_ai_choices["a6"]) == degree(tate_section_a6(t)) "Parent mismatch between given and existing Tate section a6"
    ais[5] = special_ai_choices["a6"]
  end
  return global_tate_model(base_space(t), ais; completeness_check)
end
