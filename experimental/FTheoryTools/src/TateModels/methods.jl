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


#####################################################
# 2: Tune a Tate model
#####################################################

@doc raw"""
    tune(t::GlobalTateModel, input_sections::Dict{String, <:Any}; completeness_check::Bool = true)

Tune a Tate model by fixing a special choice for the model sections.
Note that it is in particular possible to set a section to zero. We anticipate
that people might want to be able to come back from this by assigning a non-trivial
value to a section that was previously tuned to zero. This is why we keep such
trivial sections and do not delete them, say from `explicit_model_sections`
or `classes_of_model_sections`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> x1, x2, x3, x4 = gens(cox_ring(base_space(t)))
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4

julia> my_choice = Dict("a1" => x1^4, "w" => x2 - x3)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "w"  => x2 - x3
  "a1" => x1^4

julia> tuned_t = tune(t, my_choice)
Global Tate model over a concrete base

julia> tate_section_a1(tuned_t) == x1^4
true

julia> x1, x2, x3, x4 = gens(cox_ring(base_space(tuned_t)))
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4

julia> my_choice2 = Dict("a1" => zero(parent(x1)), "w" => x2 - x3)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "w"  => x2 - x3
  "a1" => 0

julia> tuned_t2 = tune(tuned_t, my_choice2)
Global Tate model over a concrete base

julia> is_zero(explicit_model_sections(tuned_t2)["a1"])
true

julia> x1, x2, x3, x4 = gens(cox_ring(base_space(tuned_t2)))
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4

julia> my_choice3 = Dict("a1" => x1^4, "w" => x2 - x3)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "w"  => x2 - x3
  "a1" => x1^4

julia> tuned_t3 = tune(tuned_t2, my_choice3)
Global Tate model over a concrete base

julia> is_zero(explicit_model_sections(tuned_t3)["a1"])
false
```
"""
function tune(t::GlobalTateModel, input_sections::Dict{String, <:Any}; completeness_check::Bool = true)
  # Consistency checks
  @req base_space(t) isa NormalToricVariety "Currently, tuning is only supported for models over concrete toric bases"
  isempty(input_sections) && return t
  secs_names = tunable_sections(t)
  tuned_secs_names = collect(keys(input_sections))
  @req all(in(secs_names), tuned_secs_names) "Provided section names are not among the tunable sections of the model"

  # 0. Prepare for computation by setting up some information
  explicit_secs = deepcopy(explicit_model_sections(t))
  def_secs_param = deepcopy(model_section_parametrization(t))
  tate_sections = ["a1", "a2", "a3", "a4", "a6"]

  # 1. Tune model sections different from Tate sections
  for x in setdiff(tuned_secs_names, tate_sections)
    @req parent(input_sections[x]) == parent(explicit_model_sections(t)[x]) "Parent mismatch between given and existing model section"
    if is_zero(input_sections[x]) == false
      @req degree(input_sections[x]) == divisor_class(classes_of_model_sections(t)[x]) "Degree mismatch between given and existing model section"
    end
    explicit_secs[x] = input_sections[x]
  end

  # 2. Use model sections to reevaluate the Tate sections via their known parametrization
  parametrization_keys = collect(keys(def_secs_param))
  if !isempty(parametrization_keys) && !isempty(secs_names)
    R = parent(def_secs_param[parametrization_keys[1]])
    S = parent(explicit_secs[secs_names[1]])
    vars = [string(k) for k in symbols(R)]
    images = [k in secs_names ? explicit_secs[k] : k == "Kbar" ? eval_poly("0", S) : eval_poly(k, S) for k in vars]
    map = hom(R, S, images)
    for section in tate_sections
      haskey(def_secs_param, section) && (explicit_secs[section] = map(eval_poly(string(def_secs_param[section]), R)))
    end
  end

  # 3. Does the user want to set some Tate sections? If so, overwrite existing choice with desired value.
  for sec in tate_sections
    if haskey(input_sections, sec)
      @req parent(input_sections[sec]) == parent(explicit_model_sections(t)[sec]) "Parent mismatch between given and existing Tate section"
      if is_zero(input_sections[sec]) == false
        @req degree(input_sections[sec]) == divisor_class(classes_of_model_sections(t)[sec]) "Degree mismatch between given and existing Tate section"
      end
      explicit_secs[sec] = input_sections[sec]
      delete!(def_secs_param, sec)
    end
  end

  # 4. There could be unused model sections...
  if !isempty(parametrization_keys)
    polys = [eval_poly(string(def_secs_param[section]), R) for section in tate_sections if haskey(def_secs_param, section)]
    all_appearing_monomials = vcat([collect(monomials(p)) for p in polys]...)
    all_appearing_exponents = [collect(exponents(m))[1] for m in all_appearing_monomials]
    potentially_redundant_sections = gens(R)
    for k in 1:length(potentially_redundant_sections)
      string(potentially_redundant_sections[k]) in tate_sections && continue
      is_used = any(all_appearing_exponents[l][k] != 0 for l in 1:length(all_appearing_exponents))
      is_used || delete!(explicit_secs, string(potentially_redundant_sections[k]))
    end
  end
  
  # 5. After removing some sections, we must go over the parametrization again and adjust the ring in which the parametrization is given.
  if !isempty(def_secs_param)
    naive_vars = string.(gens(parent(first(values(def_secs_param)))))
    filtered_vars = filter(x -> haskey(explicit_secs, x), naive_vars)
    desired_ring, _ = polynomial_ring(QQ, filtered_vars, cached = false)
    for (key, value) in def_secs_param
      def_secs_param[key] = eval_poly(string(value), desired_ring)
    end
  end
  
  # 6. Build the new model
  resulting_model = global_tate_model(base_space(t), explicit_secs, def_secs_param; completeness_check)

  # 7. Copy the classes of model sections, but only of those sections that are used!
  new_classes_of_model_sections = Dict{String, ToricDivisorClass}()
  for key in keys(explicit_model_sections(resulting_model))
    m = divisor_class(classes_of_model_sections(t)[key]).coeff
    @req nrows(m) == 1 "Encountered inconsistency"
    new_classes_of_model_sections[key] = toric_divisor_class(base_space(resulting_model), m[1, :])
  end
  set_attribute!(resulting_model, :classes_of_model_sections => new_classes_of_model_sections)

  # 8. Return the model
  return resulting_model

end
