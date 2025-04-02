#####################################################
# 1: Tune a Weierstrass model
#####################################################

@doc raw"""
    function tune(w::WeierstrassModel, input_sections::Dict{String, <:Any}; completeness_check::Bool = true)

Tune a Weierstrass model by fixing a special choice for the model sections.
Note that it is in particular possible to set a section to zero. We anticipate
that people might want to be able to come back from this by assigning a non-trivial
value to a section that was previously tuned to zero. This is why we keep such
trivial sections and do not delete them, say from `explicit_model_sections`
or `classes_of_model_sections`.

# Examples
```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> w = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, defining_classes = Dict("b" => b), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)

julia> x1, x2, x3 = gens(cox_ring(base_space(w)))
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3

julia> my_choice = Dict("b" => x2, "c2" => x1^6)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "c2" => x1^6
  "b"  => x2

julia> tuned_w = tune(w, my_choice)
Weierstrass model over a concrete base

julia> explicit_model_sections(tuned_w)["c2"]
x1^6

julia> x1, x2, x3 = gens(cox_ring(base_space(tuned_w)))
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3

julia> my_choice2 = Dict("b" => x2, "c2" => zero(parent(x1)))
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "c2" => 0
  "b"  => x2

julia> tuned_w2 = tune(tuned_w, my_choice2)
Weierstrass model over a concrete base

julia> is_zero(explicit_model_sections(tuned_w2)["c2"])
true

julia> x1, x2, x3 = gens(cox_ring(base_space(tuned_w2)))
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3

julia> my_choice3 = Dict("b" => x2, "c2" => x1^6)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  "c2" => x1^6
  "b"  => x2

julia> tuned_w3 = tune(tuned_w2, my_choice3)
Weierstrass model over a concrete base

julia> is_zero(explicit_model_sections(tuned_w3)["c2"])
false
```
"""
function tune(w::WeierstrassModel, input_sections::Dict{String, <:Any}; completeness_check::Bool = true)
  # Consistency checks
  @req base_space(w) isa NormalToricVariety "Currently, tuning is only supported for models over concrete toric bases"
  isempty(input_sections) && return w
  secs_names = tunable_sections(w)
  tuned_secs_names = collect(keys(input_sections))
  @req all(in(secs_names), tuned_secs_names) "Provided section names are not among the tunable sections of the model"

  # 0. Prepare for computation by setting up some information
  explicit_secs = deepcopy(explicit_model_sections(w))
  def_secs_param = deepcopy(model_section_parametrization(w))
  weierstrass_sections = ["f", "g"]

  # 1. Tune model sections different from Weierstrass sections
  for x in setdiff(tuned_secs_names, weierstrass_sections)
    @req parent(input_sections[x]) == parent(explicit_model_sections(w)[x]) "Parent mismatch between given and existing model section"
    if is_zero(input_sections[x]) == false
      @req degree(input_sections[x]) == divisor_class(classes_of_model_sections(w)[x]) "Degree mismatch between given and existing model section"
    end
    explicit_secs[x] = input_sections[x]
  end

  # 2. Use model sections to reevaluate the Weierstrass sections via their known parametrization
  parametrization_keys = collect(keys(def_secs_param))
  if !isempty(parametrization_keys) && !isempty(secs_names)
    R = parent(def_secs_param[parametrization_keys[1]])
    S = parent(explicit_secs[secs_names[1]])
    vars = [string(k) for k in symbols(R)]
    images = [k in secs_names ? explicit_secs[k] : k == "Kbar" ? eval_poly("0", S) : eval_poly(k, S) for k in vars]
    map = hom(R, S, images)
    for section in weierstrass_sections
      haskey(def_secs_param, section) && (explicit_secs[section] = map(eval_poly(string(def_secs_param[section]), R)))
    end
  end

  # 3. Does the user want to set some Weierstrass sections? If so, overwrite existing choice with desired value.
  for sec in weierstrass_sections
    if haskey(input_sections, sec)
      @req parent(input_sections[sec]) == parent(explicit_model_sections(w)[sec]) "Parent mismatch between given and existing Weierstrass section"
      if is_zero(input_sections[sec]) == false
        @req degree(input_sections[sec]) == divisor_class(classes_of_model_sections(w)[sec]) "Degree mismatch between given and existing Weierstrass section"
      end
      explicit_secs[sec] = input_sections[sec]
      delete!(def_secs_param, sec)
    end
  end

  # 4. There could be unused model sections...
  if !isempty(parametrization_keys)
    polys = [eval_poly(string(def_secs_param[section]), R) for section in weierstrass_sections if haskey(def_secs_param, section)]
    all_appearing_monomials = vcat([collect(monomials(p)) for p in polys]...)
    all_appearing_exponents = [collect(exponents(m))[1] for m in all_appearing_monomials]
    potentially_redundant_sections = gens(R)
    for k in 1:length(potentially_redundant_sections)
      string(potentially_redundant_sections[k]) in weierstrass_sections && continue
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
  resulting_model = weierstrass_model(base_space(w), explicit_secs, def_secs_param; completeness_check)

  # 7. Copy the classes of model sections
  new_classes_of_model_sections = Dict{String, ToricDivisorClass}()
  for (key, value) in classes_of_model_sections(w)
    m = divisor_class(value).coeff
    @req nrows(m) == 1 "Encountered inconsistency"
    new_classes_of_model_sections[key] = toric_divisor_class(base_space(resulting_model), m[1, :])
  end
  set_attribute!(resulting_model, :classes_of_model_sections => new_classes_of_model_sections)

  # 8. Return the model
  return resulting_model
end
