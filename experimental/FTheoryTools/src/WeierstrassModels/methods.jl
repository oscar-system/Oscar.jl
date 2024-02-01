#####################################################
# 1: Tune a Weierstrass model
#####################################################

@doc raw"""
    function tune(w::WeierstrassModel, input_sections::Dict{String, <:Any}; completeness_check::Bool = true)

Tune a Weierstrass model  by fixing a special choice for the model sections.

# Examples
```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> w = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, model_sections = Dict("b" => b), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)

julia> x1, x2, x3 = gens(cox_ring(base_space(w)))
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3

julia> my_choice = Dict("f" => x1^12, "b" => x2, "c2" => x1^6)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 3 entries:
  "c2" => x1^6
  "f"  => x1^12
  "b"  => x2

julia> tuned_w = tune(w, my_choice)
Weierstrass model over a concrete base

julia> weierstrass_section_f(tuned_w) == my_choice["f"]
true
```
"""
function tune(w::WeierstrassModel, input_sections::Dict{String, <:Any}; completeness_check::Bool = true)
  # Consistency checks
  @req base_space(w) isa NormalToricVariety "Currently, tuning is only supported for models over concrete toric bases"
  isempty(input_sections) && return w
  secs_names = collect(keys(explicit_model_sections(w)))
  tuned_secs_names = collect(keys(input_sections))
  @req all(x -> x in secs_names, tuned_secs_names) "Provided section name not recognized"

  # 0. Prepare for computation by setting up some information
  explicit_secs = deepcopy(explicit_model_sections(w))
  def_secs_param = deepcopy(defining_section_parametrization(w))
  weierstrass_sections = ["f", "g"]

  # 1. Tune model sections different from Weierstrass sections
  for x in setdiff(tuned_secs_names, weierstrass_sections)
    section_parent = parent(input_sections[x])
    section_degree = degree(input_sections[x])
    @req section_parent == parent(explicit_model_sections(w)[x]) "Parent mismatch between given and existing model section"
    @req section_degree == degree(explicit_model_sections(w)[x]) "Degree mismatch between given and existing model section"
    explicit_secs[x] = input_sections[x]
  end

  # 2. Use model sections to reevaluate the Weierstrass sections via their known parametrization
  parametrization_keys = collect(keys(def_secs_param))
  if !isempty(parametrization_keys) && !isempty(secs_names)
    R = parent(def_secs_param[parametrization_keys[1]])
    S = parent(explicit_secs[secs_names[1]])
    vars = [string(k) for k in gens(R)]
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
      @req degree(input_sections[sec]) == degree(explicit_model_sections(w)[sec]) "Degree mismatch between given and existing Weierstrass section"
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
  
  # 5. Build the new model
  return weierstrass_model(base_space(w), explicit_secs, def_secs_param; completeness_check)
end
