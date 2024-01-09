#####################################################
# 1: Tune a Weierstrass model
#####################################################

@doc raw"""
    tune(w::WeierstrassModel, special_section_choices::Dict{String, <:MPolyRingElem}; completeness_check::Bool = true)

Tune a Weierstrass model. For this, a dictionary is presented as second argument,
which specifies the desired specific choices for the Weierstrass sections.

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

julia> my_choice = Dict("f" => x1^12)
Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 1 entry:
  "f" => x1^12

julia> tuned_w = tune(w, my_choice)
Weierstrass model over a concrete base

julia> weierstrass_section_f(tuned_w) == my_choice["f"]
true
```
"""
function tune(w::WeierstrassModel, special_section_choices::Dict{String, <:MPolyRingElem}; completeness_check::Bool = true)
  @req (typeof(base_space(w)) <: NormalToricVariety) "Currently, tuning is only possible for models over concrete toric bases"
  isempty(special_section_choices) && return w
  f = weierstrass_section_f(w)
  g = weierstrass_section_g(w)
  if haskey(special_section_choices, "f")
    @req parent(special_section_choices["f"]) == parent(weierstrass_section_f(w)) "Parent mismatch between given and existing Weierstrass section f"
    @req degree(special_section_choices["f"]) == degree(weierstrass_section_f(w)) "Parent mismatch between given and existing Weierstrass section f"
    f = special_section_choices["f"]
  end
  if haskey(special_section_choices, "g")
    @req parent(special_section_choices["g"]) == parent(weierstrass_section_g(w)) "Parent mismatch between given and existing Weierstrass section g"
    @req degree(special_section_choices["g"]) == degree(weierstrass_section_g(w)) "Parent mismatch between given and existing Weierstrass section g"
    g = special_section_choices["g"]
  end
  return weierstrass_model(base_space(w), f, g; completeness_check)
end
