#####################################################################
# 1: Constructors with toric variety as base
#####################################################################

@doc raw"""
    weierstrass_model(base::NormalToricVariety; completeness_check::Bool = true)

This method constructs a Weierstrass model over a given toric base
3-fold. The Weierstrass sections ``f`` and ``g`` are taken with (pseudo)random
coefficients.

# Examples
```jldoctest
julia> w = weierstrass_model(sample_toric_variety(); completeness_check = false)
Weierstrass model over a concrete base
```
"""
function weierstrass_model(base::NormalToricVariety; completeness_check::Bool = true)
  (f, g) = _weierstrass_sections(base)
  return weierstrass_model(base, f, g; completeness_check = completeness_check)
end


@doc raw"""
    weierstrass_model(base::NormalToricVariety, f::MPolyRingElem, g::MPolyRingElem; completeness_check::Bool = true)

This method operates analogously to `weierstrass_model(base::NormalToricVarietyType)`.
The only difference is that the Weierstrass sections ``f`` and ``g`` can be specified with non-generic values.

# Examples
```jldoctest
julia> base = sample_toric_variety()
Normal toric variety

julia> f = generic_section(anticanonical_bundle(base)^4);

julia> g = generic_section(anticanonical_bundle(base)^6);

julia> w = weierstrass_model(base, f, g; completeness_check = false)
Weierstrass model over a concrete base
```
"""
function weierstrass_model(base::NormalToricVariety, f::MPolyRingElem, g::MPolyRingElem; completeness_check::Bool = true)
  return weierstrass_model(base, Dict("f" => f, "g" => g), Dict{String, MPolyRingElem}(); completeness_check = completeness_check)
end

function weierstrass_model(base::NormalToricVariety,
                           explicit_model_sections::Dict{String, <: Union{MPolyRingElem, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}},
                           defining_section_parametrization::Dict{String, <:MPolyRingElem};
                           completeness_check::Bool = true)
  vs = collect(values(explicit_model_sections))
  @req all(x -> parent(x) == cox_ring(base), vs) "All model sections must reside in the Cox ring of the base toric variety"
  @req haskey(explicit_model_sections, "f") "Weierstrass section f must be specified"
  @req haskey(explicit_model_sections, "g") "Weierstrass section g must be specified"
  vs2 = collect(keys(defining_section_parametrization))
  @req all(x -> x in ["f", "g"], vs2) "Only the Weierstrass sections f, g must be parametrized"

  gens_base_names = [string(g) for g in gens(cox_ring(base))]
  if ("x" in gens_base_names) || ("y" in gens_base_names) || ("z" in gens_base_names)
    @vprint :FTheoryModelPrinter 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  
  if completeness_check
    @req is_complete(base) "Base space must be complete"
  end
  
  # construct the ambient space
  fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
  set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])
  D1 = 2 * anticanonical_divisor_class(base)
  D2 = 3 * anticanonical_divisor_class(base)
  D3 = trivial_divisor_class(base)
  ambient_space = _ambient_space(base, fiber_ambient_space, [D1, D2, D3])
  
  # construct the model
  pw = _weierstrass_polynomial(explicit_model_sections["f"], explicit_model_sections["g"], cox_ring(ambient_space))
  model = WeierstrassModel(explicit_model_sections, defining_section_parametrization, pw, base, ambient_space)
  set_attribute!(model, :partially_resolved, false)
  return model
end


#####################################################################
# 2: Constructors with scheme as base
#####################################################################

# Yet to come...
# This requires that the ai are stored as sections of the anticanonical bundle, and not "just" polynomials.
# -> Types to be generalized then.


#####################################################################
# 3: Constructor with toric base attempting to represent moduli space
#####################################################################

@doc raw"""
    weierstrass_model(auxiliary_base_ring::MPolyRing, auxiliary_base_grading::Matrix{Int64}, d::Int, weierstrass_f::MPolyRingElem, weierstrass_g::MPolyRingElem)

This method constructs a Weierstrass model over a base space that is not
fully specified.

Note that many studies in the literature use the class of the anticanonical bundle
in their analysis. We anticipate this by adding this class as a variable of the
auxiliary base space, unless the user already provides this grading. Our convention
is that the first grading refers to Kbar and that the homogeneous variable corresponding
to this class carries the name "Kbar".

The following example illustrates this approach.

# Examples
```jldoctest
julia> auxiliary_base_ring, (f, g, Kbar, v) = QQ["f", "g", "Kbar", "u"]
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[f, g, Kbar, u])

julia> auxiliary_base_grading = [4 6 1 0]
1×4 Matrix{Int64}:
 4  6  1  0

julia> w = weierstrass_model(auxiliary_base_ring, auxiliary_base_grading, 3, f, g)
Assuming that the first row of the given grading is the grading under Kbar

Weierstrass model over a not fully specified base
```
"""
function weierstrass_model(auxiliary_base_ring::MPolyRing, auxiliary_base_grading::Matrix{Int64}, d::Int, weierstrass_f::MPolyRingElem, weierstrass_g::MPolyRingElem)

  # Execute consistency checks
  gens_base_names = [string(g) for g in gens(auxiliary_base_ring)]
  @req ((parent(weierstrass_f) == auxiliary_base_ring) && (parent(weierstrass_g) == auxiliary_base_ring)) "All Weierstrass sections must reside in the provided auxiliary base ring"
  @req d > 0 "The dimension of the base space must be positive"
  if ("x" in gens_base_names) || ("y" in gens_base_names) || ("z" in gens_base_names)
    @vprint :FTheoryModelPrinter 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  
  # Inform about the assume Kbar grading
  @vprint :FTheoryModelPrinter 0 "Assuming that the first row of the given grading is the grading under Kbar\n\n"
  
  # Compute Weierstrass polynomial
  (S, auxiliary_base_space, auxiliary_ambient_space) = _construct_generic_sample(auxiliary_base_grading, gens_base_names, d)
  ring_map = hom(parent(weierstrass_f), S, gens(S)[1:ngens(parent(weierstrass_f))])
  (f, g) = [ring_map(weierstrass_f), ring_map(weierstrass_g)]
  pw = _weierstrass_polynomial(f, g, coordinate_ring(auxiliary_ambient_space))

  # Compute explicit model sections
  explicit_model_sections = Dict("f" => f, "g" => g)
  section_candidates = gens(S)
  for k in section_candidates
    haskey(explicit_model_sections, string(k)) || (explicit_model_sections[string(k)] = k)
  end

  # Compute defining_section_parametrization
  defining_section_parametrization = Dict{String, MPolyRingElem}()
  vars_S = [string(k) for k in gens(S)]
  if !("f" in vars_S) || (f != eval_poly("f", parent(f)))
    defining_section_parametrization["f"] = f
  end
  if !("g" in vars_S) || (g != eval_poly("g", parent(g)))
    defining_section_parametrization["g"] = g
  end
  
  model = WeierstrassModel(explicit_model_sections, defining_section_parametrization, pw, auxiliary_base_space, auxiliary_ambient_space)
  set_attribute!(model, :partially_resolved, false)
  return model
end


#####################################################################
# 4: Display
#####################################################################

function Base.show(io::IO, w::WeierstrassModel)
  properties_string = String[]
  if is_partially_resolved(w)
    push!(properties_string, "Partially resolved Weierstrass model over a")
  else
    push!(properties_string, "Weierstrass model over a")
  end
  if is_base_space_fully_specified(w)
    push!(properties_string, "concrete base")
  else
    push!(properties_string, "not fully specified base")
  end
  if has_model_description(w)
    push!(properties_string, "-- " * string(get_attribute(w, :model_description)))
    if has_model_parameters(w)
      push!(properties_string, "with parameter values (" * join(["$key = $(string(val))" for (key, val) in model_parameters(t)], ", ") * ")")
    end
  end
  if has_arxiv_id(w)
    push!(properties_string, "based on arXiv paper " * string(get_attribute(w, :arxiv_id)))
  end
  if has_arxiv_model_equation_number(w)
    push!(properties_string, "Eq. (" * string(get_attribute(w, :arxiv_model_equation_number)) * ")")
  end
  join(io, properties_string, " ")
end
