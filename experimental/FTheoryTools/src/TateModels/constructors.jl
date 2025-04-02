################################################
# 1: Constructors with toric variety as base
################################################

@doc raw"""
    global_tate_model(base::NormalToricVariety; completeness_check::Bool = true)

This method constructs a global Tate model over a given toric base
3-fold. The Tate sections ``a_i`` are taken with (pseudo) random coefficients.

# Examples
```jldoctest
julia> t = global_tate_model(sample_toric_variety(); completeness_check = false)
Global Tate model over a concrete base
```
"""
global_tate_model(base::NormalToricVariety; completeness_check::Bool = true) = global_tate_model(base, _tate_sections(base); completeness_check = completeness_check)


@doc raw"""
    global_tate_model(base::NormalToricVariety, ais::Vector{T}; completeness_check::Bool = true) where {T<:MPolyRingElem}

This method operates analogously to `global_tate_model(base::NormalToricVarietyType)`.
The only difference is that the Tate sections ``a_i`` can be specified with non-generic values.

# Examples
```jldoctest
julia> chosen_base = sample_toric_variety()
Normal toric variety

julia> a1 = generic_section(anticanonical_bundle(chosen_base));

julia> a2 = generic_section(anticanonical_bundle(chosen_base)^2);

julia> a3 = generic_section(anticanonical_bundle(chosen_base)^3);

julia> a4 = generic_section(anticanonical_bundle(chosen_base)^4);

julia> a6 = generic_section(anticanonical_bundle(chosen_base)^6);

julia> t = global_tate_model(chosen_base, [a1, a2, a3, a4, a6]; completeness_check = false)
Global Tate model over a concrete base
```
"""
function global_tate_model(base::NormalToricVariety, ais::Vector{T}; completeness_check::Bool = true) where {T<:MPolyRingElem}
  @req length(ais) == 5 "All the Tate sections a1, a2, a3, a4, a6 must be provided"
  return global_tate_model(base, Dict("a1" => ais[1], "a2" => ais[2], "a3" => ais[3], "a4" => ais[4], "a6" => ais[5]), Dict{String, MPolyRingElem}(); completeness_check = completeness_check)
end

function global_tate_model(base::NormalToricVariety,
                           explicit_model_sections::Dict{String, <: Union{MPolyRingElem, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}},
                           model_section_parametrization::Dict{String, MPolyRingElem};
                           completeness_check::Bool = true)
  vs = collect(values(explicit_model_sections))
  @req all(k -> parent(k) == cox_ring(base), vs) "All Tate sections must reside in the Cox ring of the base toric variety"
  @req haskey(explicit_model_sections, "a1") "Tate section a1 must be specified"
  @req haskey(explicit_model_sections, "a2") "Tate section a2 must be specified"
  @req haskey(explicit_model_sections, "a3") "Tate section a3 must be specified"
  @req haskey(explicit_model_sections, "a4") "Tate section a4 must be specified"
  @req haskey(explicit_model_sections, "a6") "Tate section a6 must be specified"
  vs2 = collect(keys(model_section_parametrization))
  @req all(in(["a1", "a2", "a3", "a4", "a6"]), vs2) "Only the Tate sections a1, a2, a3, a4, a6 must be parametrized"
  
  gens_base_names = symbols(cox_ring(base))
  if (:x in gens_base_names) || (:y in gens_base_names) || (:z in gens_base_names)
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
  ais = [explicit_model_sections["a1"], explicit_model_sections["a2"], explicit_model_sections["a3"], explicit_model_sections["a4"], explicit_model_sections["a6"]]
  pt = _tate_polynomial(ais, cox_ring(ambient_space))
  model = GlobalTateModel(explicit_model_sections, model_section_parametrization, pt, base, ambient_space)
  set_attribute!(model, :partially_resolved, false)
  return model
end


################################################
# 2: Constructors with scheme as base
################################################

# Yet to come...
# This requires that the ai are stored as sections of the anticanonical bundle, and not "just" polynomials.
# -> Types to be generalized then.


################################################
# 3: Constructors without specified base
################################################

@doc raw"""
    global_tate_model(auxiliary_base_ring::MPolyRing, auxiliary_base_grading::Matrix{Int64}, d::Int, ais::Vector{T}) where {T<:MPolyRingElem}

This method constructs a global Tate model over a base space that is not
fully specified.

Note that many studies in the literature use the class of the anticanonical bundle
in their analysis. We anticipate this by adding this class as a variable of the
auxiliary base space, unless the user already provides this grading. Our convention
is that the first grading refers to Kbar and that the homogeneous variable corresponding
to this class carries the name "Kbar".

The following code exemplifies this approach.

# Examples
```jldoctest
julia> auxiliary_base_ring, (a10, a21, a32, a43, a65, w) = QQ[:a10, :a21, :a32, :a43, :a65, :w];

julia> auxiliary_base_grading = [1 2 3 4 6 0; 0 -1 -2 -3 -5 1]
2×6 Matrix{Int64}:
 1   2   3   4   6  0
 0  -1  -2  -3  -5  1

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = global_tate_model(auxiliary_base_ring, auxiliary_base_grading, 3, ais)
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base
```
"""
function global_tate_model(auxiliary_base_ring::MPolyRing, auxiliary_base_grading::Matrix{Int64}, d::Int, ais::Vector{T}) where {T<:MPolyRingElem}
  
  # Execute consistency checks
  gens_base_names = [string(g) for g in symbols(auxiliary_base_ring)]
  @req length(ais) == 5 "We expect exactly 5 Tate sections"
  @req all(k -> parent(k) == auxiliary_base_ring, ais) "All Tate sections must reside in the provided auxiliary base ring"
  @req d > 0 "The dimension of the base space must be positive"
  if ("x" in gens_base_names) || ("y" in gens_base_names) || ("z" in gens_base_names)
    @vprint :FTheoryModelPrinter 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  
  # inform about the assume Kbar grading
  @vprint :FTheoryModelPrinter 0 "Assuming that the first row of the given grading is the grading under Kbar\n\n"
  
  # Compute the Tate polynomial
  (S, auxiliary_base_space, auxiliary_ambient_space) = _construct_generic_sample(auxiliary_base_grading, gens_base_names, d)
  ring_map = hom(parent(ais[1]), S, gens(S)[1:ngens(parent(ais[1]))])
  (a1, a2, a3, a4, a6) = [ring_map(k) for k in ais]
  pt = _tate_polynomial([a1, a2, a3, a4, a6], coordinate_ring(auxiliary_ambient_space))

  # Compute explicit model sections
  explicit_model_sections = Dict("a1" => a1, "a2" => a2, "a3" => a3, "a4" => a4, "a6" => a6)
  section_candidates = gens(S)
  for k in section_candidates
    haskey(explicit_model_sections, string(k)) || (explicit_model_sections[string(k)] = k)
  end

  # Compute model_section_parametrization
  model_section_parametrization = Dict{String, MPolyRingElem}()
  vars_S = [string(k) for k in symbols(S)]
  if !("a1" in vars_S) || (a1 != eval_poly("a1", parent(a1)))
    model_section_parametrization["a1"] = a1
  end
  if !("a2" in vars_S) || (a2 != eval_poly("a2", parent(a2)))
    model_section_parametrization["a2"] = a2
  end
  if !("a3" in vars_S) || (a3 != eval_poly("a3", parent(a3)))
    model_section_parametrization["a3"] = a3
  end
  if !("a4" in vars_S) || (a4 != eval_poly("a4", parent(a4)))
    model_section_parametrization["a4"] = a4
  end
  if !("a6" in vars_S) || (a6 != eval_poly("a6", parent(a6)))
    model_section_parametrization["a6"] = a6
  end

  # Compute model and return it
  model = GlobalTateModel(explicit_model_sections, model_section_parametrization, pt, auxiliary_base_space, auxiliary_ambient_space)
  set_attribute!(model, :partially_resolved, false)
  return model
end



################################################
# 4: Display
################################################

function Base.show(io::IO, t::GlobalTateModel)
  properties_string = String[]
  if is_partially_resolved(t)
    push!(properties_string, "Partially resolved global Tate model over a")
  else
    push!(properties_string, "Global Tate model over a")
  end
  if is_base_space_fully_specified(t)
    push!(properties_string, "concrete base")
  else
    push!(properties_string, "not fully specified base")
  end
  if has_model_description(t)
    push!(properties_string, "-- " * model_description(t))
    if has_model_parameters(t)
      push!(properties_string, "with parameter values (" * join(["$key = $(string(val))" for (key, val) in model_parameters(t)], ", ") * ")")
    end
  end
  if has_arxiv_id(t)
    push!(properties_string, "based on arXiv paper " * arxiv_id(t))
  end
  if has_arxiv_model_equation_number(t)
    push!(properties_string, "Eq. (" * arxiv_model_equation_number(t) * ")")
  end
  join(io, properties_string, " ")
end
