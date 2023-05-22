################################################
# 1: Constructors with toric variety as base
################################################

@doc raw"""
    global_tate_model(base::NormalToricVarietyType; completeness_check::Bool = true)

This method constructs a global Tate model over a given toric base
3-fold. The Tate sections ``a_i`` are taken with (pseudo) random coefficients.

# Examples
```jldoctest
julia> t = global_tate_model(sample_toric_variety(); completeness_check = false)
Global Tate model over a concrete base
```
"""
global_tate_model(base::NormalToricVarietyType; completeness_check::Bool = true) = global_tate_model(base, _tate_sections(base); completeness_check = completeness_check)


@doc raw"""
    global_tate_model(base::NormalToricVarietyType, ais::Vector{T}; completeness_check::Bool = true) where {T<:MPolyRingElem}

This method operates analogously to `global_tate_model(base::NormalToricVarietyType)`.
The only difference is that the Tate sections ``a_i`` can be specified with non-generic values.

# Examples
```jldoctest
julia> base = sample_toric_variety()
Normal toric variety

julia> a1 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base))]);

julia> a2 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^2)]);

julia> a3 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^3)]);

julia> a4 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)]);

julia> a6 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)]);

julia> t = global_tate_model(base, [a1, a2, a3, a4, a6]; completeness_check = false)
Global Tate model over a concrete base
```
"""
function global_tate_model(base::NormalToricVarietyType, ais::Vector{T}; completeness_check::Bool = true) where {T<:MPolyRingElem}
  @req length(ais) == 5 "We require exactly 5 Tate sections"
  @req all(k -> parent(k) == cox_ring(base), ais) "All Tate sections must reside in the Cox ring of the base toric variety"
  
  gens_base_names = [string(g) for g in gens(cox_ring(base))]
  if ("x" in gens_base_names) || ("y" in gens_base_names) || ("z" in gens_base_names)
    @vprint :GlobalTateModel 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  
  if completeness_check
    @req is_complete(base) "Base space must be complete"
  end
  
  # construct the ambient space
  fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
  set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])
  D1 = 2 * anticanonical_divisor_class(base)
  D2 = 3 * anticanonical_divisor_class(base)
  ambient_space = _ambient_space(base, fiber_ambient_space, D1, D2)
  
  # construct the model
  pt = _tate_polynomial(ais, cox_ring(ambient_space))
  model = GlobalTateModel(ais[1], ais[2], ais[3], ais[4], ais[5], pt, toric_covered_scheme(base), toric_covered_scheme(ambient_space))
  set_attribute!(model, :base_fully_specified, true)
  return model
end


################################################
# 2: Constructors with toric scheme as base
################################################


@doc raw"""
    global_tate_model(base::ToricCoveredScheme; completeness_check::Bool = true)

This method constructs a global Tate model over a given toric scheme base
3-fold. The Tate sections ``a_i`` are taken with (pseudo) random coefficients.

# Examples
```jldoctest
julia> t = global_tate_model(sample_toric_scheme(); completeness_check = false)
Global Tate model over a concrete base
```
"""
global_tate_model(base::ToricCoveredScheme; completeness_check::Bool = true) = global_tate_model(underlying_toric_variety(base), completeness_check = completeness_check)


@doc raw"""
    global_tate_model(base::ToricCoveredScheme, ais::Vector{T}; completeness_check::Bool = true) where {T<:MPolyRingElem}

This method operates analogously to `global_tate_model(base::ToricCoveredScheme)`.
The only difference is that the Tate sections ``a_i`` can be specified with non-generic values.

# Examples
```jldoctest
julia> base = sample_toric_scheme()
Scheme of a toric variety

julia> a1 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base)))]);

julia> a2 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base))^2)]);

julia> a3 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base))^3)]);

julia> a4 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base))^4)]);

julia> a6 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base))^6)]);

julia> t = global_tate_model(base, [a1, a2, a3, a4, a6]; completeness_check = false)
Global Tate model over a concrete base
```
"""
global_tate_model(base::ToricCoveredScheme, ais::Vector{T}; completeness_check::Bool = true) where {T<:MPolyRingElem} = global_tate_model(underlying_toric_variety(base), ais; completeness_check = completeness_check)


################################################
# 3: Constructors with scheme as base
################################################

# Yet to come...
# This requires that the ai are stored as sections of the anticanonical bundle, and not "just" polynomials.
# -> Types to be generalized then.


################################################
# 4: Constructors without specified base
################################################

@doc raw"""
    global_tate_model(auxiliary_base_ring::MPolyRing, auxiliary_base_grading::Matrix{Int64}, d::Int, ais::Vector{T}; toric_sample = true) where {T<:MPolyRingElem}

This method constructs a global Tate model over a base space that is not
fully specified.

Note that many studies in the literature use the class of the anticanonical bundle
in their analysis. We anticipate this by adding this class as a variable of the
auxiliary base space, unless the user already provides this grading. Our convention
is that the first grading refers to Kbar and that the homogeneous variable corresponding
to this class carries the name "Kbar".

The following example exemplifies this approach.

# Examples
```jldoctest
julia> auxiliary_base_ring, (a10, a21, a32, a43, a65, w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

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

julia> t = global_tate_model(auxiliary_base_ring, auxiliary_base_grading, 3, ais; toric_sample = false)
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base
```
"""
function global_tate_model(auxiliary_base_ring::MPolyRing, auxiliary_base_grading::Matrix{Int64}, d::Int, ais::Vector{T}; toric_sample = true) where {T<:MPolyRingElem}
  
  # Is there a grading [1, 0, ..., 0]?
  Kbar_grading_present = false
  for i in 1:ncols(auxiliary_base_grading)
    col = auxiliary_base_grading[:,i]
    if length(col) == 1 && col[1] == 1
      Kbar_grading_present = true
      break;
    end
    if Set(col[2:length(col)]) == Set([0]) && col[1] == 1
      Kbar_grading_present = true
      break;
    end
  end
  
  # If Kbar is not present, extend the auxiliary_base_vars accordingly as well as the grading
  auxiliary_base_vars = gens(auxiliary_base_ring)
  gens_base_names = [string(g) for g in auxiliary_base_vars]
  if Kbar_grading_present == false
    @req ("Kbar" in gens_base_names) == false "Variable Kbar used as base variable, but grading of Kbar not introduced."
    Kbar_grading = [0 for i in 1:nrows(auxiliary_base_grading)]
    Kbar_grading[1] = 1
    auxiliary_base_grading = hcat(auxiliary_base_grading, Kbar_grading)
    push!(gens_base_names, "Kbar")
  end
  
  # Execute consistency checks
  @req length(ais) == 5 "We expect exactly 5 Tate sections"
  @req all(k -> parent(k) == auxiliary_base_ring, ais) "All Tate sections must reside in the provided auxiliary base ring"
  @req d > 0 "The dimension of the base space must be positive"
  if d == 1
    @req length(gens_base_names) - nrows(auxiliary_base_grading) > d "We expect a number of base variables that is strictly greater than one plus the number of scaling relations"
  else
    @req length(gens_base_names) - nrows(auxiliary_base_grading) >= d "We expect at least as many base variables as the sum of the desired base dimension and the number of scaling relations"
  end
  if ("x" in gens_base_names) || ("y" in gens_base_names) || ("z" in gens_base_names)
    @vprint :GlobalTateModel 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  
  # inform about the assume Kbar grading
  @vprint :FTheoryConstructorInformation 0 "Assuming that the first row of the given grading is the grading under Kbar\n\n"
  
  # Construct the model
  if toric_sample
    (S, auxiliary_base_space, auxiliary_ambient_space) = _construct_toric_sample(auxiliary_base_grading, gens_base_names, d)
    R = cox_ring(auxiliary_ambient_space)
  else
    (S, auxiliary_base_space, auxiliary_ambient_space) = _construct_generic_sample(auxiliary_base_grading, gens_base_names, d)
    R = coordinate_ring(auxiliary_ambient_space)
  end
  ring_map = hom(parent(ais[1]), S, gens(S)[1:ngens(parent(ais[1]))])
  (a1, a2, a3, a4, a6) = [ring_map(k) for k in ais]
  pt = _tate_polynomial([a1, a2, a3, a4, a6], R)
  model = GlobalTateModel(a1, a2, a3, a4, a6, pt, auxiliary_base_space, auxiliary_ambient_space)
  set_attribute!(model, :base_fully_specified, false)
  return model
end



################################################
# 5: Display
################################################

function Base.show(io::IO, t::GlobalTateModel)
  properties_string = ["Global Tate model over a"]
  if base_fully_specified(t)
    push!(properties_string, "concrete base")
  else
    push!(properties_string, "not fully specified base")
  end
  if has_model_description(t)
    push!(properties_string, "-- " * string(get_attribute(t, :model_description)))
    if has_model_parameters(t)
      push!(properties_string, "with parameter values (" * join(["$key = $(string(val))" for (key, val) in model_parameters(t)], ", ") * ")")
    end
  end
  if has_arxiv_id(t)
    push!(properties_string, "based on arXiv paper " * string(get_attribute(t, :arxiv_id)))
  end
  if has_arxiv_model_equation_number(t)
    push!(properties_string, "Eq. (" * string(get_attribute(t, :arxiv_model_equation_number)) * ")")
  end
  join(io, properties_string, " ")
end
