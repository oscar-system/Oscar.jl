################################################
# 1: Constructors with toric variety as base
################################################

@doc raw"""
    function hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, p::MPolyRingElem; completeness_check::Bool = true)

Construct a hypersurface model, for which the user can specify a fiber ambient space
as well as divisor classes of the toric base space, in which the first two homogeneous
coordinates of the fiber ambient space transform.

# Examples
```jldoctest
julia> b = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> D1 = 2 * anticanonical_divisor_class(b)
Divisor class on a normal toric variety

julia> D2 = 3 * anticanonical_divisor_class(b)
Divisor class on a normal toric variety

julia> D3 = trivial_divisor_class(b)
Divisor class on a normal toric variety

julia> new_gens = string.(vcat(gens(cox_ring(b)), gens(cox_ring(fiber_ambient_space))))
6-element Vector{String}:
 "x1"
 "x2"
 "x3"
 "x"
 "y"
 "z"

julia> ambient_ring, (x1, x2, x3, x, y, z) = polynomial_ring(QQ, new_gens, cached=false)
(Multivariate polynomial ring in 6 variables over QQ, QQMPolyRingElem[x1, x2, x3, x, y, z])

julia> p = x^3 - y^2 + x1^12 * x * z^4 + x2^18 * z^6 + 13 * x3^3*x*y*z
x1^12*x*z^4 + x2^18*z^6 + 13*x3^3*x*y*z + x^3 - y^2

julia> h = hypersurface_model(b, fiber_ambient_space, [D1, D2, D3], p; completeness_check = false)
Hypersurface model over a concrete base
```
"""
function hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, p::MPolyRingElem; completeness_check::Bool = true)
  return hypersurface_model(base, fiber_ambient_space, fiber_twist_divisor_classes, string(p); completeness_check = completeness_check)
end

function hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, p::String; completeness_check::Bool = true)
  # Consistency checks
  gens_base_names = symbols(cox_ring(base))
  gens_fiber_names = symbols(cox_ring(fiber_ambient_space))
  if !isdisjoint(gens_base_names, gens_fiber_names)
    @vprint :FTheoryModelPrinter 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  if completeness_check
    @req is_complete(base) "Base space must be complete"
  end
  
  # Compute an ambient space
  ambient_space = _ambient_space(base, fiber_ambient_space, fiber_twist_divisor_classes)

  # Construct the model
  hypersurface_equation = eval_poly(p, cox_ring(ambient_space))
  @req is_homogeneous(hypersurface_equation) "Given hypersurface equation is not homogeneous"
  ds = [x.coeff for x in keys(homogeneous_components(hypersurface_equation))]
  @req length(ds) == 1 "Inconsistency in determining the degree of the hypersurface equation"
  @req ds[1] == divisor_class(anticanonical_divisor_class(ambient_space)).coeff "Degree of hypersurface equation differs from anticanonical bundle"
  explicit_model_sections = Dict{String, MPolyRingElem}()
  gens_S = gens(cox_ring(base))

  # The below code was removed because it is inconsistent with our standard use of explicit_model_sections. In particular,
  # explicit_model_sections should not list base coordinates as being named sections whose value is the base coordinate
  # For now, explicit_model_sections will be empty for hypersurface models. This will not be true ultimately, when we
  # allow for the definition of model sections for hypersurface models
  #
  # for k in 1:length(gens_S)
  #   explicit_model_sections[string(gens_S[k])] = gens_S[k]
  # end

  model = HypersurfaceModel(explicit_model_sections, hypersurface_equation, hypersurface_equation, base, ambient_space, fiber_ambient_space)
  set_attribute!(model, :partially_resolved, false)
  return model
end



################################################
# 2: Constructors with scheme as base
################################################

# Yet to come...


################################################
# 3: Constructors without specified base
################################################


@doc raw"""
    hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{Vector{Int64}}, p::MPolyRingElem)

This method constructs a hypersurface model over a base space that is not
fully specified. In the background, we construct a family of spaces to represent
the base space. This method requires the following information:
1. The names of the homogeneous coordinates of the coordinate ring of the generic member
of the family of bsae spaces.
2. The grading of the coordinate ring of the generic member of the family of base spaces.
3. The weights telling us how the fiber ambient space coordinates transform under the base.
4. The dimension of the generic member of the family of base spaces.
5. The fiber ambient space.
6. The hypersurface equation.

Note that many studies in the literature use the class of the anticanonical bundle
in their analysis. We anticipate this by adding this class as a variable of the
coordinate ring of the generic member of the family of base space, unless the user
already provides this grading. Our convention is that the first row of the grading
matrix refers to Kbar and that the homogeneous variable corresponding to this class
carries the name "Kbar".

The following example exemplifies this constructor.

# Examples
```jldoctest
julia> auxiliary_base_vars = ["a1", "a21", "a32", "a43", "a65", "w"]
6-element Vector{String}:
 "a1"
 "a21"
 "a32"
 "a43"
 "a65"
 "w"

julia> auxiliary_base_grading = [1 2 3 4 6 0; 0 -1 -2 -3 -5 1]
2×6 Matrix{Int64}:
 1   2   3   4   6  0
 0  -1  -2  -3  -5  1

julia> D1 = [4,0]
2-element Vector{Int64}:
 4
 0

julia> D2 = [6,0]
2-element Vector{Int64}:
 6
 0

julia> D3 = [0,0]
2-element Vector{Int64}:
 0
 0
 
julia> d = 3
3

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> auxiliary_ambient_ring, (a1, a21, a32, a43, a65, w, x, y, z)  = QQ[:a1, :a21, :a32, :a43, :a65, :w, :x, :y, :z]
(Multivariate polynomial ring in 9 variables over QQ, QQMPolyRingElem[a1, a21, a32, a43, a65, w, x, y, z])

julia> p = x^3 - y^2 - x * y * z * a1 + x^2 * z^2 * a21 * w - y * z^3 * a32 * w^2 + x * z^4 * a43 * w^3 + z^6 * a65 * w^5
-a1*x*y*z + a21*w*x^2*z^2 - a32*w^2*y*z^3 + a43*w^3*x*z^4 + a65*w^5*z^6 + x^3 - y^2

julia> h = hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, fiber_ambient_space, [D1, D2, D3], p)
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base
```
"""
function hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{Vector{Int64}}, p::MPolyRingElem)
  return hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, fiber_ambient_space, transpose(matrix(ZZ, fiber_twist_divisor_classes)), p)
end

function hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::ZZMatrix, p::MPolyRingElem)
  
  # Compute simple information
  gens_fiber_names = [string(g) for g in gens(cox_ring(fiber_ambient_space))]
  set_base_vars = Set(auxiliary_base_vars)
  set_fiber_vars = Set(gens_fiber_names)
  set_p_vars = Set([string(g) for g in gens(parent(p))])
  
  # Conduct simple consistency checks
  @req d > 0 "The dimension of the base space must be positive"
  @req isdisjoint(set_base_vars, set_fiber_vars) "Variable names duplicated between base and fiber coordinates."
  @req union(set_base_vars, set_fiber_vars) == set_p_vars "Variables names for polynomial p do not match variable choice for base and fiber"
  @req ncols(auxiliary_base_grading) == length(auxiliary_base_vars) "Number of base variables does not match the number of provided base gradings"
  
  # Inform about the assume Kbar grading
  @vprint :FTheoryModelPrinter 0 "Assuming that the first row of the given grading is the grading under Kbar\n\n"
  
  # Construct the spaces
  (S, auxiliary_base_space, auxiliary_ambient_space) = _construct_generic_sample(auxiliary_base_grading, auxiliary_base_vars, d, fiber_ambient_space, fiber_twist_divisor_classes)

  # Map p to coordinate ring of ambient space
  gens_S = gens(S)
  image_list = Vector{MPolyRingElem}()
  for g in gens(parent(p))
    index = findfirst(u -> string(u) == string(g), gens(S))
    push!(image_list, gens(S)[index])
  end
  hypersurface_equation = hom(parent(p), S, image_list)(p)

  # Remember the parametrization of the hypersurface equation
  hypersurface_equation_parametrization = eval_poly(string(p), coordinate_ring(auxiliary_ambient_space))

  # Remember the explicit model sections
  gens_R = gens(coordinate_ring(auxiliary_base_space))
  explicit_model_sections = Dict{String, MPolyRingElem}()
  for k in 1:length(gens_R)
    explicit_model_sections[string(gens_R[k])] = gens_R[k]
  end

  # Construct the model
  model = HypersurfaceModel(explicit_model_sections, hypersurface_equation_parametrization, hypersurface_equation, auxiliary_base_space, auxiliary_ambient_space, fiber_ambient_space)
  set_attribute!(model, :partially_resolved, false)
  return model
end


################################################
# 4: Display
################################################

# Detailed printing
function Base.show(io::IO, ::MIME"text/plain", h::HypersurfaceModel)
  io = pretty(io)
  properties_string = String[]
  if is_partially_resolved(h)
    push!(properties_string, "Partially resolved hypersurface model over a")
  else
    push!(properties_string, "Hypersurface model over a")
  end
  if is_base_space_fully_specified(h)
    push!(properties_string, "concrete base")
  else
    push!(properties_string, "not fully specified base")
  end
  join(io, properties_string, " ")
end

# Terse and one line printing
function Base.show(io::IO, h::HypersurfaceModel)
  print(io, "Hypersurface model")
end
