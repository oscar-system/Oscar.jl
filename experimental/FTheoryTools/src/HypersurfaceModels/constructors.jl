################################################
# 1: Constructors with toric variety as base
################################################

@doc raw"""
    hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, p::MPolyRingElem; completeness_check::Bool = true)

Constructs a hypersurface model where fiber coordinates transform as sections of line bundles over the base.

- If two `fiber_twist_divisor_classes` are provided, they specify the line bundle charges of the **first two** fiber coordinates. All remaining fiber coordinates are treated as base-trivial (i.e., constant over the base).
- Alternatively, a `fiber_twist_divisor_class` can be given **for each fiber coordinate**, allowing full control over how the entire fiber ambient space twists over the base.

The hypersurface equation `p` defines a section of the anti-canonical bundle of the total ambient space and may also be passed as a string for convenience.

For a full explanation of hypersurface models, see [Hypersurface Models](@ref "Hypersurface Models").

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

julia> h = hypersurface_model(b, fiber_ambient_space, [D1, D2], p; completeness_check = false)
Hypersurface model over a concrete base

julia> h2 = hypersurface_model(b, fiber_ambient_space, [D1, D2], string(p); completeness_check = false)
Hypersurface model over a concrete base
```
"""
function hypersurface_model(bs::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, p::MPolyRingElem; completeness_check::Bool = true)  
  return hypersurface_model(bs, fiber_ambient_space, fiber_twist_divisor_classes, string(p); completeness_check = completeness_check)
end

function hypersurface_model(bs::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, p::String; completeness_check::Bool = true)
  @req all(toric_variety(c) == bs for c in fiber_twist_divisor_classes) "All fiber divisor classes must belong to be base space"
  if length(fiber_twist_divisor_classes) == n_rays(fiber_ambient_space)
    return _build_hypersurface_model(bs, fiber_ambient_space, fiber_twist_divisor_classes, p; completeness_check = completeness_check)
  end
  @req length(fiber_twist_divisor_classes) == 2 "Exactly two fiber twist divisor classes must be provided"
  @req n_rays(fiber_ambient_space) >= 2 "The fiber ambient space must have at least 2 rays"
  new_fiber_twist_divisor_classes = vcat(fiber_twist_divisor_classes, [trivial_divisor_class(bs) for k in 1:(n_rays(fiber_ambient_space) - length(fiber_twist_divisor_classes))])
  return _build_hypersurface_model(bs, fiber_ambient_space, new_fiber_twist_divisor_classes, p; completeness_check = completeness_check)
end


@doc raw"""
    hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, indices::Vector{Int}, p::MPolyRingElem; completeness_check::Bool = true)

Constructs a hypersurface model where the two fiber coordinates at the given
`indices` transform as sections of the line bundles associated to the specified
base divisor classes.

The polynomial `p`, which encodes the hypersurface equation, can also be passed as a string.

For a full explanation of hypersurface models, see [Hypersurface Models](@ref "Hypersurface Models").

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

julia> h = hypersurface_model(b, fiber_ambient_space, [D1, D2], [1, 2], p; completeness_check = false)
Hypersurface model over a concrete base

julia> h2 = hypersurface_model(b, fiber_ambient_space, [D1, D2], [1, 2], string(p); completeness_check = false)
Hypersurface model over a concrete base
```
"""
function hypersurface_model(bs::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, indices::Vector{Int}, p::MPolyRingElem; completeness_check::Bool = true)
  return hypersurface_model(bs, fiber_ambient_space, fiber_twist_divisor_classes, indices, string(p); completeness_check = completeness_check)
end

function hypersurface_model(bs::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, indices::Vector{Int}, p::String; completeness_check::Bool = true)
  @req all(toric_variety(c) == bs for c in fiber_twist_divisor_classes) "All fiber divisor classes must belong to be base space"  
  @req length(fiber_twist_divisor_classes) == 2 "Exactly two fiber twist divisor classes must be provided"
  @req length(Set(indices)) == 2 "Exactly two distinct indices must be provided"
  @req 1 <= minimum(indices) "The indices must not be smaller than 1"
  @req n_rays(fiber_ambient_space) >= maximum(indices) "The indices must not exceed the number of rays of the fiber ambient space"
  @req n_rays(fiber_ambient_space) >= 2 "The fiber ambient space must have at least 2 rays"
  new_fiber_twist_divisor_classes = [trivial_divisor_class(bs) for k in 1:(n_rays(fiber_ambient_space))]
  new_fiber_twist_divisor_classes[indices[1]] = fiber_twist_divisor_classes[1]
  new_fiber_twist_divisor_classes[indices[2]] = fiber_twist_divisor_classes[2]
  return _build_hypersurface_model(bs, fiber_ambient_space, new_fiber_twist_divisor_classes, p; completeness_check = completeness_check)
end

function _build_hypersurface_model(bs::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, p::String; completeness_check::Bool = true)
  # Consistency checks
  gens_base_names = symbols(cox_ring(bs))
  gens_fiber_names = symbols(cox_ring(fiber_ambient_space))
  if !isdisjoint(gens_base_names, gens_fiber_names)
    @vprint :FTheoryModelPrinter 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  if completeness_check
    @req is_complete(bs) "Base space must be complete"
  end
  @req length(fiber_twist_divisor_classes) == n_rays(fiber_ambient_space) "Number of fiber twist divisor classes must match number of rays in fiber ambient space"

  # Compute an ambient space
  ambient_space = _ambient_space(bs, fiber_ambient_space, fiber_twist_divisor_classes)

  # Construct the model
  hypersurface_equation = eval_poly(p, cox_ring(ambient_space))
  @req is_homogeneous(hypersurface_equation) "Given hypersurface equation is not homogeneous"
  ds = [x.coeff for x in keys(homogeneous_components(hypersurface_equation))]
  @req length(ds) == 1 "Inconsistency in determining the degree of the hypersurface equation"
  @req ds[1] == divisor_class(anticanonical_divisor_class(ambient_space)).coeff "Degree of hypersurface equation differs from anticanonical bundle"
  explicit_model_sections = Dict{String, MPolyRingElem}()
  gens_S = gens(cox_ring(bs))

  # The below code was removed because it is inconsistent with our standard use of explicit_model_sections. In particular,
  # explicit_model_sections should not list base coordinates as being named sections whose value is the base coordinate
  # For now, explicit_model_sections will be empty for hypersurface models. This will not be true ultimately, when we
  # allow for the definition of model sections for hypersurface models
  #
  # for k in 1:length(gens_S)
  #   explicit_model_sections[string(gens_S[k])] = gens_S[k]
  # end

  model = HypersurfaceModel(explicit_model_sections, hypersurface_equation, hypersurface_equation, bs, ambient_space, fiber_ambient_space)
  set_attribute!(model, :partially_resolved, false)
  return model
end


################################################
# 2: Constructors with unspecified base
################################################

@doc raw"""
    hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{Vector{Int64}}, p::MPolyRingElem)

Constructs a hypersurface model over an unspecified base space by defining a **family** of base varieties
via auxiliary data.

The base is represented by:
- `auxiliary_base_vars`: Names of the homogeneous coordinates in the Cox ring of a generic base variety.
- `auxiliary_base_grading`: A grading matrix specifying the line bundle degrees of these coordinates.
- `d`: The (Krull) dimension of the base.

Additionally, the fiber data must be provided:
- `fiber_ambient_space`: A toric variety defining the ambient space for the fiber.
- **`fiber_twist_divisor_classes`** — A list of weight vectors indicating how fiber coordinates transform under the base line bundles.

    - If two vectors are provided, they define the line bundle charges of the **first two** fiber coordinates; all remaining fiber coordinates are assumed to be base-trivially fibered (i.e., constant over the base).
    - Alternatively, one weight vector can be provided **per fiber coordinate**, giving full control over how each fiber coordinate twists over the base.

- `p`: A multivariate polynomial defining the hypersurface equation in the total ambient space.

To ensure that the hypersurface is Calabi--Yau, the equation `p` must define a section of the anti-canonical bundle of the total
ambient space. If the provided grading matrix does not include this class, a new grading row corresponding to the anti-canonical
class is automatically added. By convention, this first row corresponds to ``\overline{K}_A`` and the associated variable is named `"Kbar"`.

For convenience, `fiber_twist_divisor_classes` can also be passed as a `ZZMatrix` instead of a `Vector{Vector{Int64}}`.

For a full explanation of hypersurface models, see [Hypersurface Models](@ref "Hypersurface Models").

# Examples
```jldoctest
julia> auxiliary_base_vars = ["a1", "a21", "a32", "a43", "a65", "w"];

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
 
julia> d = 3
3

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> auxiliary_ambient_ring, (a1, a21, a32, a43, a65, w, x, y, z)  = QQ[:a1, :a21, :a32, :a43, :a65, :w, :x, :y, :z]
(Multivariate polynomial ring in 9 variables over QQ, QQMPolyRingElem[a1, a21, a32, a43, a65, w, x, y, z])

julia> p = x^3 - y^2 - x * y * z * a1 + x^2 * z^2 * a21 * w - y * z^3 * a32 * w^2 + x * z^4 * a43 * w^3 + z^6 * a65 * w^5
-a1*x*y*z + a21*w*x^2*z^2 - a32*w^2*y*z^3 + a43*w^3*x*z^4 + a65*w^5*z^6 + x^3 - y^2

julia> h = hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, fiber_ambient_space, [D1, D2], p)
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base
```
"""
function hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{Vector{Int64}}, p::MPolyRingElem)
  return hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, fiber_ambient_space, transpose(matrix(ZZ, fiber_twist_divisor_classes)), p)
end

function hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::ZZMatrix, p::MPolyRingElem)
  @req nrows(fiber_twist_divisor_classes) == nrows(auxiliary_base_grading) "Twist does not match the number of base gradings"
  if ncols(fiber_twist_divisor_classes) == n_rays(fiber_ambient_space)
    return build_hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, fiber_ambient_space, fiber_twist_divisor_classes, p)  
  end
  @req ncols(fiber_twist_divisor_classes) == 2 "Exactly two fiber twist divisor classes must be specified"
  ext_mat = transpose(matrix(ZZ, [[0 for l in 1:nrows(fiber_twist_divisor_classes)] for k in 1:(n_rays(fiber_ambient_space) - ncols(fiber_twist_divisor_classes))]))
  new_fiber_twist_divisor_classes = hcat(fiber_twist_divisor_classes, ext_mat)
  return build_hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, fiber_ambient_space, new_fiber_twist_divisor_classes, p)
end


@doc raw"""
    hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{Vector{Int64}}, indices::Vector{Int}, p::MPolyRingElem)

Constructs a hypersurface model over an unspecified base space by defining a **family** of base varieties
via auxiliary data.

The base is represented by:
- `auxiliary_base_vars`: Names of the homogeneous coordinates in the Cox ring of a generic base variety.
- `auxiliary_base_grading`: A grading matrix specifying the line bundle degrees of these coordinates.
- `d`: The (Krull) dimension of the base.

Additionally, the fiber data must be provided:
- `fiber_ambient_space`: A toric variety defining the ambient space for the fiber.
- `fiber_twist_divisor_classes`: A list of two weight vectors specifying how the fiber coordinates at the given `indices` transform under the base line bundles.
- `p`: A multivariate polynomial defining the hypersurface equation in the total ambient space.

To ensure that the hypersurface is Calabi--Yau, the equation `p` must define a section of the anti-canonical bundle of the total
ambient space. If the provided grading matrix does not include this class, a new grading row corresponding to the anti-canonical
class is automatically added. By convention, this first row corresponds to ``\overline{K}_A`` and the associated variable is named `"Kbar"`.

For convenience, `fiber_twist_divisor_classes` can also be provided as `ZZMatrix`.

For a full explanation of hypersurface models, see [Hypersurface Models](@ref "Hypersurface Models").

# Examples
```jldoctest
julia> auxiliary_base_vars = ["a1", "a21", "a32", "a43", "a65", "w"];

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
 
julia> d = 3
3

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal toric variety

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> auxiliary_ambient_ring, (a1, a21, a32, a43, a65, w, x, y, z)  = QQ[:a1, :a21, :a32, :a43, :a65, :w, :x, :y, :z]
(Multivariate polynomial ring in 9 variables over QQ, QQMPolyRingElem[a1, a21, a32, a43, a65, w, x, y, z])

julia> p = x^3 - y^2 - x * y * z * a1 + x^2 * z^2 * a21 * w - y * z^3 * a32 * w^2 + x * z^4 * a43 * w^3 + z^6 * a65 * w^5
-a1*x*y*z + a21*w*x^2*z^2 - a32*w^2*y*z^3 + a43*w^3*x*z^4 + a65*w^5*z^6 + x^3 - y^2

julia> h = hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, fiber_ambient_space, [D1, D2], [1, 2], p)
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base
```
"""
function hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{Vector{Int64}}, indices::Vector{Int}, p::MPolyRingElem)
  return hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, fiber_ambient_space, transpose(matrix(ZZ, fiber_twist_divisor_classes)), indices, p)
end


function hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::ZZMatrix, indices::Vector{Int}, p::MPolyRingElem)
  @req ncols(fiber_twist_divisor_classes) == 2 "Exactly two fiber twist divisor classes must be specified"
  @req length(Set(indices)) == 2 "Exactly two distinct indices must be provided"
  @req 1 <= minimum(indices) "The indices must not be smaller than 1"
  @req n_rays(fiber_ambient_space) >= maximum(indices) "The indices must not exceed the number of rays of the fiber ambient space"
  @req nrows(fiber_twist_divisor_classes) == nrows(auxiliary_base_grading) "Twist does not match the number of base gradings"
  new_fiber_twist_divisor_classes = zero_matrix(ZZ, nrows(auxiliary_base_grading), n_rays(fiber_ambient_space))
  for l in 1:nrows(fiber_twist_divisor_classes)
    new_fiber_twist_divisor_classes[l, indices[1]] = fiber_twist_divisor_classes[l, 1]
    new_fiber_twist_divisor_classes[l, indices[2]] = fiber_twist_divisor_classes[l, 2]
  end
  return build_hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, fiber_ambient_space, new_fiber_twist_divisor_classes, p)
end


function build_hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::ZZMatrix, p::MPolyRingElem)
  
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
# 3: Display
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
