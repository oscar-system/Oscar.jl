################################################
# 1: Constructors with toric variety as base
################################################

@doc raw"""
    hypersurface_model(base::NormalToricVariety; completeness_check::Bool = true)

Construct a hypersurface model. This constructor takes $\mathbb{P}^{2,3,1}$ as fiber
ambient space with coordinates $[x:y:z]$ and ensures that $x$ transforms as
$2 \overline{K}_{B_3}$ and $y$ as $3 \overline{K}_{B_3}$.

# Examples
```jldoctest
julia> base = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> hypersurface_model(base; completeness_check = false)
Hypersurface model over a concrete base
```
"""
function hypersurface_model(base::NormalToricVariety; completeness_check::Bool = true)
  fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
  set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])
  D1 = 2 * anticanonical_divisor_class(base)
  D2 = 3 * anticanonical_divisor_class(base)
  return hypersurface_model(base, fiber_ambient_space, D1, D2; completeness_check = completeness_check)
end


@doc raw"""
    hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, D1::ToricDivisorClass, D2::ToricDivisorClass; completeness_check::Bool = true)

Construct a hypersurface model, for which the user can specify a fiber ambient space
as well as divisor classes of the toric base space, in which the first two homogeneous
coordinates of the fiber ambient space transform.

# Examples
```jldoctest
julia> base = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal, non-affine, simplicial, projective, 2-dimensional toric variety without torusfactor

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> D1 = 2 * anticanonical_divisor_class(base)
Divisor class on a normal toric variety

julia> D2 = 3 * anticanonical_divisor_class(base)
Divisor class on a normal toric variety

julia> hypersurface_model(base, fiber_ambient_space, D1, D2; completeness_check = false)
Hypersurface model over a concrete base
```
"""
function hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, D1::ToricDivisorClass, D2::ToricDivisorClass; completeness_check::Bool = true)
  
  # Consistency checks
  gens_base_names = [string(g) for g in gens(cox_ring(base))]
  gens_fiber_names = [string(g) for g in gens(cox_ring(fiber_ambient_space))]
  if intersect(Set(gens_base_names), Set(gens_fiber_names)) != Set()
    @vprint :HypersurfaceModel 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  if completeness_check
    @req is_complete(base) "Base space must be complete"
  end
  
  # Compute an ambient space
  ambient_space = _ambient_space(base, fiber_ambient_space, D1, D2)
  
  # Construct the model
  hypersurface_equation = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(ambient_space))])
  model = HypersurfaceModel(base, ambient_space, fiber_ambient_space, hypersurface_equation)
  set_attribute!(model, :base_fully_specified, true)
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
    hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, D1::Vector{Int64}, D2::Vector{Int64}, p::MPolyRingElem; toric_sample = true)

This method constructs a hypersurface model over a base space that is not
fully specified. In the background, we construct an auxiliary toric base space.
This method requires the following information:
1. The names of the homogeneous coordinates of the auxiliary toric base space.
2. The grading of the Cox ring of the auxiliary toric base space.
3. The weights corresponding to the divisor class `D_1` of the auxiliary toric base space under which the first fiber coordinate transforms.
4. The weights corresponding to the divisor class `D_2` of the auxiliary toric base space under which the first fiber coordinate transforms.
5. The dimension of the auxiliary toric base space.
6. The fiber ambient space.
7. The hypersurface equation.

Note that many studies in the literature use the class of the anticanonical bundle
in their analysis. We anticipate this by adding this class as a variable of the
auxiliary base space, unless the user already provides this grading. Our convention
is that the first grading refers to Kbar and that the homogeneous variable corresponding
to this class carries the name "Kbar".

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

julia> d = 3
3

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal, non-affine, simplicial, projective, 2-dimensional toric variety without torusfactor

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> auxiliary_ambient_ring, (a1, a21, a32, a43, a65, w, x, y, z)  = QQ["a1", "a21", "a32", "a43", "a65", "w", "x", "y", "z"]
(Multivariate polynomial ring in 9 variables over QQ, QQMPolyRingElem[a1, a21, a32, a43, a65, w, x, y, z])

julia> p = x^3 - y^2 - x * y * z * a1 + x^2 * z^2 * a21 * w - y * z^3 * a32 * w^2 + x * z^4 * a43 * w^3 + z^6 * a65 * w^5
-a1*x*y*z + a21*w*x^2*z^2 - a32*w^2*y*z^3 + a43*w^3*x*z^4 + a65*w^5*z^6 + x^3 - y^2

julia> h = hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, fiber_ambient_space, D1, D2, p)
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base

julia> h = hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, d, fiber_ambient_space, D1, D2, p; toric_sample = false)
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base
```
"""
function hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, D1::Vector{Int64}, D2::Vector{Int64}, p::MPolyRingElem; toric_sample = true)
  
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
  if Kbar_grading_present == false
    @req ("Kbar" in auxiliary_base_vars) == false "Variable Kbar used as base variable, but grading of Kbar not introduced."
    Kbar_grading = [0 for i in 1:nrows(auxiliary_base_grading)]
    Kbar_grading[1] = 1
    auxiliary_base_grading = hcat(auxiliary_base_grading, Kbar_grading)
    push!(auxiliary_base_vars, "Kbar")
  end
  
  # Compute simple information
  gens_fiber_names = [string(g) for g in gens(cox_ring(fiber_ambient_space))]
  set_base_vars = Set(auxiliary_base_vars)
  set_fiber_vars = Set(gens_fiber_names)
  set_p_vars = Set([string(g) for g in gens(parent(p))])
  
  # Conduct simple consistency checks
  @req d > 0 "The dimension of the base space must be positive"
  if d == 1
    @req length(auxiliary_base_vars) - nrows(auxiliary_base_grading) > d "We expect a number of base variables that is strictly greater than one plus the number of scaling relations"
  else
    @req length(auxiliary_base_vars) - nrows(auxiliary_base_grading) >= d "We expect at least as many base variables as the sum of the desired base dimension and the number of scaling relations"
  end
  @req intersect(set_base_vars, set_fiber_vars) == Set() "Variable names duplicated between base and fiber coordinates."
  @req union(set_base_vars, set_fiber_vars) == set_p_vars "Variables names for polynomial p do not match variable choice for base and fiber"
  @req ncols(auxiliary_base_grading) == length(auxiliary_base_vars) "Number of base variables does not match the number of provided base gradings"
  
  # Inform about the assume Kbar grading
  @vprint :FTheoryConstructorInformation 0 "Assuming that the first row of the given grading is the grading under Kbar\n\n"
  
  # Construct the spaces
  if toric_sample
    (S, auxiliary_base_space, auxiliary_ambient_space) = _construct_toric_sample(auxiliary_base_grading, auxiliary_base_vars, d, fiber_ambient_space, D1, D2, p)
  else
    (S, auxiliary_base_space, auxiliary_ambient_space) = _construct_generic_sample(auxiliary_base_grading, auxiliary_base_vars, d, fiber_ambient_space, D1, D2, p)
  end
  
  # Map p to coordinate ring of ambient space
  gens_S = gens(S)
  image_list = Vector{MPolyRingElem}()
  for g in gens(parent(p))
    index = findfirst(u -> string(u) == string(g), gens(S))
    push!(image_list, gens(S)[index])
  end
  hypersurface_equation = hom(parent(p), S, image_list)(p)

  # Construct the model
  model = HypersurfaceModel(auxiliary_base_space, auxiliary_ambient_space, fiber_ambient_space, hypersurface_equation)
  set_attribute!(model, :base_fully_specified, false)
  return model
end


################################################
# 4: Display
################################################

function Base.show(io::IO, h::HypersurfaceModel)
  if base_fully_specified(h)
    print(io, "Hypersurface model over a concrete base")
  else
    print(io, "Hypersurface model over a not fully specified base")
  end
end
