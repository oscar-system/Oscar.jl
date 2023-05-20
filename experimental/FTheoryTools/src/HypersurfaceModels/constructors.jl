################################################
# 1: Constructors with toric variety as base
################################################

@doc raw"""
    hypersurface_model(base::AbstractNormalToricVariety; completeness_check::Bool = true)

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
function hypersurface_model(base::AbstractNormalToricVariety; completeness_check::Bool = true)
  fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
  set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])
  D1 = 2 * anticanonical_divisor_class(base)
  D2 = 3 * anticanonical_divisor_class(base)
  return hypersurface_model(base, fiber_ambient_space, D1, D2; completeness_check = completeness_check)
end


@doc raw"""
    hypersurface_model(base::AbstractNormalToricVariety, fiber_ambient_space::AbstractNormalToricVariety, D1::ToricDivisorClass, D2::ToricDivisorClass; completeness_check::Bool = true)

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
function hypersurface_model(base::AbstractNormalToricVariety, fiber_ambient_space::AbstractNormalToricVariety, D1::ToricDivisorClass, D2::ToricDivisorClass; completeness_check::Bool = true)
  
  # Consistency checks
  gens_base_names = [string(g) for g in gens(cox_ring(base))]
  gens_fiber_names = [string(g) for g in gens(cox_ring(fiber_ambient_space))]
  if length(findall(in(gens_base_names), gens_fiber_names)) != 0
    @vprint :HypersurfaceModel 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  if completeness_check
    @req is_complete(base) "Base space must be complete"
  end
  
  # Compute an ambient space
  ambient_space = _ambient_space(base, fiber_ambient_space, D1, D2)
  
  # Construct the model
  hypersurface_equation = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(ambient_space))])
  model = HypersurfaceModel(toric_covered_scheme(base), toric_covered_scheme(ambient_space), toric_covered_scheme(fiber_ambient_space), hypersurface_equation)
  set_attribute!(model, :base_fully_specified, true)
  return model
end


################################################
# 2: Constructors with toric scheme as base
################################################

@doc raw"""
    hypersurface_model(base::ToricCoveredScheme; completeness_check::Bool = true)

Construct a hypersurface model. This constructor takes $\mathbb{P}^{2,3,1}$ as fiber
ambient space with coordinates $[x:y:z]$ and ensures that $x$ transforms as
$2 \overline{K}_{B_3}$ and $y$ as $3 \overline{K}_{B_3}$.

# Examples
```jldoctest
julia> base = projective_space(ToricCoveredScheme, 2)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0], [0, 1], [-1, -1]]

julia> hypersurface_model(base; completeness_check = false)
Hypersurface model over a concrete base
```
"""
hypersurface_model(base::ToricCoveredScheme; completeness_check::Bool = true) = hypersurface_model(underlying_toric_variety(base); completeness_check = completeness_check)


@doc raw"""
    hypersurface_model(base::ToricCoveredScheme, fiber_ambient_space::ToricCoveredScheme, D1::ToricDivisorClass, D2::ToricDivisorClass; completeness_check::Bool = true)

Construct a hypersurface model, for which the user can specify a fiber ambient space
as well as divisor classes of the toric base space, in which the first two homogeneous
coordinates of the fiber ambient space transform.

# Examples
```jldoctest
julia> base = projective_space(ToricCoveredScheme, 2)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0], [0, 1], [-1, -1]]

julia> fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
Normal, non-affine, simplicial, projective, 2-dimensional toric variety without torusfactor

julia> set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])

julia> fiber_ambient_space = ToricCoveredScheme(fiber_ambient_space)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[-1, 1//3], [1, -1//2], [0, 1]]

julia> D1 = 2 * anticanonical_divisor_class(underlying_toric_variety(base))
Divisor class on a normal toric variety

julia> D2 = 3 * anticanonical_divisor_class(underlying_toric_variety(base))
Divisor class on a normal toric variety

julia> h = hypersurface_model(base, fiber_ambient_space, D1, D2; completeness_check = false)
Hypersurface model over a concrete base
```
"""
hypersurface_model(base::ToricCoveredScheme, fiber_ambient_space::ToricCoveredScheme, D1::ToricDivisorClass, D2::ToricDivisorClass; completeness_check::Bool = true) = hypersurface_model(underlying_toric_variety(base), underlying_toric_variety(fiber_ambient_space), D1, D2; completeness_check = completeness_check)


################################################
# 3: Constructors with scheme as base
################################################

# Yet to come...


################################################
# 4: Constructors without specified base
################################################


@doc raw"""
    hypersurface_model(p::MPolyRingElem, auxiliary_base_ring::MPolyRing, d::Int)

This method constructs a hypersurface model over a base space that is not
fully specified. The following example exemplifies this approach.

# Examples
```jldoctest
julia> auxiliary_base_ring, (a10, a21, a32, a43, a65, w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = global_tate_model(ais, auxiliary_base_ring, 3)
Global Tate model over a not fully specified base
```
"""
function hypersurface_model(p::MPolyRingElem, auxiliary_base_ring::MPolyRing, d::Int)
  # Do something cool!
  #=
  @req length(ais) == 5 "We expect exactly 5 Tate sections"
  @req all(k -> parent(k) == auxiliary_base_ring, ais) "All Tate sections must reside in the provided auxiliary base ring"
  @req d > 0 "The dimension of the base space must be positive"
  @req ngens(auxiliary_base_ring) >= d "We expect at least as many base variables as the desired base dimension"
  gens_base_names = [string(g) for g in gens(auxiliary_base_ring)]
  if ("x" in gens_base_names) || ("y" in gens_base_names) || ("z" in gens_base_names)
    @vprint :GlobalTateModel 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  
  # convert Tate sections into polynomials of the auxiliary base
  auxiliary_base_space = _auxiliary_base_space([string(k) for k in gens(auxiliary_base_ring)], d)
  S = cox_ring(auxiliary_base_space)
  ring_map = hom(auxiliary_base_ring, S, gens(S))
  (a1, a2, a3, a4, a6) = [ring_map(k) for k in ais]
  
  # construct model
  auxiliary_ambient_space = _ambient_space_from_base(auxiliary_base_space)
  pt = _tate_polynomial([a1, a2, a3, a4, a6], cox_ring(auxiliary_ambient_space))
  model = GlobalTateModel(a1, a2, a3, a4, a6, pt, toric_covered_scheme(auxiliary_base_space), toric_covered_scheme(auxiliary_ambient_space))
  set_attribute!(model, :base_fully_specified, false)
  return model=#
end


################################################
# 5: Display
################################################

function Base.show(io::IO, h::HypersurfaceModel)
  if base_fully_specified(h)
    print(io, "Hypersurface model over a concrete base")
  else
    print(io, "Hypersurface model over a not fully specified base")
  end
end
