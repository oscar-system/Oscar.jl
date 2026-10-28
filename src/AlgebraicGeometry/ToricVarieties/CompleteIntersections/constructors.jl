# Toric hypersurfaces and chosen complete-intersection presentations.

# A chosen regular-sequence presentation, not an arbitrary Cox ideal.
abstract type AbsCompleteIntersectionInToricVariety <: AbsCoveredScheme{QQField} end

@attributes mutable struct CompleteIntersectionInToricVariety <:
                           AbsCompleteIntersectionInToricVariety
  ambient::NormalToricVariety
  equations::Vector{MPolyDecRingElem{QQFieldElem,QQMPolyRingElem}}

  function CompleteIntersectionInToricVariety(
    X::NormalToricVariety, equations::AbstractVector;
    check_codimension::Bool=true,
  )
    fan = polyhedral_fan(X)
    @req is_pure(fan) && is_fulldimensional(fan) "The ambient fan must have only full-dimensional maximal cones"
    @req 1 <= length(equations) <= dim(X) "Provide between one and $(dim(X)) equations"
    R = cox_ring(X)
    @req all(f -> f isa MPolyDecRingElem && parent(f) === R, equations) "Equations must belong to the ambient Cox ring"
    @req all(f -> !iszero(f), equations) "Defining equations must be nonzero"
    @req all(is_homogeneous, equations) "Defining equations must be homogeneous in the Cox grading"

    # Copy the input polynomials to keep the stored presentation independent
    # of subsequent mutations of the caller's equations.
    Y = new(X, elem_type(R)[copy(f) for f in equations])
    _check_complete_intersection_presentation(Y; check_codimension)
    return Y
  end
end

@attributes mutable struct HypersurfaceInToricVariety <:
                           AbsCompleteIntersectionInToricVariety
  presentation::CompleteIntersectionInToricVariety

  function HypersurfaceInToricVariety(presentation::CompleteIntersectionInToricVariety)
    @req length(defining_equations(presentation)) == 1 "A hypersurface requires one equation"
    return new(presentation)
  end
end

@doc raw"""
    complete_intersection_in_toric_variety(X::NormalToricVariety, equations::AbstractVector;
                                           check_codimension=true)

Construct a complete-intersection presentation in a normal toric variety
from homogeneous Cox-ring equations. The provided input equations define the
presentation. Different equation lists may define the same closed subscheme.
The ambient fan must have only full-dimensional maximal cones; this is checked
even when `check_codimension=false`.

By default, check that the resulting scheme is empty or has the expected
codimension. Set `check_codimension=false` to skip this potentially expensive
check; this does not certify that the equations form a regular sequence.
For nonempty schemes on singular ambient varieties, certification also
requires Cartier defining divisor classes. Neither smoothness nor reducedness
is required. A nonzero constant equation defines an empty zero locus.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 2);

julia> x, y, z = gens(cox_ring(X));

julia> C = complete_intersection_in_toric_variety(X, [x, y]);

julia> C isa CompleteIntersectionInToricVariety
true
```
"""
function complete_intersection_in_toric_variety(
  X::NormalToricVariety, equations::AbstractVector;
  check_codimension::Bool=true,
)
  return CompleteIntersectionInToricVariety(X, equations; check_codimension)
end

@doc raw"""
    hypersurface_in_toric_variety(X::NormalToricVariety, equation::MPolyDecRingElem;
                                  check_codimension=true)

Construct a hypersurface in a normal toric variety from one homogeneous
Cox-ring equation. The ambient fan must have only full-dimensional maximal
cones. Singular and nonreduced hypersurfaces are supported.
By default the hypersurface is checked for the expected codimension;
`check_codimension=false` skips this potentially expensive check.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 2);

julia> x = first(gens(cox_ring(X)));

julia> H = hypersurface_in_toric_variety(X, x);

julia> H isa HypersurfaceInToricVariety
true
```
"""
function hypersurface_in_toric_variety(
  X::NormalToricVariety, equation::MPolyDecRingElem;
  check_codimension::Bool=true,
)
  return HypersurfaceInToricVariety(
    complete_intersection_in_toric_variety(
      X, [equation]; check_codimension
    ),
  )
end

Base.show(io::IO, ::HypersurfaceInToricVariety) =
  print(io, "Hypersurface in a normal toric variety")
Base.show(io::IO, Y::CompleteIntersectionInToricVariety) =
  print(io, "Complete-intersection presentation with ",
    length(defining_equations(Y)), " equations in a normal toric variety")

# Preserve inexpensive printing; the generic covered-scheme display would
# calculate affine charts solely to show the presentation.
Base.show(io::IO, ::MIME"text/plain", Y::AbsCompleteIntersectionInToricVariety) =
  show(io, Y)
