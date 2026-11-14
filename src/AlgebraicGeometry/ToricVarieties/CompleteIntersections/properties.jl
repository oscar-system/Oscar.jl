# Properties and validation of toric complete-intersection presentations.

@doc raw"""
    has_expected_codimension(Y::CompleteIntersectionInToricVariety)

Test whether the associated scheme is empty or has the expected codimension.
Use the standard scheme methods to check emptiness and compute dimension.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 2);

julia> x, y, z = gens(cox_ring(X));

julia> C = complete_intersection_in_toric_variety(X, [x, y]);

julia> has_expected_codimension(C)
true
```
"""
@attr Bool function has_expected_codimension(Y::CompleteIntersectionInToricVariety)
  return is_empty(Y) ||
         dim(Y) == dim(toric_variety(Y)) - length(defining_equations(Y))
end
has_expected_codimension(Y::HypersurfaceInToricVariety) =
  has_expected_codimension(_complete_intersection_presentation(Y))

function _check_complete_intersection_presentation(
  Y::AbsCompleteIntersectionInToricVariety;
  check_codimension::Bool=false,
)
  if check_codimension
    @req has_expected_codimension(Y) "The equations do not have the expected codimension"
    if !is_empty(Y)
      X = toric_variety(Y)
      @req is_smooth(X) ||
        all(D -> is_cartier(toric_divisor(D)), defining_divisor_classes(Y)) "Certified complete intersections on singular toric varieties require Cartier defining classes"
    end
  end
  return nothing
end

# All standard geometric scheme properties (is_empty, is_smooth, is_reduced,
# is_irreducible, is_integral, etc.) are inherited through AbsCoveredScheme.
# Only the presentation-specific codimension predicate is defined above.
