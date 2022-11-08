export empty_covered_scheme

########################################################################
# Constructors for CoveredScheme                                       #
########################################################################

### The standard constructor
@Markdown.doc """
    CoveredScheme(C::Covering)

Return a `CoveredScheme` ``X`` with `C` as its `default_covering`.

TODO: Add an example.
"""
function CoveredScheme(C::Covering)
  refinements = Dict{Tuple{Covering, Covering}, CoveringMorphism}()
  X = CoveredScheme([C], refinements)
  set_attribute!(X, :seed_covering, C)
  return X
end

### Conversion of an affine scheme into a covered scheme
CoveredScheme(X::AbsSpec) = CoveredScheme(Covering(X))

### Construct the empty covered scheme over the ring R
function empty_covered_scheme(R::RT) where {RT<:AbstractAlgebra.Ring}
  return CoveredScheme(R)
end

