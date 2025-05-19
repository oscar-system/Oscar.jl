#####################################################################
# 1: Some simple constructors
#####################################################################

@doc raw"""
    family_of_spaces(coordinate_ring::MPolyRing, grading::Matrix{Int64}, dim::Int)

Return a family of spaces that is (currently) used to build
families of F-theory models, defined by using a family of
base spaces for an elliptic fibration. It is specified by
the following data:
* A polynomial ring. This may be thought of as the coordinate
  ring of the generic member in this family of spaces.
* A grading for this polynomial ring. This may be thought of as
  specifying certain line bundles on the generic member in this
  family of spaces. Of particular interest for F-theory is always
  the canonical bundle. For this reason, the first row is always
  thought of as grading the coordinate ring with regard to the
  canonical bundle. Or put differently, if a variable shares the
  weight 1 in the first row, then this means that we should think
  of this coordinate as a section of 1 times the canonical bundle.
* An integer, which specifies the dimension of the generic member
  within this family of spaces.

Note that the coordinate ring can have strictly more variables than
the dimension. This is a desired feature for most, if not all,
F-theory literature constructions.

```jldoctest
julia> coord_ring, _ = QQ[:f, :g, :Kbar, :u];

julia> grading = [4 6 1 0; 0 0 0 1]
2×4 Matrix{Int64}:
 4  6  1  0
 0  0  0  1

julia> d = 3
3

julia> f = family_of_spaces(coord_ring, grading, d)
Family of spaces of dimension d = 3
```
"""
family_of_spaces(coordinate_ring::MPolyRing, grading::Matrix{Int64}, dim::Int) = FamilyOfSpaces(coordinate_ring, grading, dim)


#####################################################################
# 2: Display
#####################################################################

# Detailed printing
function Base.show(io::IO, ::MIME"text/plain", f::FamilyOfSpaces)
  io = pretty(io)
  print(io, "Family of spaces of dimension d = $(dim(f))")
end

# Terse and one line printing
function Base.show(io::IO, f::FamilyOfSpaces)
  print(io, "Family of spaces")
end
