#####################################################################
# 1: Some simple constructors
#####################################################################

family_of_spaces(coordinate_ring::MPolyRing, grading::Matrix{Int64}, dim::Int) = FamilyOfSpaces(coordinate_ring, grading, dim)


#####################################################################
# 2: Display
#####################################################################

function Base.show(io::IO, f::FamilyOfSpaces)
  print(io, "A family of spaces of dimension d = $(dim(f))")
end
