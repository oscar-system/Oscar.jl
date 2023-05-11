########################################################################
# This file is to allow all functionality for toric divisors on 
# normal toric varieties also for toric schemes. 
#
# The functionality is simply forwarded from the `underlying_toric_variety` 
# and for now returns an object which lives only in the toric world. 
#
# In the long run, we wish to eventually replace the respective return 
# values by new objects which do contain the original toric object in the 
# internals and implement all the toric functionality inherited from them, 
# but, at the same time, are placed within the schemes branch of the 
# type tree and allow to construct the geometric content of the object. 
########################################################################

function divisor_of_character(X::ToricCoveredScheme, character::Vector{T}) where {T<:Integer}
  return divisor_of_character(underlying_toric_variety(X), character)
end

function toric_divisor(X::ToricCoveredScheme, coeffs::Vector{T}) where {T<:IntegerUnion}
  return toric_divisor(underlying_toric_variety(X), coeffs)
end

function trivial_divisor(X::ToricCoveredScheme)
  return trivial_divisor(underlying_toric_variety(X))
end

function anticanonical_divisor(X::ToricCoveredScheme)
  return anticanonical_divisor(underlying_toric_variety(X))
end

function canonical_divisor(X::ToricCoveredScheme)
  return canonical_divisor(underlying_toric_variety(X))
end

function toric_divisor_class(X::ToricCoveredScheme, class::GrpAbFinGenElem)
  return toric_divisor_class(underlying_toric_variety(X), class)
end

function toric_divisor_class(X::ToricCoveredScheme, coeffs::Vector{T}) where {T<:IntegerUnion}
  return toric_divisor_class(underlying_toric_variety(X), coeffs)
end

function trivial_divisor_class(X::ToricCoveredScheme)
  return trivial_divisor_class(underlying_toric_variety(X))
end

function canonical_divisor_class(X::ToricCoveredScheme)
  return canonical_divisor_class(underlying_toric_variety(X))
end

function anticanonical_divisor_class(X::ToricCoveredScheme)
  return anticanonical_divisor_class(underlying_toric_variety(X))
end

function toric_line_bundle(X::ToricCoveredScheme, class::GrpAbFinGenElem)
  return toric_line_bundle(underlying_toric_variety(X), class)
end

function toric_line_bundle(X::ToricCoveredScheme, c::Vector{T}) where {T<:IntegerUnion}
  return toric_line_bundle(underlying_toric_variety(X), c)
end

function toric_line_bundle(X::ToricCoveredScheme, d::ToricDivisor)
  return toric_line_bundle(underlying_toric_variety(X), d)
end

function anticanonical_bundle(X::ToricCoveredScheme)
  return anticanonical_bundle(underlying_toric_variety(X))
end

function canonical_bundle(X::ToricCoveredScheme)
  return canonical_bundle(underlying_toric_variety(X))
end

# TODO: This is actually a problem, because structure_sheaf is already implemented for 
# AbsCoveredScheme and should not be overwritten. We have to separate or fuse the two. 
function structure_sheaf(X::ToricCoveredScheme)
  return structure_sheaf(underlying_toric_variety(X))
end
