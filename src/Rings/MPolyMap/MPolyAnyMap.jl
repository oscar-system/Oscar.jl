# A morphism from a polynomial ring phi : R[x_1,...,x_n] -> S into any
# (noncommutative) ring is uniquely determined by the restriction of phi to K
# and the images of the x_i under phi.
#
# Thus, the datum for the morphism from a polynomial ring into a ring consists of
#
# - the morphism restricted to the coefficient ring (coeff_map)
# - the images of the generators (img_gens)
#
# To ease the usage we also add the domain and the codomain.
#
# In case there is a canonical embedding R -> S (which happens for example
# if R === S, or R === ZZ and S arbitrary, we store the coefficient ring map
# as nothing.
#
# Evaluating a multivariate polynomial f  under the map is thus either
#
# evaluate((coeff_map(f), img_gens)
#
# or just
#
# evaluate(f, img_gens)
# 
# Note that the evaluate function is designed to work whenever there are
# canonical embeddings in place (aka S(...))
#
# Finally note that the same strategy works for
# - MPolyRing
# - MPolyQuo
# - MPolyRing_dec

const _DomainTypes = Union{MPolyRing, MPolyQuo}

@attributes mutable struct MPolyAnyMap{
    D <: _DomainTypes,
    C <: NCRing,
    U,
    V} <: Map{D, C, Hecke.Hecke.Map, MPolyAnyMap}

  domain::D
  codomain::C
  coeff_map::U
  img_gens::Vector{V}
  temp_ring           # temporary ring used when evaluating maps

  function MPolyAnyMap{D, C, U, V}(domain::D,
                                codomain::C,
                                coeff_map::U,
                                img_gens::Vector{V}) where {D, C, U, V}
      @assert V === elem_type(C)
      for g in img_gens
        @assert parent(g) === codomain "elements does not have the correct parent"
      end
    return new{D, C, U, V}(domain, codomain, coeff_map, img_gens)
  end
end

function MPolyAnyMap(d::D, c::C, cm::U, ig::Vector{V}) where {D, C, U, V}
  return MPolyAnyMap{D, C, U, V}(d, c, cm, ig)
end

################################################################################
#
#  Field access
#
################################################################################

# domain, codomain are automatically defined by the general map interface

# Not sure if we want to expose the following function to the user.
# It might be `nothing`. We could return `identity in the `nothing` case.
coefficient_map(f::MPolyAnyMap) = f.coeff_map

_images(f::MPolyAnyMap) = f.img_gens

################################################################################
#
#  String I/O
#
################################################################################

# There is some default printing for maps.
# We might want to hijack this if we want to indicate whether there is
# a non-trivial coefficient ring map.

################################################################################
#
#  Helper
#
################################################################################

# # Since we want to allow specifiying images in a "subring", we need to coerce
# if necessary. For example, hom(Qx, Qx, [1, 1]), should work, although
# 1 is not an element of the codomain.
function _coerce(S, img_gens)
  if eltype(img_gens) === elem_type(S)
    return img_gens::Vector{elem_type(S)}
  else
    _img_gens = S.(img_gens)
    if eltype(_img_gens) !== elem_type(S)
      @req false "Elements cannot be coerced into the codomain"
    end
    return _img_gens::Vector{elem_type(S)}
  end
end

# The following is used when the domain is graded

# if the codomain is graded, the images must be homogenous?!
_isgraded(::Ring) = false
_isgraded(R::MPolyRing_dec) = isgraded(R) # filtered not considered graded
_isgraded(R::MPolyQuo) = _isgraded(R.R)

function _check_homo(S::Ring, images)
  if _isgraded(S)
    for i in images
      @req ishomogeneous(i) "Images must be homogenous"
    end
  end
end

# When evaluating the map phi at a polynomial f, we first construct the polynomial
# map_coefficients(phi.coeff_map, f), which is a polynomial over
# codomain(phi.coeff_map).

function temp_ring(f::MPolyAnyMap{<:Any, <: Any, <: Map{D, C}}) where {D, C}
  if isdefined(f, :temp_ring)
    return f.temp_ring::mpoly_ring_type(C)
  end

  S, = PolynomialRing(codomain(coefficient_map(f)), nvars(domain(f)))
  f.temp_ring = S
  return S
end

function temp_ring(f::MPolyAnyMap{<:Any, <: Any})
  return nothing
end

################################################################################
#
#  Kernel
#
################################################################################

function kernel(f::MPolyAnyMap)
  error("Cannot compute this kernel!")
end

################################################################################
#
#  Injectivity
#
################################################################################

function isinjective(F::MPolyAnyMap)
  error("Cannot decide injectivity!")
end

################################################################################
#
#  Surjectivity
#
################################################################################

function issurjective(F::MPolyAnyMap)
  error("Cannot decide surjectivity!")
end

################################################################################
#
#  Bijectivity
#
################################################################################

function isbijective(F::MPolyAnyMap)
  error("Cannot decide bijectivity!")
end
