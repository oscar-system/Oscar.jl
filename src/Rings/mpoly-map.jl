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

@attributes mutable struct MPolyMap{D <: MPolyRing, C <: NCRing, U, V} <: Map{D, C, Hecke.Hecke.Map, MPolyMap}
  domain::D
  codomain::C
  coeff_map::U
  img_gens::Vector{V}

  function MPolyMap{D, C, U, V}(domain::D,
                                codomain::C,
                                coeff_map::U,
                                img_gens::Vector{V}) where {D, C, U, V}
      @assert V === elem_type(C)
    return new{D, C, U, V}(domain, codomain, coeff_map, img_gens)
  end
end

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
#  Constructos
#
################################################################################

function hom(R::MPolyRing, S::NCRing, img_gens::Vector; check::Bool = true)
  n = ngens(R)
  @req n == length(img_gens) "Number of images must be $n"
  for i in 1:n
    for j in 1:(i - 1)
      @req img_gens[i] * img_gens[j] ==
                               img_gens[j] * img_gens[i] "Images do not commute"
    end
  end

  return MPolyMap{typeof(R), typeof(S), Nothing, eltype(img_gens)}(R, S, nothing, img_gens)
end

function hom(R::MPolyRing, S::Ring, img_gens::Vector)
  n = ngens(R)
  @req n == length(img_gens) "Number of images must be $n"

  if eltype(img_gens) !== elem_type(S)
    _img_gens = S.(img_gens)::Vector{elem_type(S)} # promote to S
  else
    _img_gens = img_gens::Vector{elem_type(S)}
  end

  return MPolyMap{typeof(R), typeof(S), Nothing, elem_type(S)}(R, S, nothing, _img_gens)
end

function hom(R::MPolyRing, S::Ring, coeff_map, img_gens::Vector)
  if eltype(img_gens) !== elem_type(S)
    _img_gens = S.(img_gens)::Vector{elem_type(S)} # promote to S
  else
    _img_gens = img_gens::Vector{elem_type(S)}
  end

  return MPolyMap{typeof(R), typeof(S), typeof(coeff_map), elem_type(S)}(R, S, coeff_map, _img_gens)
end

################################################################################
#
#  Evaluation functions
#
################################################################################

function _evaluate_plain(F::MPolyMap, u)
  return evaluate(u, F.img_gens)
end

function _evaluate_general(F::MPolyMap, u)
  return evaluate(map_coefficients(F.coeff_map, u), F.img_gens)
end

function (F::MPolyMap{<: Any, <: Any, U, <:Any})(g) where U
  if U === Nothing
    return _evaluate_plain(F, g)
  else
    return _evaluate_general(F, g)
  end
end

################################################################################
#
#  More complicated functions
#
################################################################################

# If you want to implement a kernel function for morphisms QQ[x] -> QQ[x],
# this will be implemented using dispatch as follows:

function kernel(f::MPolyMap{D}) where {D <: MPolyRing{fmpq}}
  # Now do something
  #
  # If you need to call Singular, use the get_attribute! stuff to
  # save some Singular data etc.
end
