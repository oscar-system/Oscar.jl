################################################################################
#
#  Domains of type MPolyQuoRing
#
################################################################################

function _check_imgs_quo(R, S::NCRing, imgs, coeff_map = nothing)
  n = length(imgs)
  for i in 1:n, j in 1:(i - 1)
    @req imgs[i] * imgs[j] == imgs[j] * imgs[i] "Images $i and $j do not commute"
  end
  gensI = gens(modulus(R))
  for g in gensI
    if coeff_map !== nothing
      @req is_zero(evaluate(map_coefficients(coeff_map, g), imgs)) "Morphism is not well-defined"
    else
      @req iszero(evaluate(g, imgs)) "Morphism is not well-defined"
    end
  end
  return nothing
end

# no check for commutative codomains
function _check_imgs_quo(R, S::Ring, imgs, coeff_map = nothing)
  gensI = gens(modulus(R))
  for g in gensI
    if coeff_map !== nothing
      @req is_zero(evaluate(map_coefficients(coeff_map, g), imgs)) "Morphism is not well-defined"
    else
      @req iszero(evaluate(g, imgs)) "Morphism is not well-defined"
    end
  end
end

@doc raw"""
    hom(A::MPolyQuoRing, S::NCRing, coeff_map, images::Vector; check::Bool = true)

    hom(A::MPolyQuoRing, S::NCRing, images::Vector; check::Bool = true)

Given a homomorphism `coeff_map` from `C` to `S`, where `C` is the 
coefficient ring of the base ring of `A`, and given a vector `images` of `ngens(A)` 
elements of `S`, return the homomorphism `A` $\to$ `S` whose restriction 
to `C` is `coeff_map`, and which sends the `i`-th generator of `A` to the 
`i`-th entry of `images`.
 
If no coefficient map is entered, invoke a canonical homomorphism from `C`
to `S`, if such a homomorphism exists, and throw an error, otherwise.

!!! note
    The function returns a well-defined homomorphism `A` $\to$ `S` iff the
    given data defines a homomorphism `base_ring(A)` $\to$ `S` whose
    kernel contains the modulus of `A`. This condition is checked by the 
    function in case `check = true` (default).

!!! note
    In case `check = true` (default), the function also checks the conditions below:
    - If `S` is graded, the assigned images must be homogeneous with respect to the given grading.
    - If `S` is noncommutative, the assigned images must pairwise commute. 
 
# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]));

julia> S, (s, t) = polynomial_ring(QQ, [:s, :t]);

julia> F = hom(A, S, [s, s^2, s^3])
Ring homomorphism
  from quotient of multivariate polynomial ring by ideal (-x^2 + y, -x^3 + z)
  to multivariate polynomial ring in 2 variables over QQ
defined by
  x -> s
  y -> s^2
  z -> s^3
```
"""
function hom(R::MPolyQuoRing, S::NCRing, coeff_map, images::Vector; check::Bool = true)
  n = ngens(R)
  @req n == length(images) "Number of images must be $n"
  # Now coerce into S or throw an error if not possible
  imgs = _coerce(S, images)
  @check begin
    _check_imgs_quo(R, S, imgs, coeff_map)
    _check_homo(S, imgs)
  end

  return MPolyAnyMap(R, S, coeff_map, copy(imgs)) # copy because of #655
end

function hom(R::MPolyQuoRing, S::NCRing, images::Vector; check::Bool = true)
  n = ngens(R)
  @req n == length(images) "Number of images must be $n"
  # Now coerce into S or throw an error if not possible
  imgs = _coerce(S, images)
  @check begin
    _check_imgs_quo(R, S, imgs)
    _check_homo(S, imgs)
  end
  return MPolyAnyMap(R, S, nothing, copy(imgs)) # copy because of #655
end

################################################################################
#
#  Evaluation
#
################################################################################

function _evaluate_plain(F::MPolyAnyMap{<: MPolyQuoRing}, u)
  return evaluate(lift(u), _images(F))
end

function _evaluate_plain(F::MPolyAnyMap{<:MPolyQuoRing, <:MPolyQuoRing}, u)
  # This workaround deals with the fact that arithmetic in quotient rings is TERRIBLY slow.
  # All the simplify calls make it unusable in this case, probably due to the fact that
  # setting `is_reduced` flags does not pay off in such iterative procedures.
  A = codomain(F)
  # Shortcut needed because the `evaluate` code can not deal with empty lists
  is_empty(_images(F)) && return codomain(F)(_constant_coefficient(u)) 
  v = evaluate(lift(u), lift.(_images(F)))
  return simplify(A(v))
end

_constant_coefficient(f::MPolyRingElem) = constant_coefficient(f)
_constant_coefficient(f::MPolyQuoRingElem) = constant_coefficient(lift(f))

function _evaluate_general(F::MPolyAnyMap{<: MPolyQuoRing}, u)
  if domain(F) === codomain(F) && coefficient_map(F) === nothing
    return evaluate(map_coefficients(coefficient_map(F), lift(u),
                                     parent = domain(F)), F.img_gens)
  else
    S = temp_ring(F)
    if S !== nothing
      return evaluate(map_coefficients(coefficient_map(F), lift(u),
                                                        parent = S), F.img_gens)
    else
      return evaluate(map_coefficients(coefficient_map(F), lift(u)), F.img_gens)
    end
  end
end

# one more intermediate function

function _evaluate_help(F::MPolyAnyMap{<: MPolyQuoRing, <: Any, Nothing}, g)
  return _evaluate_plain(F, g)
end

function _evaluate_help(F::MPolyAnyMap{<: MPolyQuoRing}, g)
  return _evaluate_general(F, g)
end

# The two main evaluation methods
function (F::MPolyAnyMap{<: MPolyQuoRing})(g)
  if g isa elem_type(domain(F))
    @req parent(g) === domain(F) "Element not in domain"
    if coefficient_map(F) === nothing
      return _evaluate_plain(F, g)
    else 
      return _evaluate_general(F, g)
    end
  else 
    gg = domain(F)(g)
    @assert parent(gg) === domain(F)
    return F(gg)
  end
end
