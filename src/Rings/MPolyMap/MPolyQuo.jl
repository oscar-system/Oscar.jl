################################################################################
#
#  Domains of type MPolyQuo
#
################################################################################

function _check_imgs_quo(R, S::NCRing, imgs)
  n = length(imgs)
  for i in 1:n, j in 1:(i - 1)
    @req imgs[i] * imgs[j] == imgs[j] * imgs[i] "Images $i and $j do not commute"
  end
  gensI = gens(modulus(R))
  for g in gensI
    @req iszero(evaluate(g, imgs)) "Morphism is not well-defined"
  end
  return nothing
end

# no check for commutative codomains
function _check_imgs_quo(R, S::Ring, imgs) 
  gensI = gens(modulus(R))
  for g in gensI
    @req iszero(evaluate(g, imgs)) "Morphism is not well-defined"
  end
end

@doc Markdown.doc"""
    hom(A::MPolyQuo, S::NCRing, coeff_map, images::Vector; check::Bool = true)

    hom(A::MPolyQuo, S::NCRing, images::Vector; check::Bool = true)
    
Given a homomorphism `coeff_map` from `C` to `S`, where `C` is the 
coefficient ring of the base ring of `A`, and given a vector `images` of `ngens(A)` 
elements of `S`, return the homomorphism `A` $\to$ `S` whose restriction 
to `C` is `coeff_map`, and which sends the `i`-th generator of `A` to the 
`i`-th entry of `images`.
 
If no coefficient map is entered, invoke a canonical homomorphism of `C`
to `S`, if such a homomorphism exists, and throw an error, otherwise.

!!! note
    The function returns a well-defined homomorphism `A` $\to$ `S` iff the
    given data defines a homomorphism from the base ring of `A` to `S` whose
    kernel contains the modulus of `A`. This condition is checked by the 
    function in case `check = true` (default).

!!! note
    If `S` is noncommutative, the assigned images must pairwise commute.
    This is checked by the function in case `check = true` (default).

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"] );

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]));

julia> S, (s, t) = PolynomialRing(QQ, ["s", "t"]);

julia> F = hom(A, S, [s, s^2, s^3])
Map with following data
Domain:
=======
Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z)
Codomain:
=========
Multivariate Polynomial Ring in s, t over Rational Field
```
"""
function hom(R::MPolyQuo, S::NCRing, coeff_map, images::Vector; check::Bool = true)
  n = ngens(R)
  @req n == length(images) "Number of images must be $n"
  # Now coerce into S or throw an error if not possible
  imgs = _coerce(S, images)
  if check
    _check_imgs_quo(R, S, imgs)
    _check_homo(S, imgs)
  end

  return MPolyAnyMap(R, S, coeff_map, copy(imgs)) # copy because of #655
end

function hom(R::MPolyQuo, S::NCRing, images::Vector; check::Bool = true)
  n = ngens(R)
  @req n == length(images) "Number of images must be $n"
  # Now coerce into S or throw an error if not possible
  imgs = _coerce(S, images)
  if check
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

function _evaluate_plain(F::MPolyAnyMap{<: MPolyQuo}, u)
  return evaluate(lift(u), _images(F))
end

function _evaluate_general(F::MPolyAnyMap{<: MPolyQuo}, u)
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

function _evaluate_help(F::MPolyAnyMap{<: MPolyQuo, <: Any, Nothing}, g)
  return _evaluate_plain(F, g)
end

function _evaluate_help(F::MPolyAnyMap{<: MPolyQuo}, g)
  return _evaluate_general(F, g)
end

# The two main evaluation methods
function (F::MPolyAnyMap{<: MPolyQuo})(g)
  if g isa elem_type(domain(F))
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
