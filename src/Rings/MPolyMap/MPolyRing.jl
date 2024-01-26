# Since coeff_map is at the third position, we wannot use a variable
# with default argument


################################################################################
#
#  Domains of type MPolyRing
#
################################################################################

################################################################################
#
#  Constructos
#
################################################################################

function _check_imgs(S::NCRing, imgs)
  n = length(imgs)
  for i in 2:n, j in 1:(i - 1)
    @req imgs[i] * imgs[j] == imgs[j] * imgs[i] "Images $i and $j do not commute"
  end
  return nothing
end

# no check for commutative codomains
_check_imgs(S::Ring, imgs) = nothing 

@doc raw"""
    hom(R::MPolyRing, S::NCRing, coeff_map, images::Vector; check::Bool = true)

    hom(R::MPolyRing, S::NCRing, images::Vector; check::Bool = true)
    
Given a homomorphism `coeff_map` from `C` to `S`, where `C` is the 
coefficient ring of `R`, and given a vector `images` of `nvars(R)` 
elements of `S`, return the homomorphism `R` $\to$ `S` whose restriction 
to `C` is `coeff_map`, and which sends the `i`-th variable of `R` to the 
`i`-th entry of `images`.
 
If no coefficient map is entered, invoke a canonical homomorphism of `C`
to `S`, if such a homomorphism exists, and throw an error, otherwise.

!!! note
    In case `check = true` (default), the function checks the conditions below:
    - If `S` is graded, the assigned images must be homogeneous with respect to the given grading.
    - If `S` is noncommutative, the assigned images must pairwise commute. 
    

# Examples
```jldoctest
julia> K, a = finite_field(2, 2, "a");

julia> R, (x, y) = polynomial_ring(K, ["x", "y"]);

julia> F = hom(R, R, z -> z^2, [y, x])
Ring homomorphism
  from multivariate polynomial ring in 2 variables over GF(2, 2)
  to multivariate polynomial ring in 2 variables over GF(2, 2)
defined by
  x -> y
  y -> x
with map on coefficients
  #1

julia> F(a * y)
(a + 1)*x

julia> Qi, i = quadratic_field(-1)
(Imaginary quadratic field defined by x^2 + 1, sqrt(-1))

julia> S, (x, y) = polynomial_ring(Qi, ["x", "y"]);

julia> G = hom(S, S, hom(Qi, Qi, -i), [x^2, y^2])
Ring homomorphism
  from multivariate polynomial ring in 2 variables over imaginary quadratic field defined by x^2 + 1
  to multivariate polynomial ring in 2 variables over imaginary quadratic field defined by x^2 + 1
defined by
  x -> x^2
  y -> y^2
with map on coefficients
  Map: imaginary quadratic field defined by x^2 + 1 -> imaginary quadratic field defined by x^2 + 1

julia> G(x+i*y)
x^2 - sqrt(-1)*y^2

julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"]);

julia> f = 3*x^2+2*x+1;

julia> S, (x, y) = polynomial_ring(GF(2), ["x", "y"]);

julia> H = hom(R, S, gens(S))
Ring homomorphism
  from multivariate polynomial ring in 2 variables over ZZ
  to multivariate polynomial ring in 2 variables over GF(2)
defined by
  x -> x
  y -> y

julia> H(f)
x^2 + 1
```
"""
function hom(R::MPolyRing, S::NCRing, coeff_map, images::Vector; check::Bool = true)
  n = ngens(R)
  @req n == length(images) "Number of images must be $n"
  # Now coerce into S or throw an error if not possible
  imgs = _coerce(S, images)
  @check begin
    _check_imgs(S, imgs)
    _check_homo(S, imgs) # defined in MPolyAnyMap.jl
  end
  return MPolyAnyMap(R, S, coeff_map, copy(imgs)) # copy because of #655
end

function hom(R::MPolyRing, S::NCRing, images::Vector; check::Bool = true)
  n = ngens(R)
  @req n == length(images) "Number of images must be $n"
  # Now coerce into S or throw an error if not possible
  imgs = _coerce(S, images)
  @check begin
    _check_imgs(S, imgs)
    _check_homo(S, imgs) # defined in MPolyAnyMap.jl
  end
  return MPolyAnyMap(R, S, nothing, copy(imgs)) # copy because of #655
end



################################################################################
#
#  Evaluation functions
#
################################################################################

function _evaluate_plain(F::MPolyAnyMap{<: MPolyRing}, u)
  return evaluate(u, F.img_gens)
end

function _evaluate_general(F::MPolyAnyMap{<: MPolyRing}, u)
  if domain(F) === codomain(F) && coefficient_map(F) === nothing
    return evaluate(map_coefficients(coefficient_map(F), u,
                                     parent = domain(F)), F.img_gens)
  else
    S = temp_ring(F)
    if S !== nothing
      return evaluate(map_coefficients(coefficient_map(F), u,
                                parent = S), F.img_gens)
    else
      return evaluate(map_coefficients(coefficient_map(F), u), F.img_gens)
    end
  end
end

# one more intermediate function

function _evaluate_help(F::MPolyAnyMap{<: MPolyRing, <: Any, Nothing}, g)
  return _evaluate_plain(F, g)
end

function _evaluate_help(F::MPolyAnyMap{<: MPolyRing}, g)
  return _evaluate_general(F, g)
end

function (F::MPolyAnyMap{<: MPolyRing})(g)
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
