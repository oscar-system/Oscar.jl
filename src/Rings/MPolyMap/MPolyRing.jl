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
```jldoctest; filter = r"\#\d+" => "#"
julia> K, a = finite_field(2, 2, "a");

julia> R, (x, y) = polynomial_ring(K, [:x, :y]);

julia> F = hom(R, R, z -> z^2, [y, x])
Ring homomorphism
  from multivariate polynomial ring in 2 variables over K
  to multivariate polynomial ring in 2 variables over K
defined by
  x -> y
  y -> x
with map on coefficients
  #1

julia> F(a * y)
(a + 1)*x

julia> Qi, i = quadratic_field(-1)
(Imaginary quadratic field defined by x^2 + 1, sqrt(-1))

julia> S, (x, y) = polynomial_ring(Qi, [:x, :y]);

julia> G = hom(S, S, hom(Qi, Qi, -i), [x^2, y^2])
Ring homomorphism
  from multivariate polynomial ring in 2 variables over Qi
  to multivariate polynomial ring in 2 variables over Qi
defined by
  x -> x^2
  y -> y^2
with map on coefficients
  Map: Qi -> Qi

julia> G(x+i*y)
x^2 - sqrt(-1)*y^2

julia> R, (x, y) = polynomial_ring(ZZ, [:x, :y]);

julia> f = 3*x^2+2*x+1;

julia> S, (x, y) = polynomial_ring(GF(2), [:x, :y]);

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

# Some additional methods needed for the test in the constructor for MPolyAnyMap
_is_gen(x::MPolyQuoRingElem) = _is_gen(lift(x))
_is_gen(x::MPolyDecRingElem) = is_gen(forget_grading(x))
_is_gen(x::MPolyRingElem) = is_gen(x)
# default method; overwrite if you want this to work for your rings.
_is_gen(x::NCRingElem) = false

# In case there is a type of ring elements for which hashing is correctly implemented 
# and does not throw an error, this gives the opportunity to overwrite the `allunique`
# to be used within the constructor for maps. 
function _allunique(lst::Vector{T}) where {T<:MPolyRingElem}
  return allunique(lst)
end

# We have a lot of rings which do/can not implement correct hashing. 
# So we make the following the default.
function _allunique(lst::Vector{T}) where {T<:RingElem}
  return all(!(x in lst[i+1:end]) for (i, x) in enumerate(lst))
end

function _allunique(lst::Vector{T}) where {T<:MPolyQuoRingElem}
  rep_list = lift.(lst)
  return _allunique(rep_list)
end

function _build_poly(u::MPolyRingElem, indices::Vector{Int}, S::MPolyRing)
  kk = coefficient_ring(S)
  r = ngens(S)
  ctx = MPolyBuildCtx(S)
  for (c, e) in zip(AbstractAlgebra.coefficients(u), AbstractAlgebra.exponent_vectors(u))
    ee = [0 for _ in 1:r]
    for (i, k) in enumerate(e)
      ee[indices[i]] = k
    end
    push_term!(ctx, kk(c), ee)
  end
  return finish(ctx)
end

function _evaluate_plain(F::MPolyAnyMap{<:MPolyRing, <:MPolyRing}, u)
  fl, var_ind = _maps_variables_to_variables(F)
  if fl
    return _build_poly(u, var_ind, codomain(F))
  end
  return evaluate(u, F.img_gens)
end

function _evaluate_plain(F::MPolyAnyMap{<: MPolyRing}, u)
  return evaluate(u, F.img_gens)
end

# See the comment in MPolyQuo.jl
function _evaluate_plain(F::MPolyAnyMap{<:MPolyRing, <:MPolyQuoRing}, u)
  A = codomain(F)
  R = base_ring(A)
  fl, var_ind = _maps_variables_to_variables(F)
  if fl
    return A(_build_poly(u, var_ind, R))
  end
  v = evaluate(lift(u), lift.(_images(F)))
  return simplify(A(v))
end

# The following assumes `p` to be in `S[x₁,…,xₙ]` where `S` is the 
# actual codomain of the map.
function _evaluate_with_build_ctx(
    p::MPolyRingElem, ind::Vector{Int}, cod_ring::MPolyRing
  )
  @assert cod_ring === coefficient_ring(parent(p))
  r = ngens(cod_ring)
  kk = coefficient_ring(cod_ring)
  ctx = MPolyBuildCtx(cod_ring)
  for (q, e) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
    ee = [0 for _ in 1:r]
    for (i, k) in enumerate(e)
      ee[ind[i]] = k
    end
    for (c, d) in zip(AbstractAlgebra.coefficients(q), AbstractAlgebra.exponent_vectors(q))
      push_term!(ctx, kk(c), ee+d)
    end
  end
  return finish(ctx)
end

function _evaluate_general(F::MPolyAnyMap{<:MPolyRing, <:MPolyRing}, u)
  if domain(F) === codomain(F) && coefficient_map(F) === nothing
    return evaluate(map_coefficients(coefficient_map(F), u,
                                     parent = domain(F)), F.img_gens)
  else
    S = temp_ring(F)
    fl, var_ind = _maps_variables_to_variables(F)
    if S !== nothing
      if !fl || coefficient_ring(S) !== codomain(F)
        return evaluate(map_coefficients(coefficient_map(F), u,
                                         parent = S), F.img_gens)
      else
        tmp_poly = map_coefficients(coefficient_map(F), u, parent = S)
        return _evaluate_with_build_ctx(tmp_poly, var_ind, codomain(F))
      end
    else
      if !fl
        return evaluate(map_coefficients(coefficient_map(F), u), F.img_gens)
      else
        # For the case where we can recycle the method above, do so.
        tmp_poly = map_coefficients(coefficient_map(F), u)
        coefficient_ring(parent(tmp_poly)) === codomain(F) && return _evaluate_with_build_ctx(
                   tmp_poly,
                   var_ind,
                   codomain(F)
                 )
        # Otherwise default to the standard evaluation for the time being.
        return evaluate(tmp_poly, F.img_gens)
      end
    end
  end
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
  parent(g) === domain(F) || return F(domain(F)(g))
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
