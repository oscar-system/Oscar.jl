export hasse_derivatives

### Implementation of Hasse-Schmidt derivatives as seen in 
###
###     Fruehbis-Krueger, Ristau, Schober: 'Embedded desingularization for arithmetic surfaces -- toward a parallel implementation'


################################################################################
### HASSE-SCHMIDT derivatives for single polynomials

@doc raw"""
    hasse_derivatives(f::MPolyRingElem)

Return Hasse-Schmidt derivatives of `f`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"]);

julia> f = R(5*x^2 + 3*y^5);

julia> hasse_derivatives(f)
9-element Vector{ZZMPolyRingElem}:
 5*x^2 + 3*y^5
 15*y^4
 30*y^3
 30*y^2
 15*y
 3
 10*x
 5
```
"""
function hasse_derivatives(f::MPolyRingElem)
  R = parent(f)
  n = ngens(R)
  # degreef = maximum(degrees(f))
  # define new ring with more variables: R[x1, ..., xn] -> R[x1, ..., xn, t1, ..., tn]
  Rtemp, _ = polynomial_ring(R, "y" => 1:n, "t" => 1:n)
  # replace f(x_i) -> f(y_i + t_i)
  F = evaluate(f, gens(Rtemp)[1:n] + gens(Rtemp)[n+1:2n])
  # i = 1 # counter to iterate though degrees of monomials
  # HasseDerivativesList = empty([f])
  HasseDerivativesList = [f]
  varR = vcat(gens(R), ones(typeof(base_ring(R)(1)), n))
  # getting hasse derivs without extra attention on ordering
  for term in terms(F)
    if sum(degrees(term)[n+1:2n]) != 0
      # hasse derivatives are the factors in front of the monomial in t
      push!(HasseDerivativesList, evaluate(term, varR))
    end
  end
  return HasseDerivativesList
end

function hasse_derivatives(f::MPolyQuoRingElem)
  error("Not implemented. For experts, however, there is an internal function called _hasse_derivatives, which works for elements of type MPolyQuoRingElem")
end

function hasse_derivatives(f::Oscar.MPolyLocRingElem)
  error("Not implemented. For experts, however, there is an internal function called _hasse_derivatives, which works for elements of type Oscar.MPolyLocRingElem")
end

function hasse_derivatives(f::Oscar.MPolyQuoLocRingElem)
  error("Not implemented. For experts, however, there is an internal function called _hasse_derivatives, which works for elements of type Oscar.MPolyQuoLocRingElem")
end


################################################################################
### HASSE-SCHMIDT derivatives for a list of polynomials

@doc raw"""
    hasse_derivatives(v::Vector{MPolyRingElem})

For every `f` in `v`: Return a list of all Hasse-Schmidt derivatives of the lifted numerator of `f`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"]);

julia> f1 = R(3*x^2);

julia> f2 = R(6*y^3);

julia> hasse_derivatives([f1, f2])
2-element Vector{Vector{ZZMPolyRingElem}}:
 [3*x^2, 6*x, 3]
 [6*y^3, 18*y^2, 18*y, 6]
```
"""
function hasse_derivatives(v::Vector{<:MPolyRingElem})
  return hasse_derivatives.(v)
end




################################################################################
### internal functions for expert use

# MPolyQuoRingElem (internal, expert use only)
@doc raw"""
    _hasse_derivatives(f::MPolyQuoRingElem)

Return Hasse-Schmidt derivatives of lifted numerator of `f`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(ZZ, 4, "x");

julia> I = ideal(R, [x[2] - 1]);

julia> RQ, _ = quo(R, I);

julia> f = RQ(2*x[2]^4);

julia> _hasse_derivatives(f)
5-element Vector{ZZMPolyRingElem}:
 2*x2^4
 8*x2^3
 12*x2^2
 8*x2
 2
```
"""
function _hasse_derivatives(f::MPolyQuoRingElem)
  return hasse_derivatives(lifted_numerator(f)) 
end

# Oscar.MPolyLocRingElem (internal, expert use only)
@doc raw"""
    _hasse_derivatives(f::Oscar.MPolyLocRingElem)

Return Hasse-Schmidt derivatives of numerator of `f`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, 4, "x");

julia> m = ideal(R, [x[1] - 3, x[2] - 2, x[3] + 2, x[4]]);

julia> U = complement_of_prime_ideal(m);

julia> Rloc, _ = localization(R, U);

julia> f = Rloc(2*x[4]^5);

julia> _hasse_derivatives(f)
9-element Vector{QQMPolyRingElem}:
 2*x4^5
 10*x4^4
 20*x4^3
 20*x4^2
 10*x4
 2
```
"""
function _hasse_derivatives(f::Oscar.MPolyLocRingElem)
  return hasse_derivatives(numerator(f)) 
end

# Oscar.MPolyQuoLocRingElem (internal, expert use only)
@doc raw"""
    _hasse_derivatives(f::Oscar.MPolyQuoLocRingElem)

Return Hasse-Schmidt derivatives of lifted numerator of `f`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, 4, "x");

julia> I = ideal(R, [x[1]^3 - 1]);

julia> RQ, _ = quo(R, I);

julia> p = ideal(R, [x[2]]);

julia> U = complement_of_prime_ideal(p);

julia> RQL, _ = localization(RQ, U);

julia> f = RQL(4*x[3]^3);

julia> _hasse_derivatives(f)
4-element Vector{QQMPolyRingElem}:
 4*x3^3
 12*x3^2
 12*x3
 4
```
"""
function _hasse_derivatives(f::Oscar.MPolyQuoLocRingElem)
  return hasse_derivatives(lifted_numerator(f))
end

# for a list of elements (internal, expert use only)
@doc raw"""
    _hasse_derivatives(v::Vector{RingElem}})

For every `f` in `v`: Return a list of all Hasse-Schmidt derivatives of the lifted numerator of `f`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, 4, "x");

julia> I = ideal(R, [x[1]^3 - 1]);

julia> RQ, _ = quo(R, I);

julia> p = ideal(R, [x[2]]);

julia> U = complement_of_prime_ideal(p);

julia> RQL, _ = localization(RQ, U);

julia> f1 = RQ(4*x[3]^3);

julia> f2 = RQL(3*x[2]^2);

julia> _hasse_derivatives([f1, f2])
2-element Vector{Vector{QQMPolyRingElem}}:
 [4*x3^3, 12*x3^2, 12*x3, 4]
 [3*x2^2, 6*x2, 3]
```
"""
function _hasse_derivatives(v::Vector{RingElem}) # type of input has not to be very strict, because it gets checked for each element of v
  return _hasse_derivatives.(v)
end
