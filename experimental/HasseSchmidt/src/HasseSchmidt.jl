export hasse_derivatives

### We consider Hasse-Schmidt derivatives of polynomials as seen in
###
###   	[FRS21](@cite) Fruehbis-Krueger, Ristau, Schober: 'Embedded desingularization for arithmetic surfaces -- toward a parallel implementation'
###			
### This is a special case of a more general definition of a Hasse-Schmidt derivative. These more general and rigorous definitions can be found in the following sources:
###
###			[Cut04](@cite) Cutkosky: 'Resolution of Singularities'
###			[Haz12](@cite) Michiel Hazewinkel: 'Hasse-Schmidt derivations and the Hopf algebra of noncommutative symmetric functions'
###

################################################################################
### HASSE-SCHMIDT derivatives for single polynomials

@doc raw"""
    hasse_derivatives(f::MPolyRingElem)

Return a list of Hasse-Schmidt derivatives of `f`, each with a multiindex `[a_1, ..., a_n]`, where `a_i` describes the number of times `f` was derived w.r.t. the `i`-th variable. 

Hasse-Schmidt derivatives as seen in [FRS21](@cite).
For more general and rigorous definition see [Cut04](@cite) or [Haz12](@cite).

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"]);

julia> f = 5*x^2 + 3*y^5;

julia> hasse_derivatives(f)
8-element Vector{Vector{Any}}:
 [[0, 5], 3]
 [[0, 4], 15*y]
 [[0, 3], 30*y^2]
 [[2, 0], 5]
 [[0, 2], 30*y^3]
 [[1, 0], 10*x]
 [[0, 1], 15*y^4]
 [[0, 0], 5*x^2 + 3*y^5]
```
"""
function hasse_derivatives(f::MPolyRingElem)
  R = parent(f)
  Rtemp, t = polynomial_ring(R, :t => 1:ngens(R))
  F = evaluate(f, gens(R) + t)
  return [[degrees(monomial(term, 1)), coeff(term, 1)] for term in terms(F)]
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
### internal functions for expert use

# MPolyQuoRingElem (internal, expert use only)
@doc raw"""
    _hasse_derivatives(f::MPolyQuoRingElem)

Return a list of Hasse-Schmidt derivatives of lift of `f`, each with a multiindex `[a_1, ..., a_n]`, where `a_i` describes the number of times `f` was derived w.r.t. the `i`-th variable.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"]);

julia> I = ideal(R, [x - 1]);

julia> RQ, phi = quo(R, I);

julia> f = phi(2*y^4);

julia> Oscar._hasse_derivatives(f)
5-element Vector{Vector{Any}}:
 [[0, 4], 2]
 [[0, 3], 8*y]
 [[0, 2], 12*y^2]
 [[0, 1], 8*y^3]
 [[0, 0], 2*y^4]
```
"""
function _hasse_derivatives(f::MPolyQuoRingElem)
  return hasse_derivatives(lift(f)) 
end

# Oscar.MPolyLocRingElem (internal, expert use only)
@doc raw"""
    _hasse_derivatives(f::Oscar.MPolyLocRingElem)

Return a list of Hasse-Schmidt derivatives of numerator of `f`, each with a multiindex `[a_1, ..., a_n]`, where `a_i` describes the number of times `f` was derived w.r.t. the `i`-th variable.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> m = ideal(R, [x - 3, y - 2, z + 1]);

julia> U = complement_of_prime_ideal(m);

julia> Rloc, phi = localization(R, U);

julia> f = phi(2*z^5);

julia> Oscar._hasse_derivatives(f)
6-element Vector{Vector{Any}}:
 [[0, 0, 5], 2]
 [[0, 0, 4], 10*z]
 [[0, 0, 3], 20*z^2]
 [[0, 0, 2], 20*z^3]
 [[0, 0, 1], 10*z^4]
 [[0, 0, 0], 2*z^5]
```
"""
function _hasse_derivatives(f::Oscar.MPolyLocRingElem)
  return hasse_derivatives(numerator(f)) 
end

# Oscar.MPolyQuoLocRingElem (internal, expert use only)
@doc raw"""
    _hasse_derivatives(f::Oscar.MPolyQuoLocRingElem)

Return a list of Hasse-Schmidt derivatives of lifted numerator of `f`, each with a multiindex `[a_1, ..., a_n]`, where `a_i` describes the number of times `f` was derived w.r.t. the `i`-th variable.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> I = ideal(R, [x^3 - 1]);

julia> RQ, phi = quo(R, I);

julia> p = ideal(R, [z]);

julia> U = complement_of_prime_ideal(p);

julia> RQL, iota = localization(RQ, U);

julia> f = iota(phi(4*y^3));

julia> Oscar._hasse_derivatives(f)
4-element Vector{Vector{Any}}:
 [[0, 3, 0], 4]
 [[0, 2, 0], 12*y]
 [[0, 1, 0], 12*y^2]
 [[0, 0, 0], 4*y^3]
```
"""
function _hasse_derivatives(f::Oscar.MPolyQuoLocRingElem)
  return hasse_derivatives(lifted_numerator(f))
end
