export hasse_derivatives


################################################################################
### HASSE-SCHMIDT derivatives for single polynomials

@doc raw"""
    hasse_derivatives(f::MPolyRingElem)

Return the Hasse-Schmidt derivatives of `f`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"]);

julia> f = R(5*x^2 + 3*y^5);

julia> hasse_derivatives(f)
9-element Vector{ZZMPolyRingElem}:
 3*y^5
 15*y^4
 30*y^3
 30*y^2
 15*y
 3
 5*x^2
 10*x
 5
```
"""

#  MPolyRingElem
function hasse_derivatives(f::MPolyRingElem)
  R = parent(f)
  n = ngens(R)
  # degreef = maximum(degrees(f))
  # define new ring with more variables: R[x1, ..., xn] -> R[x1, ..., xn, t1, ..., tn]
  Rtemp, _ = polynomial_ring(R, "y" => 1:n, "t" => 1:n)
  # replace f(x_i) -> f(y_i + t_i)
  F = evaluate(f, gens(Rtemp)[1:n] + gens(Rtemp)[n+1:2n])
  # i = 1 # counter to iterate though degrees of monomials
  HasseDerivativesList = empty([f])
  varR = vcat(gens(R), ones(typeof(base_ring(R)(1)), n))
  # getting hasse derivs without extra attention on ordering
  for term in terms(F)
    # hasse derivatives are the factors in front of the monomial in t
    push!(HasseDerivativesList, evaluate(term, varR))
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

function hasse_derivatives(v::Vector)
  return hasse_derivatives.(v)
end




################################################################################
### internal functions for expert use

# MPolyQuoRingElem (internal, expert use only)
function _hasse_derivatives(f::MPolyQuoRingElem)
  return hasse_derivatives(lifted_numerator(f)) 
end

# Oscar.MPolyLocRingElem (internal, expert use only)
function _hasse_derivatives(f::Oscar.MPolyLocRingElem)
  return hasse_derivatives(numerator(f)) 
end

# Oscar.MPolyQuoLocRingElem (internal, expert use only)
function _hasse_derivatives(f::Oscar.MPolyQuoLocRingElem)
  return hasse_derivatives(lifted_numerator(f))
end

# for a list of elements (internal, expert use only)
function _hasse_derivatives(v::Vector)
  return _hasse_derivatives.(v)
end