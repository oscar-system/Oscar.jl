
###############################################################################
#
#  Cyclotomic polynomials
#
###############################################################################

function _cyclotomic_polynomial(n::Int64)
  @assert n > 0
  _, x = QQ["x"]
  return Hecke.cyclotomic(n, x)
end

function _is_cyclotomic_polynomial(p::Union{fmpz_poly, fmpq_poly})
  n = degree(p)
  R = parent(p)
  x = gen(R)
  list_cyc = union(Int64[k for k in euler_phi_inv(n)], [1])
  list_poly = [Hecke.cyclotomic(k, x) for k in list_cyc]
  return any(q -> R(collect(coefficients(q))) == p, list_poly)
end

###############################################################################
#
#  Exponent of fmpq/fmpz_mat
#
###############################################################################

function _is_of_finite_exponent(f::Union{fmpq_mat, fmpz_mat})
  !Hecke.is_squarefree(minpoly(f)) && return false
  chi = charpoly(f)
  fact = collect(factor(chi))
  return all(p -> _is_cyclotomic_polynomial(p[1]), fact)
end

function _exponent(f::Union{fmpq_mat, fmpz_mat})
  !_is_of_finite_exponent(f) && return -1
  degs = unique(degree.([p[1] for p in collect(factor(minpoly(f)))]))
  exps = euler_phi_inv(degs[1])::Vector{fmpz}
  for i in 2:length(degs)
    union!(exps, euler_phi_inv(degs[i]))
  end
  maxdeg = lcm(exps)
  divmd = divisors(maxdeg)
  n = findfirst(k -> isone(f^k), divmd)
  @assert n !== nothing
  return return divmd[n]
end
