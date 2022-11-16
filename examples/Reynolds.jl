module Reynolds

using Oscar

import Oscar: Hecke, Nemo, AbstractAlgebra
import Base:^

export ContinuedFraction, continued_fraction, convergents, reynolds

function ^(f::fmpq_mpoly, m::fmpq_mat)
  @assert nrows(m) == ncols(m) == ngens(parent(f))
  g = gens(parent(f))
  h = typeof(f)[]
  for i=1:length(g)
    push!(h, sum(m[i,j]*g[j] for j=1:length(g)))
  end
  z = evaluate(f, h)
  return z
end

function is_approx_integral(f::fmpq_mpoly, o::fmpz)
  g = MPolyBuildCtx(parent(f))
  for (c, v) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
    r = c*o
    if abs(round(r)-r) < 1//100
      push_term!(g, round(r)//o, v)
    else
      @show v, continued_fraction(c), 
      return false, finish(g)
    end
  end
  return true, finish(g)
end

function is_approx_integral(f::fmpq_mat, o::fmpz)
  g = 0*f
  for i = 1:nrows(f)
    for j = i:ncols(f)
      c = f[i, j]
      r = c*o
      if abs(round(r)-r) < 1//100
        g[i, j] = round(r)//o
      else
        @show v, continued_fraction(c), 
        return false, g
      end
    end
  end
  return true, g
end

struct ContinuedFraction
  cf::Vector{fmpz}
  function ContinuedFraction(a::fmpq, n::Int = typemax(Int))
    iszero(a) && return new(fmpz[0])
    cf = fmpz[]
    while true
      i = floor(fmpz, a)
      push!(cf, i)
      length(cf) >= n && return new(cf)
      a = a-i
      iszero(a) && return new(cf)
      a = inv(a)
    end
  end

  function ContinuedFraction(a::BigFloat, n::Int)
    iszero(a) && return new(fmpz[0])
    cf = fmpz[]
    while true
      i = floor(fmpz, a)
      push!(cf, i)
      n -= 1
      n == 0 && return new(cf)
      a = a-BigFloat(i)
      iszero(a) && return new(cf)
      a = inv(a)
    end
  end

  function ContinuedFraction(a::Vector{fmpz})
    return new(a)
  end

end

function Base.show(io::IO, C::ContinuedFraction)
  print(io, "<", C.cf[1])
  if length(C.cf) > 1
    print(io, ";", C.cf[2])
    for i=3:length(C.cf)
      print(io, ", ", C.cf[i])
    end
  end
  print(io, ">")
end

continued_fraction(a::fmpq) = ContinuedFraction(a)
continued_fraction(a::Rational) = ContinuedFraction(fmpq(a))

function best_approximation(a::fmpq, B::fmpz)
 C = continued_fraction(a)
 return convergents(C, limit = B)[end]
end

function best_approximation(a::fmpq_mpoly, B::fmpz)
  b = MPolyBuildCtx(parent(a))
  for (c, v) in zip(AbstractAlgebra.coefficients(a), AbstractAlgebra.exponent_vectors(a))
    push_term!(b, best_approximation(c, B), v)
  end
  return finish(b)
end

function convergents(C::ContinuedFraction; limit::fmpz = fmpz(-1), limit_terms::Int = -1)
  h = fmpz[C.cf[1]]
  k = fmpz[1]

  length(C.cf) == 1 && return [fmpq(h[i], k[i]) for i=1:length(h)]

  if (limit > -1 && C.cf[2] > limit) ||
    (limit_terms > -1 && length(h) >= limit_terms)
    return [fmpq(h[i], k[i]) for i=1:length(h)]
  end

  push!(h, C.cf[1]*C.cf[2]+1)
  push!(k, C.cf[2])
  
  for i=3:length(C.cf)
    if limit_terms > -1 && length(h) >= limit_terms
      break
    end
    nk = C.cf[i] * k[end]+k[end-1]
    if limit > -1 && nk > limit
      break
    end
    nh = C.cf[i] * h[end]+h[end-1]
    push!(h, nh)
    push!(k, nk)
  end
  return [fmpq(h[i], k[i]) for i=1:length(h)]
end

function bound_precision!(f::fmpq_mpoly, B::fmpz)
  for i=1:length(f)
    c = coeff(f, i)
    d = best_approximation(c, B)
    d = round(c*B)//B
    @assert !iszero(d)
    setcoeff!(f, i, d)
  end
end

#z, zi = make_integral...
#then [x*x*zi for x = M] is an integral rep
function make_integral(M::Vector{fmpq_mat})
  I = identity_matrix(ZZ, nrows(M[1]))
  d = [lcm(collect(map(denominator, m))) for m = M]
  N = [map(numerator, d[i]*M[i]) for i=1:length(M)]
  n = nrows(M[1])

  while true
    J = I
    for i=1:length(M)
      J = vcat(J*d[i], J*N[i])
      H = hnf(J)[1:n, :]
      c = content(H)
      J = Nemo.divexact(H, c)
    end
    if I == J
      I = lll(I)
      z = map_entries(QQ, I)
      zi = inv(z)
      return z, zi
    end
    I = J
  end
end

function reynolds(f::fmpq_mpoly, M::Vector{<:MatElem}, ord::fmpz)

  ord *= lcm(map(denominator, M))

  I = one(parent(M[1]))
  MM = []
  for m = M
    push!(MM, [I, m])
  end
  i = 1
  fl = true
  d = fmpz(1)
  while true
    i += 1
    f_new = f
    for J in MM
      d *= length(J)
      f_new = sum([f_new^m for m = J])
    end
#    g = best_approximation(f_new, ord)
    fl, g = is_approx_integral(f_new*fmpq(1, d), ord)
    if fl
      if all(x->g^x == g, M)
        @show "used ", i, " iterations"
        return g
      end
    else
#      @show "not approx"
    end
    f = f_new
#    bound_precision!(f, ord)
    if i>665 error("no conv"); end #not correct, but a failsafe
  end
  return f
end

#TODO not optimal at all - but  works..
function Oscar.monomials(R::MPolyRing, d::Int)
  S = Set(1:ngens(R))
  m = eltype(R)[]
  @assert d >0
  if d == 1
    return Set(gens(R))
  end
  n = monomials(R, d-1)
  nn = gen(R, 1) .* n
  for g = gens(R)[2:end]
    union!(nn, g .* n)
  end

  return nn
end

function action_on_monomials(R::MPolyRing, d::Int, A::Vector{fmpq_mat})
  g = gens(R)
  m = sort(monomials(R, d))

  B = fmpq_mat[]
  for a = A
    bb = Vector{fmpq}[]
    h = [x^a for x = g]
    for x = m
      b = zeros(QQ, length(m))
      y = evaluate(x, h)
      for (c,v) = zip(AbstractAlgebra.coefficients(y), AbstractAlgebra.monomials(y))
        b[findfirst(l -> l==v, m)] = c
      end
      push!(bb, b)
    end
    push!(B, matrix(QQ, bb))
  end
  return B, m
end

function reynolds(f::fmpq_mat, M::Vector{<:MatElem}, ord::fmpz)
  ord *= lcm(map(denominator, M))

  I = one(parent(M[1]))
  MM = []
  for m = M
    push!(MM, [I, m])
  end
  i = 1
  fl = true
  d = fmpz(1)
  while true
    i += 1
    f_new = f
    for J in MM
      d *= length(J)
      f_new = sum([f_new*m for m = J])
    end
#    g = best_approximation(f_new, ord)
    fl, g = is_approx_integral(f_new*fmpq(1, d), ord)
    if fl
      if all(x->g*x == g, M)
        @show "used ", i, " iterations"
        return g
      end
    else
#      @show "not approx"
    end
    f = f_new
#    bound_precision!(f, ord)
    if i>665 error("no conv"); end #not correct, but a failsafe
  end
  return f
end

function reynolds_via_linalg(f::fmpq_mpoly, A::Vector{fmpq_mat}, ord::fmpz)
  @assert length(f) == 1 #currently technical probs

  h, m = action_on_monomials(parent(f), total_degree(f), A)
  F = zero_matrix(QQ, 1, length(m))
  F[1, findfirst(x->x == monomial(f, 1), m)] = coeff(f, 1)
  r = reynolds(F, h, ord)

  return sum([r[1, i] * m[i] for i=1:length(m)])
end


end #module

using .Reynolds


