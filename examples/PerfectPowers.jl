module PerfectPowers

using Oscar

import Oscar: Nemo

mutable struct fmpz_limbs
  d::Int
  s::Int
  l::Ptr{UInt}
  
  function fmpz_limbs(a::fmpz)
    r = new()
    r.d = a.d
    if Nemo._fmpz_is_small(a)
      r.s = -1
      return r
    end
    r.s = Nemo.__fmpz_limbs(a.d)
    @assert r.s > 0 #sign and length
    r.l = unsafe_load(convert(Ptr{Ptr{UInt}}, unsigned(a.d) << 2) + 2*sizeof(Cint))
    return r
  end
end

function Base.iterate(a::fmpz_limbs, i::Int = 0)
  if a.s == -1
    i == 0 && return UInt(a.d), 1
    return nothing
  end
  if i< a.s
    return unsafe_load(a.l, i+1), i+1
  end
  return nothing
end

Base.length(a::fmpz_limbs) = a.s < 0 ? 1 : a.s

function ispower_exact(a::fmpz, p::Int)
  #return the p-th root, or garbage of the correct size...

  @assert isodd(a)

  if iseven(p)
    fl, a = Hecke.ispower(a, p)
    return fl ? a : fmpz(1)
  end
  n = nbits(a)
  s = div(n, p)+1 # bits in p-th root
  w = div(s, 64)+1 #words in p-th root

  A = fmpz_limbs(a)
  d = iterate(A)[1]
  dd = d

  l = invmod(p, UInt(2^32))
  l = l*(2-l*p) # inv mod 2^64 - Julia cannot do this directly...

  b = d^(p-1)
  for i=2:7
    d = (1+l)*d-l*b*d^(p+1)
  end
  if s < 64
    return d*dd
  end

  B = powmod(a, p-1, fmpz(2)^(64*w+1))
  L = fmpz(l)
  D = fmpz(d)

  M = fmpz(2)^(64)
  i = 1
  two = fmpz(2)
  #TODO; better chain of exponents, use mpn? more inline?
  #TODO: better mod 2^n, use mullow?
  while i < w
    Nemo.mul!(M, M, M)
    T = L*p
    Nemo.sub!(T, two, T)
    Nemo.mul!(L, L, T)
    Hecke.mod!(L, L, M)

    T = powmod(D, p+1, M)
    Nemo.mul!(T, T, B)
    Hecke.mod!(T, T, M)

    Nemo.mul!(T, T, L)
    Hecke.mod!(T, T, M)
    S = 1+L
    Nemo.mul!(S, S, D)
    Hecke.mod!(S, S, M)
    Nemo.sub!(D, S, T)
    i *= 2
  end

  Nemo.mul!(D, D, a)
  Hecke.mod!(D, D, M)
  if nbits(D) > s
    return fmpz(1)
  end

  return D
end

function ispower_bernstein(a::fmpz)
#https://cr.yp.to/lineartime/powers2-20060914-ams.pdf
  k, a = remove(a, fmpz(2))

  c = [a]
  cl = clog(a, 2)
  if k>0
    cl = min(cl, k)
  end
  for p = PrimesSet(2, cl)
    if k % p != 0
      continue
    end
    pp = p
    while pp <= cl
      push!(c, ispower_exact(a, pp))
      isone(c[end]) && break
      pp *= p
    end
  end

  c = coprime_base(c)
  k = [valuation(a, p) for p = c]
  return gcd(k)
end

export ispower_exact

end

using .PerfectPowers

export ispower_exact
