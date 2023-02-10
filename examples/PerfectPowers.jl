module PerfectPowers

using Oscar
import Oscar: Nemo, Hecke
import Base: powermod

function root_exact(a::fmpz, n::Int)
  @assert n>0
  n == 1 && return a
  k, a = remove(a, fmpz(2))
  if k % n != 0
    return fmpz(1)
  end
  kn = div(k,n) #since we change n
  isone(a) && return fmpz(2)^kn

  while !isone(a) && iseven(n)
    a = _root_exact(a, Val(2))
    n = div(n, 2)
  end
  (isone(a) ||isone(n)) && return a
  a = _root_exact(a, n)
  isone(a) && return a
  return a*fmpz(2)^kn
end

#TODO
#  write the _root_exat in mpn
#  prealloc temp
#  understand why the top bit is wrong for 7^10 = 16807 and p = 2
function _root_exact(a::fmpz, p::Val{2})
  #return the p-th root, or garbage of the correct size...

  @assert isodd(a)

  n = nbits(a)
  s = div(n, 2)+1 # bits in p-th root

  A = Hecke.Limbs(a, MSW=false)
  d = A[1]
  if d % 8 != 1
    return fmpz(1)
  end
  dd = d

  b = d
  for i=1:6
    _d = d
    d = d*(1+(div((1-b*d*d), 2)))
  end
  if s < 63
    d = (d<<1)>>1  # top bit is wrong for d = 16807^2
    #there are 2 roots in this case ? and -? so we need to try both
    #since we do >> 1 (or div(, 2)) we loose the high bit.
    d *= dd
    nbits(d) <= s && return fmpz(d)
    d = -d #we could have found the negative one by accident...
    nbits(d) <= s && return fmpz(d)
    return fmpz(1)
  end

  B = mod(a, fmpz(2)^(s+3))
  D = fmpz(d)

  M = fmpz(2)^(64)
  i = 64 
  one = fmpz(1)
  #TODO; better chain of exponents, use mpn? more inline?
  #TODO: better mod 2^n, use mullow?
  while i < s
    Nemo.mul!(M, M, M)
    ccall((:fmpz_fdiv_q_2exp, Nemo.libflint), Nothing, (Ref{fmpz}, Ref{fmpz}, Int), M, M, 2)

    T = powermod(D, 2, M)
    Nemo.mul!(T, T, B)
    Hecke.mod!(T, T, M)
    Nemo.sub!(T, one, T)
    ccall((:fmpz_fdiv_q_2exp, Nemo.libflint), Nothing, (Ref{fmpz}, Ref{fmpz}, Int), T, T, 1)
    Nemo.add!(T, one, T)
    Nemo.mul!(D, D, T)
    Hecke.mod!(D, D, M)

    i = 2*i-2
  end
  M = fmpz(2)^(s+2)

  Nemo.mul!(D, D, a)
  Hecke.mod!(D, D, M)
  if nbits(D) <= s
    return fmpz(D)
  end
  
  Nemo.mul!(D, D, -1) # in case the iteration found the negative root
  Nemo.add!(D, D, M)
  if nbits(D) > s
    return fmpz(1)
  end

  return D
end

@inline function fmpz_trunc!(a::fmpz, i::Int)
#  @assert a >= 0
  ccall((:fmpz_fdiv_r_2exp, Nemo.libflint), Cvoid, (Ref{fmpz}, Ref{fmpz}, Clong), a, a, i)
end

function powermod_2exp(a::fmpz, p::Int, i::Int) #too slow - much worde than 
                                              #powermod directly
  @assert a > 0 && p > 0
  a = copy(a)
  while iseven(p)
    Nemo.mul!(a, a, a)
    fmpz_trunc!(a, i)
    p >>= 1
  end
  b = a
  p >>= 1
  while p>0
    while iseven(p)
      Nemo.mul!(a, a, a)
      fmpz_trunc!(a, i)
      p >>= 1
    end
    Nemo.mul!(a, a, a)
    fmpz_trunc!(a, i)
    Nemo.mul!(b, b, a)
    fmpz_trunc!(b, i)
    p >>= 1
  end

  return b
end

function _root_exact(a::fmpz, p::Int, extra_s::Int = 5, extra_w::Int = 1)
  #return the p-th root, or garbage of the correct size...

  @assert isodd(a)

  p==2 && return _root_exact(a, Val(2))

  if iseven(p)
    @show "slow", p
    fl, a = Hecke.is_power(a, p)
    return fl ? a : fmpz(1)
  end
  n = nbits(a)
  s = div(n, p)+1 # bits in p-th root
  w = div(s, 64)+1 #words in p-th root
 
  A = Hecke.Limbs(a)
  d = A[1]
  dd = d

  l = invmod(p, UInt(2^32))
  l = l*(2-l*p) # inv mod 2^64 - Julia cannot do this directly...

  b = d^(p-1)
  pr = 1
  for i=2:7
    d = (1+l)*d-l*b*d^(p+1)
    pr *= 2
    if pr > s + extra_s
      break
    end
  end
  if s + extra_s < 64
    d *= dd
    d = d % 2^min(pr, 63)
    iseven(d) && return fmpz(1)
    nbits(d) <= s && return fmpz(d)
    return fmpz(1)
  end

#  B = powermod_2exp(a, p-1, 64*w+1)
  B = powermod(a, p-1, fmpz(2)^(64*(w+extra_w)))
  L = fmpz(l)
  D = fmpz(d)

  M = fmpz(2)^(64)
  i = 1
  two = fmpz(2)
  #TODO; better chain of exponents, use mpn? more inline?
  #TODO: better mod 2^n (done), use mullow?
  #TODO: find powermod for 2^n modulus
  while i < w + extra_w
    i *= 2
    Nemo.mul!(M, M, M)
    T = L*p
    Nemo.sub!(T, two, T)
    Nemo.add!(T, T, M)  #make pos
    Nemo.mul!(L, L, T)
    _x = L % M
    fmpz_trunc!(L, i*64)
    @assert _x == L
#    Hecke.mod!(L, L, M)

#    T = powermod_2exp(D, p+1, i*64)
    T = powermod(D, p+1, M)
    Nemo.mul!(T, T, B)
    _x = T % M
    fmpz_trunc!(T, i*64)
    @assert _x == T
#    Hecke.mod!(T, T, M)

    Nemo.mul!(T, T, L)
    _x = T % M
    fmpz_trunc!(T, i*64)
    @assert _x == T

    S = 1+L
    Nemo.mul!(S, S, D)
    _x = S % M
    fmpz_trunc!(S, i*64)
    @assert S == _x
#    Hecke.mod!(S, S, M)
    Nemo.sub!(D, S, T)
    if sign(D) < 0
      Nemo.add!(D, D, M)
    end
  end

  Nemo.mul!(D, D, a)
  fmpz_trunc!(D, (w+extra_w)*64)
#  Hecke.mod!(D, D, M)
  if nbits(D) > s
    return fmpz(1)
  end

  return D
end

function is_power_bernstein(a::fmpz)
#https://cr.yp.to/lineartime/powers2-20060914-ams.pdf
  if isone(a)
    return (0, a)
  end
  k, a = remove(a, fmpz(2))
  if isone(a)
    return (k, fmpz(2))
  end
  l, a = remove(a, fmpz(3))
  g = gcd(k, l)
  if isone(a)
    # a is/was 2^k * 3^l
    return (g, fmpz(2)^divexact(k, g)*fmpz(3)^divexact(l, g))
  end
  h, a = remove(a, fmpz(5))
  g = gcd(h, g)
  if isone(a)
    return (g, fmpz(2)^divexact(k, g)*fmpz(3)^divexact(l, g)*fmpz(5)^divexact(h, g))
  end


  p_test = next_prime(2^40)
  a_test = a % p_test

  c = [a]
  cl = clog(a, 7)
  if k>0
    cl = min(cl, k)
  end
  no_p = 0
  no_s = 0

  f = 1

  for p = PrimesSet(2, cl)
    if k % p != 0
      continue
    end
    if p > cl
      break
    end
    if p>2 && p < 2^15
      aa = UInt(a % p^2)
      if aa %p != 0 && powermod(aa, p-1, UInt(p^2)) != 1
        no_s += 1
        continue
      end
    end
    no_p += 1
    pp = p
    d = _root_exact(a, p, 10, 2)
    if !isone(d)
      if powermod(d, p, fmpz(p_test)) == a_test && d^p == a
        f *= p
        cl = min(cl, flog(d, 7))
        a = d
        a_test = a % p_test
      else
         d = fmpz(1)
#        @show p, powermod(d, p, fmpz(p_test)) == a % p_test
      end
    end
    while !isone(d) && pp <= cl
      push!(c, d)
      pp *= p
      d = root_exact(d, p)
      if !isone(d)
        if d^p== a
          f *= p
          cl = min(cl, flog(d, 7))
          a = d
          a_test = a % p_test
        else
          @show 1, p
        end
      end
    end
  end
  @show no_p, no_s
#  @assert c[1] == a
  
  c = collect(Set(c))
  @time c = [gcd(a, x) for x = c] #a should be a proper divisor of c[1]
  c = [x for x = c if !isone(x)]
  #not sure still necessary
  @time c = Hecke.coprime_base(c) 
  k = [valuation(a, p) for p = c]
  @show no_p
  return gcd(k)*f
end

export ispower_exact, root_exact

end

using .PerfectPowers

