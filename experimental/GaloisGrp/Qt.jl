#= musings about Qt - before I forget
 
 f in Q(t)[x] - or better Z[t][x]

 for all t (outside some bad points), the roots are power series over C
 
   R_i(z) =: R(z) = sum a_n (z-t)^n

 By Taylor, Cauchy, Dan and his partner, 
 
 a_n = 1/2/pi/i int_{|z| = r} R(z)/(z-t)^(n+1) dz

 Now f(z, R(z)) = 0

 f = sum f_i(t) x^i

 By standart bounds on polynomials (assume f monic)

 roots of f(z)(x) are bounded in abs. value by |f_i(z)| + 1 
 (or even 2*|f_i(z)|^(1/(n-i)) or so, 
 https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots
 Donald Knuth)

 if |f|_infty = B (largest coeff in Z), then
   |f_i(z)| <= degree(f_i) B |z|^degree(f_i)

 so |R(z)| has an explicit upper bound, growing no worse than |z|^degree(f)

 so 

 |a_n| <= 1/2/pi (2 pi r) deg(f) B (|t|+r)^deg(f)/r^(n+1)
        = deg(f) B (|t|+r)^deg(f)/ r^n

 This is valid for all r where R is analytic.

 Let d = disc(f), then R should be analytic if r < min(|t-s| for d(s) = 0)
 
 if this is > 1, the a_n -> 0    

=======================================================

OK, s is a root, so the coeffs s = sum s_n (x-t)^n 
satisfy, by above, |s_n| <= B/r^n for some explicit r > 1

Lemma:
  let s, t be power series with
    |s_n| <= B (n+1)^k/r^n for some B, k
    |t_n| <= C (n+1)^l/r^n
  then
    st = sum d_i (x-t)^n with
    |d_n| <= BC (n+1)^(k+l+1)/r^n

pf:
  Cauchy product:

  |d_n| =  |sum_{i+j = n} s_i t_j|
        <= sum |s_i t_j|
        <= sum Bi^k/r^i C j^l/r^j
        <= BC sum (n+1)^k (n+1)^l/r^(i+j) = BC (n+1)^(k+l) (n+1) 1/r^n

Clear: 
  |(s+t)_n| <= ...


Lemma:
  if |s_n| <= B (n+1)^k/r^n, then the largest coeff is either at
    floor(x) or ceil(x) for x = (k/log r) - 1

pf:
 as a function of n:
   (B (n+1)^k/r^n)' = B k (n+1)^(k-1)/r^n - B (n+1)^k log r /r^n
 which is zero iff
   k - (n+1) log r = 0 or n+1 = k/log r


Possibly we can get better estimates for sum_{i+j=n} i^k j^l
on particular for powering:
  sum_{|alpha| = k} prod alpha^k

===============================================================
Also see 
https://math.uni-paderborn.de/fileadmin/mathematik/AG-Computeralgebra/Publications-klueners/function.pdf

for analysis of the denominator and the infinite valuations  
=#


function Hecke.lcm(a::Vector{<:RingElem})
  if length(a) == 0
    error("don't know the ring")
  end
  return reduce(lcm, a)
end

function Hecke.subfields(F::Generic.FunctionField{fmpq})
  f = defining_polynomial(F)
  @assert ismonic(f)
  Zxy, (x, y) = PolynomialRing(FlintZZ, 2, cached = false)
  ff = Zxy()
  d = lcm(map(denominator, coefficients(f)))
  df = f*d

  dd = lcm(map(denominator, coefficients(d)))
  dd = lcm(dd, lcm(map(x->reduce(lcm, map(denominator, coefficients(numerator(x))), init=fmpz(1)), coefficients(df))))
  @assert isone(dd) #needs fixing....
  df *= dd

  for i=0:degree(f)
    c = coeff(df, i)
    if !iszero(c)
      ff += map_coefficients(FlintZZ, numerator(c))(y)*x^i
    end
  end
  _subfields(F, ff, (10, 4))

end


function _subfields(FF::Generic.FunctionField, f::fmpz_mpoly, prec::Tuple{Int, Int} = (2,2))
  Zxt = parent(f)
  X,T = gens(Zxt)

  Zx = Hecke.Globals.Zx

  @vprint :Subfields, 1, "Starting subfield computation...\n"
  @vprint :Subfields, 2, "for $f\n"

  t = -1
  local g::fmpz_poly
  while true
    t += 1
    g = evaluate(f, [gen(Zx), Zx(t)])
    if isirreducible(g)
      break
    end
    if t > 10
      error("cannot be")
    end
  end

  @vprint :Subfields, 2, "substituting t = $t\n"


  @vprint :Subfields, 1, "Computing subfields of the number field...\n"

  K, a = number_field(g, cached = false)
  @vtime :Subfields  1  S = subfields(K)

  @vprint :Subfields, 2, "now looking for a nice prime...\n"
  p, _ = find_prime(defining_polynomial(K))

  d = lcm(map(degree, collect(keys(factor(g, GF(p)).fac))))

  ff = evaluate(f, [X, T+t]) #do the shift...
  @assert evaluate(ff, [gen(Zx), zero(Zx)]) == g

  @vprint :Subfields  1  "using prime $p and degree $d\n"
  C = GaloisCtx(f, p, d)

  @vprint :Subfields, 2, "obtaining roots with prec $prec\n"
  @vtime :Subfields  2  R = roots(C, prec)
  
  F = parent(R[1]) # should be Qq<<t>>
  Qq, mQq = ResidueField(F)
  Fq, mFq = ResidueField(Qq)

  rc = map(ComposedFunction(mFq, mQq), R)

  res = []

  for (k, mkK) = S
    @vprint :Subfields  1  "processing $k\n"
    if degree(k) == 1 || degree(k) == degree(K)
      @vprint :Subfields 2 "dropping trivial field\n"
      continue
    end

    a = mkK(gen(k))
    A = parent(K.pol)(a)
    @vprint :Subfields, 2, "embedding polynomial: $A\n"
    ts = gen(Zx)
    local bs
    r = copy(rc)
    while true
      b = map(map_coefficients(Fq, A), r)
      bs = Hecke.MPolyFact.block_system(b)
      if length(bs) == degree(k) && all(x->length(x) == length(bs[1]), bs)
        break
      end
      @vprint :Subfields 2 "need tschirni\n"
      ts = rand(Zx, 2:length(r), -4:4)
      r = map(ts, rc)
    end
    @vprint :Subfields 2 "... with block system $bs\n"

    @vprint :Subfields 2 "trying sum as primitive element...\n"
    #TODO: test over finite field first and make sure p is large enough
    con = [sum(R[b]) for b = bs]
 
    if length(Set(con)) < length(bs)
      @vprint :Subfields 2 "...failed, produt next...\n"
      con = [prod(R[b]) for b = bs]
      while length(Set(con)) < length(bs)
        @vprint :Subfields 2 "...failed, adding 1...\n"
        R .+= 1
        con = [prod(R[b]) for b = bs]
      end
    end

    @vprint :Subfields 2 "building power sums (traces)\n"
    pow = copy(con)
    @assert length(Set(pow)) == length(bs)
    fl, tt = isinteger(C, prec, sum(pow))
    if !fl
      error("prec too low (or no subfield")
    end
    tr = [tt]
    while length(tr) < degree(k)
      pow .*= con
      fl, tt = isinteger(C, prec, sum(pow))
      if !fl
        error("prec too low (or no subfield")
      end
      push!(tr, tt)
    end

    Qt = base_ring(FF)
    _t = gen(Qt)
    @vprint :Subfields 2 "Newton relations for the polynomial\n"
    @vtime :Subfields 2 ps = power_sums_to_polynomial(map(x->x(_t-t), tr))  #should be the subfield polynomial

    @vprint :Subfields 2 "Interpolate for the embedding\n"
    @vtime :Subfields 2 f = interpolate(PolynomialRing(F)[1], R, [con[findfirst(x->i in x, bs)] for i=1:length(R)])   # should be the embedding poly
    #if I read Florian's thesis correct, then f should be Z-integral in the
    #Kronnecker presentation
    # it is given here as a poly over the function field
    # step 1: poly over poly{Q}
    fff = map_coefficients(x->numerator(x)(gen(F)), defining_polynomial(FF), parent = parent(f))
    # step 2: move to Kronnecker, so that we get poly over poly{Z}
    f = (f*derivative(fff)) % fff
    # step 3: lift (and put the denominator back in)
    emb = map_coefficients(x->isinteger(C, prec, x)[2](gen(Nemo.base_ring(FF))), f, parent = parent(defining_polynomial(FF)))(gen(FF)) // derivative(defining_polynomial(FF))(gen(FF))

    @assert iszero(ps(emb)) #final verification
    ps = map_coefficients(x->parent(x)(Generic._rat_canonicalise(numerator(x), denominator(x))...), ps, parent = parent(ps))
    @vprint :Subfields 1 "have subfield, defined by $ps with embedding $emb\n"
    push!(res, (FunctionField(ps, "a", cached = false)[1], emb))
  end
  return res
end

Oscar.gen(R::Generic.RationalFunctionField) = R(gen(base_ring(R.fraction_field)))

function isinteger(G::GaloisCtx, p::Tuple{Int, Int}, r::Generic.RelSeries{qadic})
#  @show "testing", r, "against", p
  if iszero(r) 
    return true, Hecke.Globals.Qx(0)
  end
  if r.length + r.val > p[2]
    return false, Hecke.Globals.Qx(0)
  end
  f = Hecke.Globals.Qx(0)
  x = gen(parent(f))
  xpow = parent(x)(1)
  for i = 0:(r.length + r.val - 1)
    c = coeff(r, i)
    pr = prime(parent(c))
    if c.length < 2
      cc = coeff(c, 0)
      f += xpow* Hecke.mod_sym(lift(cc), pr^p[1])
    else
      return false, x
    end
    xpow *= x
  end
  return true, f
end

