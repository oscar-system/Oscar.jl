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
  _subfields(F, ff)

end


function _subfields(FF::Generic.FunctionField, f::fmpz_mpoly)
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

  @vprint :Subfields, 2, "obtaining roots with minimal prec $prec\n"
  @vtime :Subfields  2  R = roots(C, (2, 1))
  
  F = parent(R[1]) # should be Qq<<t>>
  Qq, mQq = ResidueField(F)
  Fq, mFq = ResidueField(Qq)

  rc = map(ComposedFunction(mFq, mQq), R)

  res = []

  SL = SLPolyRing(ZZ, degree(K))
  sx = gens(SL)

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
    sx = gens(SL)
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
    con = [sum(r[b]) for b = bs]
    sx = map(ts, sx)
    conI = [sum(sx[b]) for b = bs]
 
    if length(Set(con)) < length(bs)
      @vprint :Subfields 2 "...failed, product next...\n"
      con = [prod(r[b]) for b = bs]
      conI = [prod(sx[b]) for b = bs]
      while length(Set(con)) < length(bs)
        @vprint :Subfields 2 "...failed, adding 1...\n"
        r .+= 1
        sx .+= 1
        con = [prod(r[b]) for b = bs]
        conI = [prod(sx[b]) for b = bs]
      end
    end

    @vprint :Subfields 2 "now proper bounds (for subfield poly)\n"
    B = upper_bound(C, power_sum, conI, length(conI))
    @vprint :Subfields 2 "coeffs $B\n"
    prec_poly = bound_to_precision(C, B)
    @vprint :Subfields 2 "gives a precision of $prec_poly\n"

    function get_poly(prec::Tuple{Int, Int})
      R = roots(C, prec)
      con = [evaluate(c, R) for c = conI]


      @vprint :Subfields 2 "building power sums (traces)\n"
      pow = copy(con)
      @assert length(Set(pow)) == length(bs)
      fl, tt = isinteger(C, prec, sum(pow))
      if !fl
        return nothing
      end
      tr = [tt]
      while length(tr) < degree(k)
        pow .*= con
        fl, tt = isinteger(C, prec, sum(pow))
        if !fl
          return nothing
        end
        push!(tr, tt)
      end

      Qt = base_ring(FF)
      _t = gen(Qt)
      @vprint :Subfields 2 "Newton relations for the polynomial\n"
      local ps
      try
        @vtime :Subfields 2 ps = power_sums_to_polynomial(map(x->x(_t-t), tr))  #should be the subfield polynomial
      catch e
        @show e
        return nothing
      end
      ps = map_coefficients(x->parent(x)(Generic._rat_canonicalise(numerator(x), denominator(x))...), ps, parent = parent(ps))
      return ps
    end

    @vprint :Subfields 2 "Bounds for the embedding\n"
    #=we need bounds, again...
    #more fun:
    we have gamma: primitive element in large field F
            alpha: primitive element  in subfield

    all monic, so alpha is integral, alpha in IntCls(Z[t], F)
    IntCls(Z[t], F) <= 1/f' Z[t][x]/f (via Florian)

    Markus' (or Script)
    let g = f(x)/x-gamma = sum g_i x^i for some g_i in F, then
        Tr(g_i/f' * gamma^j) = delta_i,j, ie we have the trace-dual basis
    
    we want alpha = sum a_i gamma^i    for some a_i in Q(t), moving Kronecker:
            alpha = sum b_i gamma^i/f' now some b_i in Z[t]
    
       mult by g_j and taking traces:
          Tr(alpha * g_j) = sum b_i Tr(g_j * gamma^i/f') = b_j

    We need bounds for b_j, we'll find b_j through interpolation...
    so need bound for g_j:
      The coeffs of f (and g) are elem. symm. in the roots. For f in all n roots
      for g in n-1, wlog. 1..n-1
      We have (rt_n <= B, B is a bound for all roots, not just the n-th)
        sigma_n,1 = sigma_n-1,1 + rt_n, so
        sigma_n-1,1 <= sigma_n,1 + B

        sigma_n,2 = sigma_n-1,2 + rt_n * sigma_n-1,1  so
        sigma_n-1,2 <= sigma_n, 2 + B sigma_n-1,1 <= sigma_n,2 + B sigma_n,1 + B^2
        ...
        sigma_n,n-1 = sigma_n-1,n-1 + rt_n * sigma_n-1,n-2
        ...
        sigma_n-1,n-1 <= sigma_n,n-1 + B sigma_n,n-2 + ... + B^n-2 sigma_n,1 + B^n-1

        so, basically
          g_i <= div(f, x)(B)
        (one can also see it by doing the division by "hand")  
        we have (bound for) alpha (the conjugates), for g_j (above), then *n are
        bounds for b_j

        On 2nd thougts I might also just compute the g_i via division in F[x]  

        Darn: need to think again: need estimates for poly * series.
    =#
    c = coefficients(f, 1)
    # start with the degree...
    d = C.B.val[3]
    dd = fmpq(0)
    for i=2:length(c)
      if !iszero(c[i])
        dd = max(dd, degree(c[i], 2)+(i-1)*d)
      end
    end

    #power series * poly:
    #we have a degree bound, so I can use poly * poly instead
    #there we may use Beauzamy
    #n(f) = (sum 1/binom(d, j) |f_j|^2)^(1/2)
    #this is sub multiplicative: n(fg) <= n(f) n(g)
    #also: f_j <= binom(d,j) n(f)

    # poly -> power series? need B, k s.th |a_i| <= B (i+1)^k/r^i
    # |a_i|/r^i <= B (i+1)^k
    # maybe: k = average/ceil/floor log_(i+1) a_i/r^i = log_(i+1) a_i - i log_(i+1) r
    # and then B as the minimum that works?
    #
    # Easy solution for now:poly to series: f -> (deg(f)+1)) max|f_i|

    B = parent(C.B)(0)
    for i=length(c):-1:2
      if !iszero(c[i])
        B += parent(B)(evaluate(c[i], [gen(Hecke.Globals.Zx), Hecke.Globals.Zx(0)]))
      end
      B *= C.B
    end

    @vprint :Subfields 2 "coeff of embedding $B\n"
    prec_emb = bound_to_precision(C, B)
    @vprint :Subfields 2 "precision of $prec_emb\n"

    function get_emb(prec::Tuple{Int, Int})
      R = roots(C, prec)
      con = [evaluate(c, R) for c = conI]

      local ff
      try
        @vtime :Subfields 2 ff = interpolate(PolynomialRing(F)[1], R, [con[findfirst(x->i in x, bs)] for i=1:length(R)])   # should be the embedding poly
      catch e
        @show e
        return nothing
      end
      #if I read Florian's thesis correct, then f should be Z-integral in the
      #Kronnecker presentation
      # it is given here as a poly over the function field
      # step 1: poly over poly{Q}
      fff = map_coefficients(x->numerator(x)(gen(F)), defining_polynomial(FF), parent = parent(ff))
      # step 2: move to Kronnecker, so that we get poly over poly{Z}
      ff = (ff*derivative(fff)) % fff
      # step 3: lift (and put the denominator back in)
      emb = map_coefficients(x->isinteger(C, prec, x)[2](gen(Nemo.base_ring(FF))), ff, parent = parent(defining_polynomial(FF)))(gen(FF)) // derivative(defining_polynomial(FF))(gen(FF))
      return emb
    end

    pr = (1,1)
    up = (max(prec_poly[1], prec_emb[1]), max(prec_poly[2], prec_emb[2]))

    local ps, emb

    while true
      pr = (min(pr[1], up[1]), min(pr[2], up[2]))
      @vprint :Subfields 1 "trying to compute data with precision $pr\n"
      ps = get_poly(pr)
      if ps === nothing
        if pr[1] >= prec_poly[1] && pr[2] >= prec_poly[2]
          @vprint :Subfields 2 "block system did not give polynomial\n"
          break
          error("not a subfield")
        end
        pr = (2*pr[1], 2*pr[2])
        continue
      end
      emb = get_emb(pr)
      if emb === nothing || !iszero(ps(emb)) 
        if pr[1] >= prec_emb[1] && pr[2] >= prec_emb[2]
          @vprint :Subfields 2 "block system did not give (working) embedding\n"
          break
          error("not a subfield")
        end
        pr = (2*pr[1], 2*pr[2])
        continue
      end
      break
    end

    if ps === nothing || emb === nothing
      continue
    end

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

function newton_polygon(f::Generic.Poly{Generic.Rat{T}}) where {T}
  pt = Tuple{Int, Int}[]
  for i=0:degree(f)
    c = coeff(f, i)
    if !iszero(c)
      push!(pt, (i, degree(numerator(c))-degree(denominator(c))))
    end
  end
  return Hecke.lower_convex_hull(pt)
end

function valuations_of_roots(f::Generic.Poly{Generic.Rat{T}}) where {T}
  return [(-slope(l), length(l)) for l = Hecke.lines(newton_polygon(f))]
end

Hecke.lines(P::Hecke.Polygon) = P.lines
slope(l::Hecke.Line) = l.slope


