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


function _galois_init(F::Generic.FunctionField{fmpq}; tStart::Int = -1)
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
  _subfields(F, ff, tStart = tStart)
end


function _subfields(FF::Generic.FunctionField, f::fmpz_mpoly; tStart::Int = -1)
  Zxt = parent(f)
  X,T = gens(Zxt)

  Zx = Hecke.Globals.Zx

  @vprint :Subfields, 1, "Starting subfield computation...\n"
  @vprint :Subfields, 2, "for $f\n"

  d = numerator(discriminant(FF))
  rt = roots(d, ComplexField(20))
  t = tStart
  local g::fmpz_poly
  while true
    t += 1
    g = evaluate(f, [gen(Zx), Zx(t)])
    if Hecke.lower_bound(minimum([abs(x-t) for x = rt]), fmpz) >= 2 && isirreducible(g)
      break
    end
    if t > 10
      error("cannot be")
    end
  end

  @vprint :Subfields, 2, "substituting t = $t\n"


  K, a = number_field(g, cached = false)

  @vprint :Subfields, 2, "now looking for a nice prime...\n"
  p, _ = find_prime(defining_polynomial(K), pStart = 200)

  d = lcm(map(degree, collect(keys(factor(g, GF(p)).fac))))

  @assert evaluate(evaluate(f, [X, T+t]), [gen(Zx), zero(Zx)]) == g

  @vprint :Subfields  1  "using prime $p and degree $d\n"
  C = GaloisCtx(f, t, p, d)
  return C, K, p
end

function galois_group(FF::Generic.FunctionField{fmpq}; overC::Bool = false)
  tStart = -1
  tr = -1
  while true
    tr += 1
#    @show tr, tStart
    if tr > 15
      error("not plausible")
    end
    C, K, p = _galois_init(FF, tStart = tStart)

    if issquare(discriminant(K)) != issquare(discriminant(FF))
      @vprint :GaloisGroup 2 "bad evaluation point: parity changed\n"
      tStart += 1
      continue
    end

    f = C.C.f
    tStart = Int(C.data[2]) # fragile...
    @vprint :GaloisGroup 1 "specialising at t = $tStart, computing over Q\n"
    Gal, S = galois_group(K, prime = p)

    @vprint :GaloisGroup 1 "after specialisation, group is: $(transitive_identification(Gal))\n"

#    @show S.start, S.chn

    #need to map the ordering of the roots in C to the one in S
    #the roots in C should be power-series over the field used in S
    rC = roots(C, (1,1))
    rS = roots(S, 1)

    F, mF = ResidueField(parent(rC[1]))
    G, mG = ResidueField(F)
    Qt_to_G = x->G(numerator(x)(tStart))//G(denominator(x)(tStart))

    H, mH = ResidueField(parent(rS[1]))

    mp = embed(H, G)
    _rC = map(x->mG(mF(x)), rC)
    _rS = map(x->mp(mH(x)), rS)
    @assert Set(_rC) == Set(_rS)
    prm = [findfirst(x->x == y, _rS) for y = _rC]
    pr = symmetric_group(length(prm))(prm)

    # need to verify the starting group
    #either, if block-systems are there, find the subfields
    #or verify the factorisation of the resolvent:
    if isdefined(S, :start) && S.start[1] == 1 && length(S.start) > 0
      if !all(x->issubfield(FF, C, [[inv(pr)(i) for i = y] for y = x]) !== nothing, S.start[2])
        @vprint :GaloisGroup 2 "bad evaluation point: subfields don't exist\n"
        continue
      end
      C.start = (1, [[[inv(pr)(i) for i = y] for y = x] for x = S.start[2]])
    elseif isdefined(S, :start)
      @assert S.start[1] == 2
        
      if length(S.start[2]) == 1 #msum was irreducible over Q, so also over Qt
        @vprint :GaloisGroup 1 "msum irreducible over Q, starting with Sn/An\n"
      else
        O = sum_orbits(FF, Qt_to_G, map(x->mG(mF(x)), rC))
        if sort(map(length, O)) != sort(map(length, S.data[2]))
          @vprint :GaloisGroup 2 "bad evaluation point: 2-sum polynomial has wrong factorisation\n"
          continue
        end
      end
    else
      #should be in the An/Sn case
      #but the parity is already checked, so group is correct
    end

    C.data[3] = overC
    #then verify the descent-chain in S
    C.chn = typeof(S.chn)()
    more_t = false
    for (U, I, ts, cs) in S.chn
      B = upper_bound(C, I, ts)
      prc = bound_to_precision(C, B)
      prc = (ceil(Int, prc[1]*1.2)+1, ceil(Int, prc[2]*1.2)+1) #to have some space for failed evals
      act_prc = (2,4)
#      act_prc = prc
      more_prec = true
      while more_prec
#        @show act_prc, prc
        r = roots(C, act_prc)
        if ts != gen(parent(ts))
          r = map(ts, r)
        end
        for s = cs
          a = evaluate(I^(s*inv(pr)), r)

          fl, val = isinteger(C, B, a)
          if !fl
            if act_prc[1] >= prc[1] && act_prc[2] >= prc[2]
              @vprint :GaloisGroup -2 "bad evaluation point: invariant yields no root\n"
              more_t = true
              break
            end
            act_prc = (min(act_prc[1]*2, prc[1]), min(act_prc[2]*2, prc[2]))
            @vprint :GaloisGroup 2 "isinteger failed, increasing precision to $act_prc\n"
            more_prec = true
            break
          end
          if act_prc[1] < prc[1] || act_prc[2] < prc[2]
            if degree(val) > min(act_prc[2]-3, act_prc[2]*0.9) #TODO: also check on the p-adic side
              act_prc = (min(2*act_prc[1], prc[1]), min(2*act_prc[2], prc[2]))
              @vprint :GaloisGroup 2 "isinteger passed, but result too large, increasing precision to $act_prc\n"
              more_prec = true
              break
            end
            more_prec = false
          else
            more_prec = false
          end
        end
        more_t && break
        more_prec && continue
        push!(C.chn, (U, I, ts, [s*inv(pr) for s = cs]))
      end
      more_t && break
    end
    more_t && continue
    C.G = Gal^inv(pr)
    if overC
      F = GroupFilter() # no filter: need also intransitive groups
                        # could restrict (possibly) to only those
                        # cannot use short_cosets...
      descent(C, C.G, F, one(C.G), grp_id = x->order(x))
      return C.G, C, fixed_field(S, C.G^pr)
    end
    return C.G, C
  end
  #if any fails, "t" was bad and I need to find a way of restarting
end

"""
    subfields(FF:Generic.FunctionField{fmpq})

For a finite extensino of the univariate function field over the rationals, 
find all subfields. The implemented algorithm proceeds by substituting
the transcendental element to an integer, then computing the subfields
of the resulting number field and lifting this information.

It is an adaptaion of Klueners.
"""
function Hecke.subfields(FF::Generic.FunctionField{fmpq})
  C, K, p = _galois_init(FF)
  f = C.C.f

  @vtime :Subfields 1 S = subfields(K)

  Zx = Hecke.Globals.Zx

  prec = (2,1)

  @vprint :Subfields, 1, "Computing subfields of the number field...\n"
  @vprint :Subfields, 2, "obtaining roots with minimal prec $prec\n"
  @vtime :Subfields  2  R = roots(C, prec)
  
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
    fl = issubfield(FF, C, bs, ts = ts)
    if fl === nothing
      continue
    end
    ps, emb = fl
    push!(res, (FunctionField(ps, "a", cached = false)[1], emb))
  end
  return res
end 

function issubfield(FF::Generic.FunctionField, C::GaloisCtx, bs::Vector{Vector{Int}}; ts::fmpz_poly = gen(Hecke.Globals.Zx))    

  SL = SLPolyRing(ZZ, length(bs)*length(bs[1]))
  sx = gens(SL)

  @vprint :Subfields 2 "trying sum as primitive element...\n"
  r = roots(C, (2,1))
  F = parent(r[1]) # should be Qq<<t>>
  if ts != gen(parent(ts))
    r = map(ts, r)
    sx = map(ts, sx)
  end
  
  #TODO: test over finite field first and make sure p is large enough
  con = [sum(r[b]) for b = bs]
  conI = [sum(sx[b]) for b = bs]
  @assert con == [evaluate(x, roots(C, (2,1))) for x= conI]

  if length(Set(con)) < length(bs)
    @vprint :Subfields 2 "...failed, product next...\n"
    con = [prod(r[b]) for b = bs]
    conI = [prod(sx[b]) for b = bs]
    @assert con == [evaluate(x, roots(C, (2,1))) for x= conI]
    while length(Set(con)) < length(bs)
      @vprint :Subfields 2 "...failed, adding 1...\n"
      r .+= 1
      sx .+= 1
      con = [prod(r[b]) for b = bs]
      conI = [prod(sx[b]) for b = bs]
      @assert con == [evaluate(x, roots(C, (2,1))) for x= conI]
    end
  end

  @vprint :Subfields 2 "now proper bounds (for subfield poly)\n"
  B = upper_bound(C, power_sum, conI, length(conI))
  B = 5*B^2
  @vprint :Subfields 2 "coeffs $B\n"
  prec_poly = bound_to_precision(C, B)
  B_poly = B
  @vprint :Subfields 2 "gives a precision of $prec_poly\n"

  function get_poly(prec::Tuple{Int, Int})
    R = roots(C, prec)
    con = [evaluate(c, R) for c = conI]
    if length(Set(con)) < length(bs)
      return nothing
    end


    @vprint :Subfields 2 "building power sums (traces)\n"
    pow = copy(con)
    @assert length(Set(pow)) == length(bs)
    fl, tt = isinteger(C, B_poly, sum(pow))
    if !fl
      return nothing
    end
    tr = [tt]
    while length(tr) < length(bs)
      pow .*= con
      fl, tt = isinteger(C, B_poly, sum(pow))
      if !fl
        return nothing
      end
      push!(tr, tt)
    end

    Qt = base_ring(FF)
    @vprint :Subfields 2 "Newton relations for the polynomial\n"
    local ps
    try
      @vtime :Subfields 2 ps = power_sums_to_polynomial(map(x->x(gen(Qt)), tr))  #should be the subfield polynomial
    catch e
      throw(e)
      return nothing
    end
    if any(x->!isone(denominator(x)), coefficients(ps)) ||
       any(x->any(!isone(denominator(y)) for y = coefficients(numerator(x))), coefficients(ps))
      return nothing
    end
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
  c = coefficients(C.C.f, 1)
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
  B_emb = B

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
    #Kronecker presentation
    # it is given here as a poly over the function field
    # step 1: poly over poly{Q}
    fff = map_coefficients(x->numerator(x)(gen(F)), C.f, parent = parent(ff))
    
    # step 2: move to Kronecker, so that we get poly over poly{Z}
    ff = (ff*derivative(fff)) % fff
    # step 3: lift (and put the denominator back in)
    em = map(x->isinteger(C, B_emb, x), coefficients(ff))
    if any(x->!x[1], em)
      return nothing
    end
    emb = parent(defining_polynomial(FF))([x[2](gen(Nemo.base_ring(FF))) for x in em])(gen(FF)) // derivative(defining_polynomial(FF))(gen(FF))
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
        return nothing
      end
      pr = (2*pr[1], 2*pr[2])
      continue
    end
    emb = get_emb(pr)
    if emb === nothing || !iszero(ps(emb)) 
      if pr[1] >= prec_emb[1] && pr[2] >= prec_emb[2]
        @vprint :Subfields 2 "block system did not give (working) embedding\n"
        return nothing
      end
      pr = (2*pr[1], 2*pr[2])
      continue
    end
    break
  end

  if ps === nothing || emb === nothing
    return nothing
  end

  @vprint :Subfields 1 "have subfield, defined by $ps with embedding $emb\n"
  return ps, emb
end

function isinteger(G::GaloisCtx, B::BoundRingElem{Tuple{fmpz, Int, fmpq}}, r::Generic.RelSeries{qadic})
#  @show "testing", r, "against", B
  p = bound_to_precision(G, B)
  p2 = min(p[2], precision(r))

  Qx = parent(numerator(gen(base_ring(G.f))))
  if iszero(r) 
    return true, Qx(0)
  end
  if r.length + r.val > p2 
    return false, Qx(0)
  end
  f = Qx()
  x = gen(parent(f))
  xpow = parent(x)(1)

  if G.data[3]
    return true, x
  end
  
  for i = 0:(r.length + r.val - 1)
    c = coeff(r, i)
    pr = prime(parent(c))
    
    if c.length < 2 || all(x->iszero(coeff(c, x)), 1:c.length-1)
      cc = coeff(c, 0)
      l = Hecke.mod_sym(lift(cc), pr^precision(cc))
      if abs(l) > pr^p[1]
        return false, x
      end
#      if p[1] > precision(cc)
#        return false, x
#      end
      f += xpow* l
    else
      return false, x
    end
    
    xpow *= x
  end
  return true, f(gen(parent(f))-G.data[2]) #.. and unshift
end

function Hecke.newton_polygon(f::Generic.Poly{Generic.Rat{fmpq}})
  pt = Tuple{Int, Int}[]
  for i=0:degree(f)
    c = coeff(f, i)
    if !iszero(c)
      push!(pt, (i, -degree(numerator(c))+degree(denominator(c))))
    end
  end
  return Hecke.lower_convex_hull(pt)
end

function valuations_of_roots(f::Generic.Poly{Generic.Rat{T}}) where {T}
  return [(slope(l), length(l)) for l = Hecke.lines(Hecke.newton_polygon(f))]
end

Hecke.lines(P::Hecke.Polygon) = P.lines
slope(l::Hecke.Line) = l.slope


