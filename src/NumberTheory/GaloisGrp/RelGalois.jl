function galois_group(mkK::Map{AbsSimpleNumField, AbsSimpleNumField})
  k = domain(mkK)
  K = codomain(mkK)
  
  G, C = galois_group(K)
  p = C.C.p
  rt = roots(C, 5)
  f = parent(defining_polynomial(K))(mkK(gen(k)))
  s = map(f, rt)
  @assert length(Set(s)) == degree(k)
  b = findall(==(s[1]), s)
  H = stabilizer(G, b, on_sets)[1]
  h = action_homomorphism(H, b)
  return image(h)[1]
end

mutable struct recoCtx
  Ml::ZZMatrix
  pMr::Tuple{ZZMatrix, ZZRingElem, Hecke.fmpz_preinvn_struct}
  pM::Tuple{ZZMatrix, ZZRingElem}
end

function GaloisCtx(f::PolyRingElem{AbsSimpleNumFieldElem}, P::AbsSimpleNumFieldOrderIdeal)
  k = base_ring(f)
  @assert k == Hecke.nf(order(P))
  zk = order(P)
  C, mC = Hecke.completion_easy(k, P)
  setprecision!(C, 10)
  F, mF = Hecke.ResidueFieldSmallDegree1(order(P), P)
  p = Int(characteristic(F))
  mF = Hecke.extend(mF, k)
  d = reduce(lcm, keys(factor_shape(map_coefficients(mF, f))))
  D = QadicField(p, d, 10)[1]
  mCD = MapFromFunc(C, D, x->D(coeff(x, 0)))
  H = Hecke.HenselCtxQadic(map_coefficients(x->mCD(mC(x)), f))
  V = Hecke.vanHoeijCtx()
  V.P = P
  V.H = H
  C = GaloisCtx(typeof(V))
  C.prime = P
  C.f = f
  C.C = V
  if Hecke.is_maximal_order_known(k)
    den = k(1)
  else
    den = derivative(defining_polynomial(k))(gen(k))
  end
  C.data = [mC, 1, mCD, Hecke.norm_change_const(zk), den, Dict{Int, recoCtx}()]
  #data:
  # [1] map into the completion. Careful, currently Q_p is returned and used as a deg-1 extension (ie. Qq)
  #     use only coeff(, 0)
  # [2] the current precision of the roots
  # [3] the embedding of completion of k into completion of K
  # [4] as the name suggest
  # [5] the denominator used, f'(alpha), the Kronnecker rep
  #     careful: there is (might be) a den from the non-monic poly. This
  #     is "dealt with" in roots: the roots are automatically scaled by this
  #     so they are mathematically integral
  #     the den here is from any order to max order used to recognize 
  #     integral elements in isinteger
  # [6] Dict: mapping precision to reco data

  C.B = add_ring()(maximum([iroot(ceil(ZZRingElem, length(x)), 2)+1 for x = coefficients(f)]))
  return C
end

function find_prime(f::PolyRingElem{AbsSimpleNumFieldElem}, extra::Int = 5; pStart::Int = degree(f)+1, prime::Any = 0)
  if prime !== 0
    @assert isa(prime, AbsSimpleNumFieldOrderIdeal)
    return prime, Set{CycleType}()
  end
  @assert degree(f) >= 2
  k = base_ring(f)
  bp = (1,1)
  ct = Set{CycleType}()

  local zk::AbsSimpleNumFieldOrder
  local den::AbsSimpleNumFieldElem
  if Hecke.is_maximal_order_known(k)
    zk = maximal_order(k)
    if isdefined(zk, :lllO)
      zk = zk.lllO::AbsSimpleNumFieldOrder
    end
    den = k(1)
  else
    zk = any_order(k)
    den = derivative(defining_polynomial(k))(gen(k))
  end
  zk = lll(zk) # always a good option!
  p = pStart
  f *= lcm(map(denominator, coefficients(f)))
  np = 0
  _num = 0
  while true
    @vprint :PolyFactor 3 "Trying with $p\n "
    _num += 1
    p = next_prime(p)
    if !Hecke.is_prime_nice(zk, p)
      continue
    end
    P = prime_decomposition(zk, p, 1) #not quite sure... but lets stick with it
    if length(P) == 0
      continue
    end
    F, mF1 = Hecke.ResidueFieldSmallDegree1(zk::AbsSimpleNumFieldOrder, P[1][1])
    mF = extend(mF1, k)
    fp = map_coefficients(mF, f, cached = false)
    if degree(fp) < degree(f) || iszero(constant_coefficient(fp)) 
      continue
    end
    if !is_squarefree(fp)
      continue
    end
    lf = factor_shape(fp)
    c = CycleType(vec(collect(lf)))
    push!(ct, c)
    dg = order(c)
    if order(c) == 1
      continue
    end
    if bp == (1,1) 
      bp = (p, dg, P[1][1])
    elseif bp[2] > dg 
      bp = (p, dg, P[1][1])
    end
    if ceil(Int, degree(f)/4) <= bp[2] <= floor(Int, degree(f)/2) && length(ct) > 2*degree(f) || _num > 200
      #counter example: C_5: cycle_types are [1,1,1,1,1] or [5], degrees
      #                      1 and 5...
      break
    end
  end

  #next: need the prime (and the completion of k at that prime)
  P = bp[3]
  return P, ct
end

function galois_group(K::Hecke.SimpleNumField{AbsSimpleNumFieldElem}; prime::Any = 0, pStart::Int = degree(K)+1)
  f = defining_polynomial(K)

  P, ct = find_prime(f, prime = prime, pStart = pStart)
  C = GaloisCtx(f, P)

  if an_sn_by_shape(ct, degree(K))
    @vprint :GaloisGroup 1 "An/Sn by cycle type\n"
    if is_square(discriminant(K))
      G = alternating_group(degree(K))
    else
      G = symmetric_group(degree(K))
    end
    C.G = G
    return G, C
  end

  G, F, si = starting_group(C, K)
  return descent(C, G, F, si)
end

function Oscar.roots(H::Hecke.HenselCtxQadic)
  f = factor(H)
  return [-constant_coefficient(x)//leading_coefficient(x) for x = f if degree(x) == 1]
end

function Oscar.roots(V::Hecke.vanHoeijCtx)
  return roots(V.H)
end

function Base.show(io::IO, C::GaloisCtx{Hecke.vanHoeijCtx})
  print(io, "GaloisCtx for computations modulo $(C.C.P)")
end

function Oscar.roots(C::GaloisCtx{Hecke.vanHoeijCtx}, pr::Int = 5; raw::Bool = false)
  #TODO: deal with raw and denominators.
  if pr > C.data[2]
    setprecision!(codomain(C.data[1]), pr)
    setprecision!(codomain(C.data[3]), pr)
    C.C.H.f = map_coefficients(x->C.data[3](C.data[1](x)), C.f)
    Hecke.grow_prec!(C.C, pr)
    C.data[6][pr] = recoCtx(C.C.Ml, C.C.pMr, C.C.pM)
    C.data[2] = pr
  end
  if pr < 0.8*C.data[2]
    if false && !isone(C.data[5]) && !raw
      d = setprecision(coeff(C.data[1](C.data[5]), 0), pr)
      return [d*setprecision(x, pr) for x = roots(C.C)]
    else
      return [setprecision(x, pr) for x = roots(C.C)]
    end
  end
  if false && !isone(C.data[5]) && !raw
    d = coeff(C.data[1](C.data[5]), 0)
    return [d*x for x = roots(C.C)]
  else
    r = roots(C.C)
    return roots(C.C)
  end
end
  

function isinteger(C::GaloisCtx{Hecke.vanHoeijCtx}, y::BoundRingElem{ZZRingElem}, x::QadicFieldElem)
  P = C.C.P
  zk = order(P)
  if any(i->!iszero(coeff(x, i)), 1:length(x)-1)
    return false, zero(Hecke.nf(zk))
  end
  mkc = C.data[1]
  c = codomain(mkc)
  d = codomain(C.data[3])
  @assert parent(x) == d
  P = C.C.P
  zk = order(P)

  x *= map_coeff(C, C.data[5]) #the den of any in max order
  pr = precision(x)
  if !haskey(C.data[6], pr)
    #copied from grow_prec! in Hecke
    # I think the preinvn is not used
    @vtime :PolyFactor 2 X1 = C.C.P^pr
    @vtime :PolyFactor 2 X2 = basis_matrix(X1)
    @vtime :PolyFactor 2 Ml = lll(X2)
    @vtime :PolyFactor 2 F = Hecke.FakeFmpqMat(pseudo_inv(Ml))
    #(M*B)^-1 = B^-1 * M^-1, so I need basis_mat_inv(zk) * pM
    pMr = (F.num, F.den, Hecke.fmpz_preinvn_struct(2*F.den))
    C.data[6][pr] = recoCtx(Ml, pMr, (F.num, F.den))
  end

  r = C.data[6][pr]
  a = Hecke.nf(zk)(Hecke.reco(zk(preimage(mkc, c(coeff(x, 0)))), r.Ml, r.pMr))
  a = a*inv(C.data[5])
  if ceil(ZZRingElem, length(a)) <= value(y)^2
    @hassert :GaloisGroup 2 is_integral(a)
    return true, a
  else
    return false, a
  end
end

#TODO: maybe add precision argument, see comment in main file.
function map_coeff(C::GaloisCtx{Hecke.vanHoeijCtx}, x)
  return coeff(C.data[1](x), 0)
end

function bound_to_precision(C::GaloisCtx{Hecke.vanHoeijCtx}, y::BoundRingElem{ZZRingElem}, extra::Int = 10)
  #the bound is a bound on the sqrt(T_2(x)). This needs to be used with the norm_change stuff
  #and possible denominators and such. Possibly using Kronecker...
  c1, c2 = C.data[4] # the norm-change-const
  v = value(y) + iroot(ceil(ZZRingElem, length(C.data[5])), 2)+1 #correct for den
  #want to be able to detect x in Z_k of T_2(x) <= v^2
  #if zk = order(C.C.P) is (known to be) maximal, 2-norm of coeff. vector squared < c2*v^2
  #otherwise, we need tpo multiply by f'(alpha) (increasing the size), revover and divide
  #needs to be coordinated with isinteger

  #step 2: copy the precision estimate from NumField/NfAbs/PolyFactor.jl
  #step 3: be careful: x*den being integral (and small enough) does not mean x
  #        is correct. Using extra makes this virtually certain though
  #        otherwise, the minpoly (or charpoly) is called for
  P = C.C.P
  zk = order(P)
  k = Hecke.nf(zk)
  N = ceil(Int, degree(k)/2/log(norm(P))*(log(c1*c2) + 2*log(v)))
  return N + extra
end

function galois_group(f::PolyRingElem{AbsSimpleNumFieldElem}, ::QQField)
  @assert is_irreducible(f)

  g = f
  k = 0
  
  p = Hecke.p_start
  F = GF(p)

  Kx = parent(f)
  K = base_ring(Kx)

  Zx = Hecke.Globals.Zx
  local N
  @vtime :PolyFactor Np = Hecke.norm_mod(g, p, Zx)
  while is_constant(Np) || !is_squarefree(map_coefficients(F, Np))
    k = k + 1
    g = compose(f, gen(Kx) - k*gen(K))
    @vtime :PolyFactor 2 Np = Hecke.norm_mod(g, p, Zx)
  end

  @vprint :PolyFactor 2 "need to shift by $k, now the norm\n"
  if any(x -> denominator(x) > 1, coefficients(g)) ||
     !Hecke.is_defining_polynomial_nice(K)
    @vtime :PolyFactor 2 N = Hecke.Globals.Qx(norm(g))
  else
    @vtime :PolyFactor 2 N = Hecke.norm_mod(g, Zx)
    @hassert :PolyFactor 1 N == Zx(norm(g))
  end

  while is_constant(N) || !is_squarefree(N)
    error("should not happen")
  end

  return galois_group(N)
end
