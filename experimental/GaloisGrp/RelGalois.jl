function galois_group(mkK::Map{AnticNumberField, AnticNumberField})
  k = domain(mkK)
  K = codomain(mkK)
  
  G, C = galois_group(K)
  p = C.C.p
  rt = roots(C, 5)
  f = parent(defining_polynomial(K))(mkK(gen(k)))
  s = map(f, rt)
  @assert length(Set(s)) == degree(k)
  @show b = findall(x->x == s[1], s)
  @show H = stabilizer(G, b, on_sets)[1]
  h = action_homomorphism(H, b)
  return image(h)[1]
end

function galois_group(K::Hecke.SimpleNumField{nf_elem})
  f = defining_polynomial(K)
  @assert degree(f) >= 2
  k = base_ring(f)
  bp = (1,1)
  ct = CycleType[]

  local zk::NfOrd
  local den::nf_elem
  if Hecke.ismaximal_order_known(k)
    zk = maximal_order(k)
    if isdefined(zk, :lllO)
      zk = zk.lllO::NfOrd
    end
    den = k(1)
  else
    zk = any_order(k)
    den = derivative(defining_polynomial(k))(gen(k))
  end
  zk = lll(zk) # always a good option!
  p = degree(f)
  f *= lcm(map(denominator, coefficients(f)))
  np = 0
  while true
    @vprint :PolyFactor 3 "Trying with $p\n "
    p = next_prime(p)
    if !Hecke.isprime_nice(zk, p)
      continue
    end
    P = prime_decomposition(zk, p, 1) #not quite sure... but lets stick with it
    if length(P) == 0
      continue
    end
    F, mF1 = Hecke.ResidueFieldSmallDegree1(zk::NfOrd, P[1][1])
    mF = extend(mF1, k)
    fp = map_coefficients(mF, f, cached = false)
    if degree(fp) < degree(f) || iszero(constant_coefficient(fp)) 
      continue
    end
    if !issquarefree(fp)
      continue
    end
    lf = factor_shape(fp)
    push!(ct, CycleType(vec(collect(keys(lf)))))
    dg = order(ct[end])
    if order(ct[end]) == 1
      continue
    end
    if bp == (1,1) 
      bp = (p, dg, P[1][1])
    elseif bp[2] > dg 
      bp = (p, dg, P[1][1])
    end
    if ceil(Int, degree(f)/4) <= bp[2] <= floor(Int, degree(f)/2) || length(ct) > 2*degree(f)
      break
    end
  end

  #next: need the prime (and the completion of k at that prime)
  P = bp[3]
  C, mC = completion(k, P)
  setprecision!(C, 10)
  D = QadicField(Int(bp[1]), Int(bp[2]), 10)[1]
  mCD = MapFromFunc(x->D(coeff(x, 0)), C, D)
  H = Hecke.HenselCtxQadic(map_coefficients(x->mCD(mC(x)), f))
  V = Hecke.vanHoeijCtx()
  V.P = P
  V.H = H
  C = GaloisCtx(typeof(V))
  C.f = f
  C.C = V
  C.data = [mC, 1, mCD, Hecke.norm_change_const(zk), den]
  #data:
  # [1] map into the completion. Careful, currently Q_p is returned and used as a deg-1 extension (ie. Qq)
  #     use only coeff(, 0)
  # [2] the current precision of the roots
  # [3] the embedding of completion of k into completion of K
  # [4] as the namse suggest
  # [5] the denominator used, f'(alpha), the Kronnecker rep

  C.B = add_ring()(maximum([iroot(ceil(fmpz, length(x)), 2)+1 for x = coefficients(f)]))
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
  println(io, "GaloisCtx for computations modulo $(C.C.P)")
end

function Oscar.roots(C::GaloisCtx{Hecke.vanHoeijCtx}, pr::Int = 5; raw::Bool = true)
  #TODO: deal with raw and denominators.
  if pr > C.data[2]
    setprecision!(codomain(C.data[1]), pr)
    setprecision!(codomain(C.data[3]), pr)
    C.C.H.f = map_coefficients(x->C.data[3](C.data[1](x)), C.f)
    Hecke.grow_prec!(C.C, pr)
    C.data[2] = pr
  end
  return roots(C.C)
end
  

function isinteger(C::GaloisCtx{Hecke.vanHoeijCtx}, y::BoundRingElem{fmpz}, x::qadic)
  P = C.C.P
  zk = order(P)
  if any(i->!iszero(coeff(x, i)), 1:length(x)-1)
    @show :wrongDegree
    return false, zero(nf(zk))
  end
  mkc = C.data[1]
  c = codomain(mkc)
  d = codomain(C.data[3])
  @assert parent(x) == d
  P = C.C.P
  zk = order(P)
  x *= map_coeff(C, C.data[5]) #the den
  @show a = nf(zk)(Hecke.reco(zk(preimage(mkc, c(coeff(x, 0)))), C.C.Ml, C.C.pMr))
  @show a = a*inv(C.data[5])
  if ceil(fmpz, length(a)) <= value(y)^2
    return true, a
  else
    @show :tooLarge
    return false, a
  end
end

#TODO: maybe add precision argument, see comment in main file.
function map_coeff(C::GaloisCtx{Hecke.vanHoeijCtx}, x)
  return coeff(C.data[1](x), 0)
end

function bound_to_precision(C::GaloisCtx{Hecke.vanHoeijCtx}, y::BoundRingElem{fmpz}, extra::Int = 0)
  #the bound is a bound on the sqrt(T_2(x)). This needs to be used with the norm_change stuff
  #and possible denominators and such. Possibly using Kronnecker...
  c1, c2 = C.data[4] # the norm-change-const
  @show v = value(y) + iroot(ceil(fmpz, length(C.data[5])), 2)+1 #correct for den
  #want to be able to detect x in Z_k of T_2(x) <= v^2
  #if zk = order(C.C.P) is (known to be) maximal, 2-norm of coeff. vector squared < c2*v^2
  #otherwise, we need tpo multiply by f'(alpha) (increasing the size), revover and divide
  #needs to be corrdinated with isinteger

  #step 2: copy the precision estimate from NumField/NfAbs/PolyFactor.jl
  #        (and check if all logs should really be log2?)
  P = C.C.P
  zk = order(P)
  k = nf(zk)
  @show N = ceil(Int, degree(k)/2/log(norm(P))*(log2(c1*c2) + 2*nbits(v)))
  return N
end
