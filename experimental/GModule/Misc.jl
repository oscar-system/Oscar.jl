function Hecke.minpoly(a::QabElem)
  return minpoly(data(a))
end

function Hecke.number_field(::FlintRationalField, a::QabElem; cached::Bool = false)
  f = minpoly(a)
  k, b = number_field(f, check = false, cached = cached)
  return k, b
end

Hecke.minpoly(a::qqbar) = minpoly(Hecke.Globals.Qx, a)

function Hecke.number_field(::FlintRationalField, a::qqbar; cached::Bool = false)
  f = minpoly(a)
  k, b = number_field(f, check = false, cached = cached)
  Qx = parent(k.pol)
  function to_k(x::qqbar)
    if x == a
      return b
    end
    f = minpoly(x)
    r = roots(f, k)
    pr = 10
    while true
      C = AcbField(pr)
      ca = C(a)
      lp = findall(i->contains_zero(Qx(i)(ca) - C(x)), r)
      if length(lp) == 1
        return r[lp[1]]
      end
      if length(lp) == 0
        error("not in the image")
      end
      pr *= 2
      @assert pr < 2^16
    end
  end
  function to_qqbar(x::nf_elem)
    return Qx(x)(a)
  end
  #TODO: make map canonical?
  # ... and return gen(k) instead?
  return k, MapFromFunc(to_qqbar, to_k, k, parent(a))
end

Base.getindex(::FlintRationalField, a::qqbar) = number_field(QQ, a)

function Hecke.numerator(f::fmpq_poly, parent::FmpzPolyRing = Hecke.Globals.Zx)
  g = parent()
  ccall((:fmpq_poly_get_numerator, Nemo.libflint), Cvoid, (Ref{fmpz_poly}, Ref{fmpq_poly}), g, f)
  return g
end

function cyclo_fixed_group_gens(a::nf_elem)
  C = parent(a)
  fl, f = Hecke.is_cyclotomic_type(C)
  @assert fl
  if isone(f)
    return [(1,1)]
  end
  p = first(PrimesSet(1000, -1, f, 1))
  k = GF(p)
  o = k(1)
  lf = factor(f)
  test = [div(f, x) for x = keys(lf.fac)]
  while any(i->isone(o^i), test)
    o = rand(k)^divexact(p-1, f)
  end
  poly_a = numerator(parent(C.pol)(a))
  #the conjugates will be o^j for all j coprime to f
  #the plan is to check conjugates being equal mod p and then
  #verify exactly. Alternatively one could make sure that p is
  #large enough...
  bs = Dict{typeof(o), Set{Int}}()
  co = Set([j for j=1:f-1 if gcd(j, f) == 1])
  all_aut = Set([1])
  gn = Tuple{Int, Int}[]
  while length(co) > 0
    j = pop!(co)
    c = poly_a(o^j)
    if haskey(bs, c)
      d = poly_a(gen(C)^j)
      @assert d == poly_a(gen(C)^first(bs[c]))
      #now we have that conj[i] == conj[j], so 
      #I'd like all autmorphisms mapping zeta^i -> zeta^j
      #given aut: zeta -> zeta^k, so
      #zeta^i -> zeta^(ki) and ki = j mod f
      #so k = j * modinv(i, f)
      #but if aut fixes, then all powers will also fix.
      #that should deal with many conjugates...
      m = Set{Int}()
      for i = bs[c]
        k = (invmod(i, f)*j) % f
        if k in all_aut
#          @show "already have $k"
          continue
        end
        push!(gn, (k, f))
        while true
          n = [(k*x) % f for x = all_aut]
          sz = length(all_aut)
          push!(all_aut, n...)
          push!(m, n...)
          if length(all_aut) == sz
            break
          end
        end
      end
      for i = copy(bs[c])
        for k = m
          ki = (k*i) % f
          if ki in co
            pop!(co, ki)
          end
          push!(bs[c], ki)
        end
      end
    else
      bs[c] = Set([j])
    end
  end
  return gn
end

function Hecke.preimage(h::GrpAbFinGenMap, u::GrpAbFinGen, L::Hecke.GrpAbLattice = Hecke.GroupLattice)

  fl, f = is_subgroup(u, codomain(h))
  @assert fl
  k, mk = kernel(h)
  return sub(domain(h), vcat(map(mk, gens(k)), [preimage(h, x) for x = map(f, gens(u))]))
end

function cyclo_fixed_group_gens(A::AbstractArray{nf_elem})
  if length(A) == 0
    return [(1,1)]
  end
  G = map(cyclo_fixed_group_gens, A)
  F = lcm([x[1][2] for x = G])
  R, mR = unit_group(quo(ZZ, F)[1])
  Q, mQ = sub(R, gens(R))
  for s = G
    zf = quo(ZZ, s[1][2])[1]
    S, mS = unit_group(zf)
    u, mu = sub(S, [preimage(mS, zf(x[1])) for x = s])
    h = hom(R, S, [preimage(mS, mR(R[i])) for i=1:ngens(R)])
    Q = intersect(Q, preimage(h, u)[1])
  end
  s, ms = snf(Q)
  _, sR = is_subgroup(Q, R)
  Qgen = map(sR, map(ms, gens(s)))
  lf = factor(F)
  for p = keys(lf.fac)
    while true
      G = div(F, p)
      zg = quo(ZZ, G)[1]
      S, mS = unit_group(zg)
      hRS = hom(R, S, [preimage(mS, zg(mR(R[i]))) for i = 1:ngens(R)])
      if length(quo(R, Qgen)[1]) == length(quo(S, map(hRS, Qgen))[1])
        R, mR = S, mS
        Qgen = map(hRS, Qgen)
        F = G
      else
        break
      end
    end
  end
  return [(mR(sR(ms(x))), F) for x = gens(s)]
end

function Hecke.number_field(::FlintRationalField, a::AbstractVector{QabElem}; cached::Bool = false)
  if length(a) == 0
    return Hecke.rationals_as_number_field()[1]
  end
  f = lcm([Hecke.is_cyclotomic_type(parent(data(x)))[2] for x = a])
  K = cyclotomic_field(f)[1]
  k, mkK = Hecke.subfield(K, [K(data(x)) for x = a])
  return k, gen(k)
end

Base.getindex(::FlintRationalField, a::QabElem) = number_field(QQ, a)
Base.getindex(::FlintRationalField, a::Vector{QabElem}) = number_field(QQ, a)
Base.getindex(::FlintRationalField, a::QabElem...) = number_field(QQ, [x for x =a])
Base.getindex(::FlintRationalField, chi::Oscar.GAPGroupClassFunction) = number_field(QQ, chi)

