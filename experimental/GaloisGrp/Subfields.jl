
function embedding_hom(k, K)
  return MapFromFunc(x->K(x), k, K)
end

const BlockSystem_t = Vector{Vector{Int}}

#TODO: create complete lattice and record it being complete?
#      intersect SubfieldLattice with subfield to get a lattice for a smaller
#         field? For recursion in Galois?
#      make sure this works also for Gal(K[t]) and Q(t)
#      possible function bound(::GaloisCtx, ::CoeffRingElem)?
#                                           ::FieldElem)?
#      hove down to Hecke?
mutable struct SubfieldLattice{T}
  K::T
  P::POSet{BlockSystem_t}
  l::Dict{POSetElem{BlockSystem_t}, Any} # Any: Tuple{T, Map{T, T}}
  G::GaloisCtx #for the roots and such

  function SubfieldLattice(K::T, G::GaloisCtx) where {T}
    r = new{T}()
    r.K = K
    r.G = G
    n = degree(K)
    bs = [[i] for i=1:n]
    cs = [[i for i=1:n]]
    r.P = POSet([bs, cs],
                  (x,y) -> issubset(x[1], y[1]) || issubset(y[1], x[1]),
                  (x,y) -> issubset(x[1], y[1]) - issubset(y[1], x[1]))

    r.l = Dict(POSetElem(r.P, 1) =>(K, id_hom(K)),
               POSetElem(r.P, 2) => (base_ring(K), 
                           embedding_hom(base_ring(K), K)))
    return r
  end
end

function Base.show(io::IO, S::SubfieldLattice)
  print(io, "Subfield lattice for $(S.K)")
end

function field(S::SubfieldLattice)
  return S.K
end

function POSet(S::SubfieldLattice)
  return S.P
end

function GaloisCtx(S::SubfieldLattice)
  return S.G
end

function block_system(G::GaloisCtx, a::SimpleNumFieldElem)
  K = parent(a)
  kx = parent(defining_polynomial(K))
  f = kx(a)
  pr = 1
  while true
    r = roots(G, pr, raw = true)
    c = map(f, r) # TODO: use the embedding map!
    bs = Hecke.MPolyFact.block_system(c)
    if all(x->length(x) == length(bs[1]), bs) 
      return bs
    end
    pr *= 2
#    error("adada")
    @show bs, pr
    if pr > 100
      error("too bad")
    end
  end
end

mutable struct SubfieldLatticeElem{T} 
  p::SubfieldLattice{T}
  b::BlockSystem_t
end

Base.getindex(S::SubfieldLattice, i::Int) = SubfieldLatticeElem(S, Oscar.data(S.P[i]))
Base.length(S::SubfieldLattice) = length(S.P)

function Oscar.parent(A::SubfieldLatticeElem)
  return A.p
end

import Base: ==
function ==(A::SubfieldLatticeElem, B::SubfieldLatticeElem)
  S = parent(A)
  @assert parent(A).K == parent(B).K
  P = S.P
  bs = A.b
  cs = B.b
  return P.can_cmp(bs, cs) && P.cmp(bs, cs) == 0
end

function Base.:*(A::SubfieldLatticeElem, B::SubfieldLatticeElem)
  S = parent(A)
  @assert parent(A).K == parent(B).K
  bs = A.b
  cs = B.b
  ds = [intersect(b, c) for b = bs for c = cs]
  ds = [x for x = ds if length(x) > 0]
  return S(ds)
end

function Base.intersect(A::SubfieldLatticeElem, B::SubfieldLatticeElem)
  S = parent(A)
  @assert parent(A).K == parent(B).K
  bs = A.b
  cs = B.b
  ds = deepcopy(bs)
  while true
    n = length(ds[1])
    for d=ds
      i = findall(x->any(y->y in d, x), cs)
      x = Set(d)
      union!(x, cs[i]...)
      empty!(d)
      append!(d, x)
    end
    ds = vec(collect(Set(ds)))
    if length(ds[1]) == n
      return S(ds)
    end
    n = length(ds[1])
    for d=ds
      i = findall(x->any(y->y in d, x), bs)
      x = Set(d)
      union!(x, cs[i]...)
      empty!(d)
      append!(d, x)
    end
    ds = vec(collect(Set(ds)))
    if length(ds[1]) == n
      return S(ds)
    end
  end
end

function (S::SubfieldLattice{T})(mK::Map) where T
  k = domain(mK)
  g = mK(gen(k))
  bs = block_system(S.G, g)

  if !(bs in S.P)
    push!(S.P, bs)
    S.l[S.P(bs)] = (k, mK)
  end
  return SubfieldLatticeElem(S, bs)
end

function (S::SubfieldLattice{T})(bs::BlockSystem_t) where T
  if !(bs in S.P)
    push!(S.P, bs)
  end
  return SubfieldLatticeElem(S, bs)
end

function Base.push!(S::SubfieldLattice, mK::Map)
  S(mK)
  nothing
end

function Base.push!(S::SubfieldLattice, v::Vector)
  vv = map(x->S(block_system(GaloisCtx(S), x)), v)
  return prod(vv)
end

"""
For a (potential) block system `bs` either find the corresponding subfield,
thus proving the block system to be valid, or return `nothing` showing the
block system to be wrong.
"""
function subfield(S::SubfieldLattice, bs::BlockSystem_t)
  if bs in S.P && haskey(S.l, S.P(bs))
    return S.l[S.P(bs)]
  end
  pr = 5
  local func

  R = roots(GaloisCtx(S), pr, raw = true)
  s = Set{Int}()
  union!(s, bs...)

  if s != Set(1:length(R))
    return nothing
  end

  if any(x->length(x) != length(bs[1]), bs)
    return nothing
  end

  X = gens(SLPolyRing(ZZ, length(R)))

  r = [sum(R[b]) for b = bs]
  if length(Set(r)) == length(bs)
    func = [sum(X[b]) for b = bs]
  else
    k = 0
    while true
      r = [prod(R[b] .+ k) for b = bs]
      if length(Set(r)) == length(bs)
        let k = k
          func = [prod(X[b] .+ k) for b = bs]
          break
        end
      end
      k += 1
    end
  end

  G = GaloisCtx(S)
  #power sums: degree of k is length(bs), so need
  #Tr(beta^i) for i=1:length(bs)
  B = upper_bound(G, power_sum, func, length(bs))
  pr = bound_to_precision(G, B)
  R = roots(G, pr, raw = true)
  beta = [evaluate(f, R) for f = func]
  pow = copy(beta)

  K = field(S)
  k = base_ring(K)
  if !isa(k, AbstractAlgebra.Field)
    k = QQ
  end
  tr = [k(isinteger(G, B, sum(beta))[2])]
  while length(tr) < length(bs)
    pow .*= beta
    fl, v = isinteger(G, B, sum(pow))
    fl || return nothing
    push!(tr, k(v))
  end

  g = power_sums_to_polynomial(tr)

  #from Qt:
  # beta = 1/f'(alpha) sum b_i alpha^i
  # f(alpha)/(t-alpha) = sum g_i(alpha) t^i
  # Tr(beta * g_i(alpha)/f'(alpha)) = b_i (dual basis)
  Kt, t = PolynomialRing(K, "t", cached = false)
  Gk = divexact(map_coefficients(K, defining_polynomial(K), parent = Kt), t-gen(K))
  Qt = parent(defining_polynomial(K))
  Gt = [Qt(x) for x = coefficients(Gk)]
  fsa = derivative(defining_polynomial(K))(gen(K))
  fsat = Qt(fsa)
  B = length(bs)*evaluate(func[1], [G.B for x = R])*parent(B)(maximum(ceil(fmpz, length(x)) for x = coefficients(Gk)))
  pr = bound_to_precision(G, B)
  R = roots(G, pr, raw = true)
  beta = K()
  for k=0:degree(K)-1
    fl, v = isinteger(G, B, sum(evaluate(func[i], R)*sum(Gt[k+1](R[bs[i][j]]) for j=1:length(bs[1])) for i=1:length(bs)))
    fl || return nothing
    beta += gen(K)^k//fsa*v
  end
  if iszero(g(beta))
    k, a = number_field(g)
    h = hom(k, K, beta)
    if !(bs in S.P)
      push!(S.P, bs)
    end
    S.l[S.P(bs)] = (k, h)
    return k, h
  end
  return nothing
end

function subfield(S::SubfieldLatticeElem)
  return subfield(parent(S), S.b)
end

function Oscar.degree(S::SubfieldLatticeElem)
  return length(S.b)
end

###################################
# van Hoeij, Novacin, Klueners
###################################

function _subfields(K::AnticNumberField; pStart = 2*degree(K)+1, prime = 0)
  Zx = Hecke.Globals.Zx

  f = Zx(mapreduce(denominator, lcm, coefficients(defining_polynomial(K)), init = fmpz(1))*defining_polynomial(K))
  f = divexact(f, content(f))

  p, ct = find_prime(Hecke.Globals.Qx(f), pStart = pStart, prime = prime,
                                          filter_pattern = x->any(t->degree(t) == 1, keys(x.fac)))
  n = degree(K)
  if primitive_by_shape(ct, n)
    return nothing
  end
  G = GaloisCtx(f, p)
  S = SubfieldLattice(K, G)

  pr = 5
  nf = sum(x*x for x = coefficients(f))
  B = degree(f)^2*(iroot(nf, 2)+1) #from Paper: bound on the coeffs we need
  @show B = B^2 # Nemo works with norm-squared....
  @show "using ", p
  @show pr = clog(B, p) 
  @show pr *= div(n,2)
  @show pr += 2*clog(2*n, p)
  H = factor_mod_pk_init(f, p)
  @show factor_mod_pk(H, 1)

  b = basis(K)
  b .*= inv(derivative(f)(gen(K)))
  bt = [parent(defining_polynomial(K))(x) for x = b]
  bd = map(denominator, bt)
  bz = [Zx(bd[i]*bt[i]) for i=1:length(bd)]
  @assert parent(f) == parent(bz[1])

  lf = factor_mod_pk(Array, H, 1)
  #the roots in G are of course roots of the lf[i]
  #roots that belong to the same factor would give rise
  #to the same principal subfield. So we can save on LLL calls.
  r = roots(G, 1)
  F, mF = ResidueField(parent(r[1]))
  r = map(mF, r)
  rt_to_lf = [findall(x->iszero(f[1](x)), r) for f = lf]
  done = zeros(Int, length(lf))
  while true
    lf = factor_mod_pk(Array, H, pr)
    ppr = fmpz(p)^pr
    @assert parent(lf[1][1]) == parent(f)
    @assert all(x->is_monic(x[1]), lf)
    rt = findfirst(x->degree(x[1]) == 1, lf)
    done[rt] = 1
    for i=1:length(lf)
      i == rt && continue
      done[i] == 1 && continue
      di = degree(lf[i][1])
      #lf[rt][1] is linear, we need bd evaluated at the root
      #lf[i] [1] is gen, we need bd reduced mod this poly
      # then the difference, put into a matrix...
      R = -constant_coefficient(lf[rt][1])
      den = [invmod(x, ppr) for x = bd]
      M = zero_matrix(ZZ, n, di)
      for j=1:n
        id = bz[j](R)*den[j]
        phi = rem(bz[j], lf[i][1])*den[j] - id
        mod_sym!(phi, ppr)
        for k=1:di
          M[j, k] = coeff(phi, k-1)
        end
      end
      M = [M identity_matrix(ZZ, n); ppr*identity_matrix(ZZ, di) zero_matrix(ZZ, di, n)]
      D = diagonal_matrix(vcat([B for j=1:di], [fmpz(1) for j=1:n]))
#      M = M*D
      while true
        @show maximum(nbits, M), nbits(B), size(M)
        global last_M = M
        #TODO: possible scale (and round) by 1/sqrt(B) so that
        #      the lattice entries are smaller (ie like in the 
        #      van Hoeij factoring)
        @time r, M = lll_with_removal(M, B, lll_ctx(0.501, 0.75))
        M = M[1:r, :]
        @show r, i, pr

        if iszero(M[:, 1:di])
          if n % r == 0
            gens = [sum(M[k, di+j] * b[j] for j=1:n) for k=1:r]
            #need the subfield poly - or a guess...
            #this is an upper bound.
            Qt = parent(defining_polynomial(K))
            gz = [numerator(Qt(x)) for x = gens] 
            @assert all(x->!iszero(denominator(x) % p), gens)
            T = Int[]
            ggz = [x(R) for x = gz]
            for j=1:length(lf)
              hz = [x % lf[j][1]  for x = gz] .- ggz
              for k=1:length(hz)
                mod_sym!(hz[k], ppr)
              end
              if all(iszero, hz)
                push!(T, j)
              end
            end
            @assert i in T
            @assert rt in T
            @show T
            #the paper says the (unknown) (algebraic) coefficients of
            # prod(f[i] i in T)
            #have to generate the subfield and f in K[t] is is_irreducible
            #so deg(poly)*deg(subfield) = deg(K)
            #deg subfield should be r, so
            #if precision is too low, then this poly will be too
            #large (T too large), hence an exclusion
            #I think thay if the degree is correct, the subfield is correct
            #thus one should be able to get the 1st block as well
            #However, instead of blocks, T is as good an indicator for 
            #subfields
            if sum(degree(lf[t][1]) for t = T) *r != degree(K)
              @show :subfield_pol_wrong
            else
              E = push!(S, gens)
              if S.P(E.b) in keys(S.l)
                @show "field already known"
              end
              if degree(E) != length(gens)
                @show :wrong_block
              elseif subfield(E) === nothing
                @show :no_subfield
              else
                for j in T
                  for h in rt_to_lf[j]
                    done[j] = 1
                  end
                end
                done[i] = 1
              end
            end
          else
            @show :flop
          end
          break
        else
          @show :scale
          M = M*D
        end
      end
    end
    pr *= 2
    pr < 2^15 || error("ada")
    all(isequal(1), done) && return S
  end
end
