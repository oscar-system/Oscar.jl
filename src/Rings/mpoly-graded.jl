export weight, decorate, ishomogenous, homogenous_components, filtrate,
grade, homogenous_component, jacobi_matrix, jacobi_ideal,
HilbertData, hilbert_series, hilbert_series_reduced, hilbert_series_expanded, hilbert_function, hilbert_polynomial,
homogenization, dehomogenization

mutable struct MPolyRing_dec{T, S} <: AbstractAlgebra.MPolyRing{T}
  R::S
  D::GrpAbFinGen
  d::Vector{GrpAbFinGenElem}
  lt
  Hecke.@declare_other
  function MPolyRing_dec(R::S, d::Array{GrpAbFinGenElem, 1}) where {S}
    r = new{elem_type(base_ring(R)), S}()
    r.R = R
    r.D = parent(d[1])
    r.d = d
    return r
  end
  function MPolyRing_dec(R::S, d::Array{GrpAbFinGenElem, 1}, lt) where {S}
    r = new{elem_type(base_ring(R)), S}()
    r.R = R
    r.D = parent(d[1])
    r.d = d
    r.lt = lt
    return r
  end
end

isgraded(W::MPolyRing_dec) = !isdefined(W, :lt)
isfiltered(W::MPolyRing_dec) = isdefined(W, :lt)

function show(io::IO, W::MPolyRing_dec)
  Hecke.@show_name(io, W)
  Hecke.@show_special(io, W)
  if isfiltered(W)
    println(io, "$(W.R) filtrated by ")
  else
    println(io, "$(W.R) graded by ")
  end
  R = W.R
  g = gens(R)
  for i = 1:ngens(R)
    if i == ngens(R)
       print(io, "\t$(g[i]) -> $(W.d[i].coeff)")
    else  
       println(io, "\t$(g[i]) -> $(W.d[i].coeff)")
    end  
  end
#  println(IOContext(io, :compact => true, ), W.d)
end


function decorate(R::MPolyRing)
  A = abelian_group([0])
  S = MPolyRing_dec(R, [1*A[1] for i = 1: ngens(R)], (x,y) -> x[1] < y[1])
  return S, map(R, gens(R))
end

@doc Markdown.doc"""
    grade(R::MPolyRing, v::Array{Int, 1})

Grade `R` by assigning weights to the variables according to the entries of `v`. 

    grade(R::MPolyRing)

Grade `R` by assigning weight 1 to each variable. 

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> v = [1, 2, 3]
3-element Array{Int64,1}:
 1
 2
 3

julia> grade(R, v)
Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
	x -> [1]
	y -> [2]
	z -> [3]


julia> S = grade(R)
Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
	x -> [1]
	y -> [1]
	z -> [1]


julia> typeof(x)
fmpq_mpoly

julia> typeof(S(x))
Oscar.MPolyElem_dec{fmpq}
```
"""
function grade(R::MPolyRing, v::Array{Int, 1})
  A = abelian_group([0])
  Hecke.set_special(A, :show_elem => show_special_elem_grad) 
  return MPolyRing_dec(R, [i*A[1] for i = v])
end
function grade(R::MPolyRing)
  A = abelian_group([0])
  S = MPolyRing_dec(R, [1*A[1] for i = 1: ngens(R)])
  return S, map(S, gens(R))
end

filtrate(R::MPolyRing) = decorate(R)

function show_special_elem_grad(io::IO, a::GrpAbFinGenElem)
  if get(io, :compact, false)
    print(io, a.coeff)
  else
    print(io, "graded by $(a.coeff)")
  end
end

function filtrate(R::MPolyRing, v::Array{Int, 1})
  A = abelian_group([0])
  Hecke.set_special(A, :show_elem => show_special_elem_grad) 
  S = MPolyRing_dec(R, [i*A[1] for i = v], (x,y) -> x[1] < y[1])
  return S, map(S, gens(R))
end

function grade(R::MPolyRing, v::Array{Int, 1})
  A = abelian_group([0])
  Hecke.set_special(A, :show_elem => show_special_elem_grad) 
  S = MPolyRing_dec(R, [i*A[1] for i = v])
  return S, map(S, gens(R))
end

function filtrate(R::MPolyRing, v::Array{GrpAbFinGenElem, 1}, lt)
  S = MPolyRing_dec(R, v, lt)
  return S, map(S, gens(R))
end

function grade(R::MPolyRing, v::Array{GrpAbFinGenElem, 1})
  S = MPolyRing_dec(R, v)
  return S, map(S, gens(R))
end

struct MPolyElem_dec{T, S} <: MPolyElem{T}
  f::S
  parent
  function MPolyElem_dec(f::S, p) where {S}
    r = new{elem_type(base_ring(f)), S}(f, p)
#    if isgraded(p) && length(r) > 1
#      if !ishomogenous(r)
#        error("element not homogenous")
#      end
#both wrong and undesired.
#    end
    return r
  end
end

function show(io::IO, w::MPolyElem_dec)
  show(io, w.f)
end

parent(a::MPolyElem_dec{T, S}) where {T, S} = a.parent::MPolyRing_dec{T, parent_type(S)}

Nemo.symbols(R::MPolyRing_dec) = symbols(R.R)
Nemo.nvars(R::MPolyRing_dec) = nvars(R.R)

elem_type(::MPolyRing_dec{T, S}) where {T, S} = MPolyElem_dec{T, elem_type(S)}
elem_type(::Type{MPolyRing_dec{T, S}}) where {T, S} = MPolyElem_dec{T, elem_type(S)}
parent_type(::Type{MPolyElem_dec{T, S}}) where {T, S} = MPolyRing_dec{T, parent_type(S)}

(W::MPolyRing_dec)() = MPolyElem_dec(W.R(), W)
(W::MPolyRing_dec)(i::Int) = MPolyElem_dec(W.R(i), W)
(W::MPolyRing_dec)(i::RingElem) = MPolyElem_dec(W.R(i), W)
(W::MPolyRing_dec)(f::Singular.spoly) = MPolyElem_dec(W.R(f), W)
(W::MPolyRing_dec)(f::MPolyElem) = MPolyElem_dec(f, W)
(W::MPolyRing_dec)(g::MPolyElem_dec) = MPolyElem_dec(g.f, W)
one(W::MPolyRing_dec) = MPolyElem_dec(one(W.R), W)
zero(W::MPolyRing_dec) = MPolyElem_dec(zero(W.R), W)

################################################################################
#
#  Binary operations
#
################################################################################

for T in [:(+), :(-), :(*), :divexact]
  @eval ($T)(a::MPolyElem_dec,
             b::MPolyElem_dec) = MPolyElem_dec($T(a.f, b.f), parent(a))
end

################################################################################
#
#  Unitary operations
#
################################################################################

-(a::MPolyElem_dec)   = MPolyElem_dec(-a.f, parent(a))

################################################################################
#
#  Binary ad hoc operations
#
################################################################################

divexact(a::MPolyElem_dec, b::RingElem) = MPolyElem_dec(divexact(a.f, b), parent(a))

divexact(a::MPolyElem_dec, b::Integer) = MPolyElem_dec(divexact(a.f, b), parent(a))

divexact(a::MPolyElem_dec, b::Rational) = MPolyElem_dec(divexact(a.f, b), parent(a))

for T in [:(-), :(+)]
  @eval ($T)(a::MPolyElem_dec,
             b::RingElem) = MPolyElem_dec($(T)(a.poly, b), parent(a))

  @eval ($T)(a::MPolyElem_dec,
             b::Integer) = MPolyElem_dec($(T)(a.poly, b), parent(a))

  @eval ($T)(a::MPolyElem_dec,
             b::Rational) = MPolyElem_dec($(T)(a.poly, b), parent(a))

  @eval ($T)(a::RingElem,
             b::MPolyElem_dec) = MPolyElem_dec($(T)(a, b.poly), b.parent)

  @eval ($T)(a::Integer,
             b::MPolyElem_dec) = MPolyElem_dec($(T)(a, b.poly), b.parent)

  @eval ($T)(a::Rational,
             b::MPolyElem_dec) = MPolyElem_dec($(T)(a, b.poly), b.parent)
end

################################################################################
#
#  Equality
#
################################################################################

function factor(x::Oscar.MPolyElem_dec)
  R = parent(x)
  D = Dict{elem_type(R), Int64}()
  F = factor(x.f)
  n=length(F.fac)
  #if n == 1
  #  return Fac(R(F.unit), D)
  #else
    for i in keys(F.fac)
     push!(D, R(i) => Int64(F[i]))
    end
  return Fac(R(F.unit), D)
  #end
end


function gcd(x::Oscar.MPolyElem_dec, y::Oscar.MPolyElem_dec)
  R = parent(x)
  return R(gcd(x.f, y.f))
end

function div(x::Oscar.MPolyElem_dec, y::Oscar.MPolyElem_dec)
  R = parent(x)
  return R(div(x.f, y.f))
end


==(a::MPolyElem_dec, b::MPolyElem_dec) = a.f == b.f

^(a::MPolyElem_dec, i::Int) = MPolyElem_dec(a.f^i, parent(a))

function mul!(a::MPolyElem_dec, b::MPolyElem_dec, c::MPolyElem_dec)
  return b*c
end

function addeq!(a::MPolyElem_dec, b::MPolyElem_dec)
  return a+b
end

length(a::MPolyElem_dec) = length(a.f)
monomial(a::MPolyElem_dec, i::Int) = parent(a)(monomial(a.f, i))
coeff(a::MPolyElem_dec, i::Int) = coeff(a.f, i)

function singular_ring(R::MPolyRing_dec; keep_ordering::Bool = false)
  return singular_ring(R.R, keep_ordering = keep_ordering)
end

MPolyCoeffs(f::MPolyElem_dec) = MPolyCoeffs(f.f)
MPolyExponentVectors(f::MPolyElem_dec) = MPolyExponentVectors(f.f)

function push_term!(M::MPolyBuildCtx{<:MPolyElem_dec{T, S}}, c::T, expv::Vector{Int}) where {T <: RingElement, S}
  if iszero(c)
    return M
  end
  len = length(M.poly.f) + 1
  set_exponent_vector!(M.poly.f, len, expv)
  setcoeff!(M.poly.f, len, c)
  return M
end

function finish(M::MPolyBuildCtx{<:MPolyElem_dec})
  f = sort_terms!(M.poly.f)
  f = combine_like_terms!(M.poly.f)
  return parent(M.poly)(f)
end

function ideal(g::Array{T, 1}) where {T <: MPolyElem_dec}
  if isgraded(parent(g[1]))
    @assert all(ishomogenous, g)
  end
  return MPolyIdeal(g)
end

function jacobi_matrix(f::MPolyElem_dec)
  R = parent(f)
  n = nvars(R)
  return matrix(R, n, 1, [derivative(f, i) for i=1:n])
end

function jacobi_ideal(f::MPolyElem_dec)
  R = parent(f)
  n = nvars(R)
  return ideal(R, [derivative(f, i) for i=1:n])
end

function jacobi_matrix(g::Array{<:MPolyElem_dec, 1})
  R = parent(g[1])
  n = nvars(R)
  @assert all(x->parent(x) == R, g)
  return matrix(R, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

function degree(a::MPolyElem_dec)
  W = parent(a)
  w = W.D[0]
  first = true
  d = W.d
  for c = MPolyExponentVectors(a.f)
    u = W.D[0]
    for i=1:length(c)
      u += c[i]*d[i]
    end
    if isfiltered(W)
      w = W.lt(w, u) ? u : w
    elseif first
      first = false
      w = u
    else
      w == u || error("element not homogenous")
    end
  end
  return w
end

function ishomogenous(a::MPolyElem_dec)
  D = parent(a).D
  d = parent(a).d
  S = Set{elem_type(D)}()
  for c = MPolyExponentVectors(a.f)
    u = parent(a).D[0]
    for i=1:length(c)
      u += c[i]*d[i]
    end
    push!(S, u)
    if length(S) > 1
      return false
    end
  end
  return true
end

function homogenous_components(a::MPolyElem_dec{T, S}) where {T, S}
  D = parent(a).D
  d = parent(a).d
  h = Dict{elem_type(D), typeof(a)}()
  W = parent(a)
  R = W.R
  # First assemble the homogenous components into the build contexts.
  # Afterwards compute the polynomials.
  hh = Dict{elem_type(D), MPolyBuildCtx{S, DataType}}()
  dmat = vcat([d[i].coeff for i in 1:length(d)])
  tmat = zero_matrix(ZZ, 1, nvars(R))
  res_mat = zero_matrix(ZZ, 1, ncols(dmat))
  for (c, e) = Base.Iterators.zip(coefficients(a.f), exponent_vectors(a.f))
    # this is non-allocating
    for i in 1:length(e)
      tmat[1, i] = e[i]
    end
    mul!(res_mat, tmat, dmat)
    u = GrpAbFinGenElem(D, res_mat)
    if haskey(hh, u)
      ctx = hh[u]
      push_term!(ctx, c, e)
    else
      # We put u in the dictionary
      # Make a fresh res_mat, which can be used the for the next u
      res_mat = deepcopy(res_mat)
      ctx = MPolyBuildCtx(R)
      push_term!(ctx, c, e)
      hh[u] = ctx
    end
  end
  hhh = Dict{elem_type(D), typeof(a)}()
  for (u, C) in hh
    hhh[u] = W(finish(C))
  end

  return hhh
end

function homogenous_component(a::MPolyElem_dec, g::GrpAbFinGenElem)
  R = parent(a).R
  r = R(0)
  d = parent(a).d
  for (c, m) = Base.Iterators.zip(MPolyCoeffs(a.f), Generic.MPolyMonomials(a.f))
    e = exponent_vector(m, 1)
    u = parent(a).D[0]
    for i=1:length(e)
      u += e[i]*d[i]
    end
    if u == g
      r += c*m
    end
  end
  return parent(a)(r)
end

function degree(a::MPolyQuoElem{<:MPolyElem_dec})
  simplify!(a)
  return degree(a.f)
end

isfiltered(q::MPolyQuo) = isfiltered(q.R)
isgraded(q::MPolyQuo) = isgraded(q.R)

function homogenous_component(a::MPolyQuoElem{<:MPolyElem_dec}, d::GrpAbFinGenElem)
  simplify!(a)
  return homogenous_component(a.f, d)
end

function homogenous_components(a::MPolyQuoElem{<:MPolyElem_dec})
  simplify!(a)
  return homogenous_components(a.f)
end

function ishomogenous(a::MPolyQuoElem{<:MPolyElem_dec})
  simplify!(a)
  return ishomogenous(a.f)
end

decoration(q::MPolyQuo{<:MPolyElem_dec}) = decoration(q.R)

base_ring(W::MPolyRing_dec) = base_ring(W.R)
Nemo.ngens(W::MPolyRing_dec) = Nemo.ngens(W.R)
Nemo.ngens(R::MPolyRing) = Nemo.nvars(R)
Nemo.gens(W::MPolyRing_dec) = map(W, gens(W.R))
Nemo.gen(W::MPolyRing_dec, i::Int) = W(gen(W.R, i))
Base.getindex(W::MPolyRing_dec, i::Int) = W(W.R[i])

base_ring(f::MPolyElem_dec) = base_ring(f.f)

*(r::fmpq, w::MPolyElem_dec) = parent(w)(r*w.f)

function show_homo_comp(io::IO, M)
  (W, d) = Hecke.get_special(M, :data)
  n = Hecke.get_special(W, :name)
  if n != nothing
    print(io, "$(n)_$(d.coeff) of dim $(dim(M))")
  else
    println(io, "homogenous component of $W of degree $d")
  end
end

function homogenous_component(W::MPolyRing_dec, d::GrpAbFinGenElem)
  #TODO: lazy: ie. no enumeration of points
  #      aparently it is possible to get the number of points faster than the points
  D = W.D
  h = hom(free_abelian_group(ngens(W)), W.d)
  fl, p = haspreimage(h, d)
  R = base_ring(W)
  @assert fl
  k, im = kernel(h)
  #need the positive elements in there...
  #Ax = b, Cx >= 0
  C = identity_matrix(FlintZZ, ngens(W))
  A = vcat([x.coeff for x = W.d])
  k = solve_mixed(A', d.coeff', C)
  B = elem_type(W)[]
  for ee = 1:nrows(k)
    e = k[ee, :]
    a = MPolyBuildCtx(W.R)
    push_term!(a, R(1), [Int(e[i]) for i in 1:length(e)])
    push!(B, W(finish(a)))
  end
  M, h = vector_space(R, B, target = W)
  Hecke.set_special(M, :show => show_homo_comp, :data => (W, d))
  add_relshp(M, W, x -> sum(x[i] * B[i] for i=1:length(B)))
#  add_relshp(W, M, g)
  return M, h
end

base_ring(W::MPolyQuo) = W.R
modulus(W::MPolyQuo) = W.I

function hash(w::MPolyQuoElem, u::UInt)
  simplify!(w)
  return hash(w.f, u)
end

function homogenous_component(W::MPolyQuo{<:MPolyElem_dec}, d::GrpAbFinGenElem)
  #TODO: lazy: ie. no enumeration of points
  #      aparently it is possible to get the number of points faster than the points
  D = parent(d)
  @assert D == decoration(W)
  R = base_ring(W)
  I = modulus(W)

  H, mH = homogenous_component(R, d)
  B = Set{elem_type(W)}()
  for h = basis(H)
    b = W(mH(h))
    if !iszero(b)
      push!(B, b)
    end
  end
  B = [x for x = B]

  M, h = vector_space(base_ring(R), B, target = W)
  Hecke.set_special(M, :show => show_homo_comp, :data => (W, d))
  return M, h
end

function vector_space(K::AbstractAlgebra.Field, e::Array{T, 1}; target = nothing) where {T <:MPolyElem}
  local R
  if length(e) == 0
    R = target
    @assert R !== nothing
  else
    R = parent(e[1])
  end
  @assert base_ring(R) == K
  mon = Dict{elem_type(R), Int}()
  mon_idx = Array{elem_type(R), 1}()
  M = sparse_matrix(K)
  last_pos = 1
  for i = e
    pos = Array{Int, 1}()
    val = Array{elem_type(K), 1}()
    for (c, m) = Base.Iterators.zip(coefficients(i), monomials(i))
      if haskey(mon, m)
        push!(pos, mon[m])
        push!(val, c)
      else
        push!(mon_idx, m)
        mon[m] = last_pos
        push!(pos, last_pos)
        push!(val, c)
        last_pos += 1
      end
    end
    push!(M, sparse_row(K, pos, val))
  end
  Hecke.echelon!(M, complete = true)

  b = Array{elem_type(R), 1}()
  for i=1:nrows(M)
    s = zero(e[1])
    for (k,v) = M[i]
      s += v*mon_idx[k]
    end
    push!(b, s)
  end

  F = FreeModule(K, length(b), cached = false)
  function g(x::T)
    @assert parent(x) == R
    v = zero(F)
    for (c, m) = Base.Iterators.zip(coefficients(x), monomials(x))
      if !haskey(mon, m)
        error("not in image")
      end
      v += c*gen(F, mon[m])
    end
    return v
  end
  h = MapFromFunc(x -> sum(x[i] * b[i] for i=1:length(b)), g, F, R)

  return F, h
end

###########################################
# needs re-thought
function (W::MPolyRing_dec)(m::Generic.FreeModuleElem) 
  h = hasrelshp(parent(m), W)
  if h !== nothing
    return h(m)
  end
  error("no coercion possible")
end

#########################################
function add_relshp(R, S, h)
  #this assumes that h is essentially a canonical map from R -> S
  if !isdefined(R, :other)
    R.other = dict{Symbol, Any}()
  end
  if !haskey(R.other, :relshp)
    R.other[:relshp] = Dict{Any, Any}()
  end
  if haskey(R.other[:relshp], S)
    error("try to add double")
  end
  R.other[:relshp][S] = h
end

function hasrelshp(R, S)
  r = Hecke.get_special(R, :relshp)
  if r === nothing
    return r
  end
  if haskey(r, S)
    return r[S]
  end
  #now the hard bit: traverse the graph, not falling into cycles.
end

############################################################################
############################################################################
############################################################################

function sing_hilb(I::Singular.sideal)
  a = Array{Int32, 1}()
  @assert I.isGB
  Singular.libSingular.scHilb(I.ptr, base_ring(I).ptr, a)
  return a
end

mutable struct HilbertData
  data::Array{Int32, 1}
  I::MPolyIdeal
  function HilbertData(I::MPolyIdeal)
    Oscar.groebner_assure(I)
    h = sing_hilb(I.gb.S)
    return new(h, I)
  end
  function HilbertData(B::BiPolyArray)
    return HilbertData(Oscar.MPolyIdeal(B))
  end
end

function hilbert_series(H::HilbertData, i::Int= 1)
  Zt, t = ZZ["t"]
  if i==1
    return Zt(map(fmpz, H.data[1:end-1])), (1-gen(Zt))^(ngens(base_ring(H.I)))
  elseif i==2
    h = hilbert_series(H, 1)[1]
    return divexact(h, (1-gen(Zt))^(ngens(base_ring(H.I))-dim(H.I))), (1-gen(Zt))^dim(H.I)
  end
  error("2nd parameter must be 1 or 2")
end

#Decker-Lossen, p23/24
function hilbert_polynomial(H::HilbertData)
  q, dn = hilbert_series(H, 2)
  a = fmpq[]
  nf = fmpq(1)
  d = degree(dn)-1
  for i=0:d
    push!(a, q(1)//nf)
    if i>0
      nf *= i
    end
    q = derivative(q)
  end
  Qt, t = QQ["t"]
  t = gen(Qt)
  bin = one(parent(t))
  b = fmpq_poly[]
  if d==-1 return zero(parent(t)) end
  for i=0:d
    push!(b, (-1)^(d-i)*a[d-i+1]*bin)
    bin *= (t+i+1)*fmpq(1, i+1)
  end
  return sum(b)
end

function Oscar.degree(H::HilbertData)
  P = hilbert_polynomial(H)
  if P==zero(parent(P))
     q, _ = hilbert_series(H, 2)
     return q(1)
  end
  return leading_coefficient(P)*factorial(degree(P))
end

function (P::FmpqRelSeriesRing)(H::HilbertData)
  n = hilbert_series(H, 1)[1]
  d = (1-gen(parent(n)))^ngens(base_ring(H.I))
  g = gcd(n, d)
  n = divexact(n, g)
  d = divexact(d, g)
  Qt, t = QQ["t"]
  nn = map_coefficients(QQ, n, parent = Qt)
  dd = map_coefficients(QQ, d, parent = Qt)
  gg, ee, _ = gcdx(dd, gen(Qt)^max_precision(P))
  @assert isone(gg)
  nn = Hecke.mullow(nn, ee, max_precision(P)+1)
  c = collect(coefficients(nn))
  return P(map(fmpq, c), length(c), max_precision(P), 0)
end

function hilbert_series_expanded(H::HilbertData, d::Int)
   T, t = PowerSeriesRing(QQ, d, "t")   
   return T(H)
end

function hilbert_function(H::HilbertData, d::Int)
   HS = hilbert_series_expanded(H,d)
   return coeff(hilbert_series_expanded(H, d), d)
end

function Base.show(io::IO, h::HilbertData)
  print(io, "Hilbert Series for $(h.I), data: $(h.data)")
end

############################################################################
### Homogenization and Dehomogenization
############################################################################

function homogenization(f::MPolyElem, S::MPolyRing_dec, pos::Int = 1)
  d = total_degree(f)
  B = MPolyBuildCtx(S)
  for (c,e) = zip(coefficients(f), exponent_vectors(f))
    insert!(e, pos, d-sum(e))
    push_term!(B, c, e)
  end
  return finish(B)
end

@doc Markdown.doc"""
    homogenization(f::MPolyElem, var::String, pos::Int = 1)

    homogenization(V::Vector{T}, var::String, pos::Int = 1) where {T <: MPolyElem}

    homogenization(I::MPolyIdeal{T}, var::String, pos::Int = 1; ordering::Symbol = :degrevlex) where {T <: MPolyElem}

Return the homogenization of `f`, `V`, or `I` in a graded ring with additional variable `var` at position `pos`.

CAVEAT: Homogenizing an ideal requires a Gröbner basis computation. This may take some time.
"""
function homogenization(f::MPolyElem, var::String, pos::Int = 1)
  R = parent(f)
  A = String.(symbols(R))
  l = length(A)
  if (pos > l+1) ||  (pos <1)
      throw(ArgumentError("Index out of range."))
  end
  insert!(A, pos, var)
  L, _ = PolynomialRing(R.base_ring, A)
  S = grade(L)
  return homogenization(f, S, pos)
end
function homogenization(V::Vector{T}, var::String, pos::Int = 1) where {T <: MPolyElem}
  @assert all(x->parent(x) == parent(V[1]), V)
  R = parent(V[1])
  A = String.(symbols(R))
  l = length(A)
  if (pos > l+1) ||  (pos <1)
      throw(ArgumentError("Index out of range."))
  end
  insert!(A, pos, var)
  L, _ = PolynomialRing(R.base_ring, A)
  S = grade(L)
  l = length(V)
  return [homogenization(V[i], S, pos) for i=1:l]
end
function homogenization(I::MPolyIdeal{T}, var::String, pos::Int = 1; ordering::Symbol = :degrevlex) where {T <: MPolyElem}
  return ideal(homogenization(groebner_basis(I, ordering), var, pos))
end

function dehomogenization(F::MPolyElem_dec, R::MPolyRing, pos::Int)
  B = MPolyBuildCtx(R)
  for (c,e) = zip(coefficients(F), exponent_vectors(F))
    deleteat!(e, pos)
    push_term!(B, c, e)
  end
  return finish(B)
end

@doc Markdown.doc"""
    dehomogenization(F::MPolyElem_dec, pos::Int)

    dehomogenization(V::Vector{T}, pos::Int) where {T <: MPolyElem_dec}

    dehomogenization(I::MPolyIdeal{T}, pos::Int) where {T <: MPolyElem_dec}

Return the dehomogenization of `F`, `V`, or `I` in a ring not depending on the variable at position `pos`.
"""
function dehomogenization(F::MPolyElem_dec, pos::Int)
  S = parent(F)
  A = String.(symbols(S))
  l = length(A)
  if (pos > l+1) ||  (pos <1)
      throw(ArgumentError("Index out of range."))
  end
  deleteat!(A, pos)
  R, _ = PolynomialRing(base_ring(S), A)
  return dehomogenization(F, R, pos)
end
function dehomogenization(V::Vector{T}, pos::Int) where {T <: MPolyElem_dec}
  @assert all(x->parent(x) == parent(V[1]), V)
  S = parent(V[1])
  A = String.(symbols(S))
  l = length(A)
  if (pos > l+1) ||  (pos <1)
      throw(ArgumentError("Index out of range."))
  end
  deleteat!(A, pos)
  R, _ = PolynomialRing(base_ring(S), A)
  l = length(V)
  return [dehomogenization(V[i], R, pos) for i=1:l]
end
function dehomogenization(I::MPolyIdeal{T}, pos::Int) where {T <: MPolyElem_dec}
  return ideal(dehomogenization(gens(I), pos))
end


