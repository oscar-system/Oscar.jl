export weight, decorate, ishomogenous, homogenous_components, filtrate,
grade, homogenous_component, jacobi_matrix, jacobi_ideal

mutable struct MPolyRing_dec{T} <: AbstractAlgebra.MPolyRing{T}
  R::MPolyRing{T}
  D::GrpAbFinGen
  d::Array{GrpAbFinGenElem}
  lt
  Hecke.@declare_other
  function MPolyRing_dec(R :: MPolyRing{S}, d::Array{GrpAbFinGenElem, 1}) where {S}
    r = new{S}()
    r.R = R
    r.D = parent(d[1])
    r.d = d
    return r
  end
  function MPolyRing_dec(R::MPolyRing{T}, d::Array{GrpAbFinGenElem, 1}, lt) where {T}
    r = new{T}()
    r.R = R
    r.D = parent(d[1])
    r.d = d
    r.lt = lt
    return r
  end
end

isgraded(W::MPolyRing_dec) = !isdefined(W, :lt)
isfiltrated(W::MPolyRing_dec) = isdefined(W, :lt)

function show(io::IO, W::MPolyRing_dec)
  Hecke.@show_name(io, W)
  Hecke.@show_special(io, W)
  if isfiltrated(W)
    println(io, "$(W.R) filtrated by ")
  else
    println(io, "$(W.R) graded by ")
  end
  R = W.R
  g = gens(R)
  for i = 1:ngens(R)
    println(io, "\t$(g[i]) -> $(W.d[i].coeff)")
  end
#  println(IOContext(io, :compact => true, ), W.d)
end

function decorate(R::MPolyRing)
  A = abelian_group([0])
  return MPolyRing_dec(R, [1*A[1] for i = 1: ngens(R)], (x,y) -> x[1] < y[1])
end

grade(R::MPolyRing) = decorate(R)
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
  return MPolyRing_dec(R, [i*A[1] for i = v], (x,y) -> x[1] < y[1])
end

function grade(R::MPolyRing, v::Array{Int, 1})
  A = abelian_group([0])
  Hecke.set_special(A, :show_elem => show_special_elem_grad) 
  return MPolyRing_dec(R, [i*A[1] for i = v])
end

function filtrate(R::MPolyRing, v::Array{GrpAbFinGenElem, 1}, lt)
  return MPolyRing_dec(R, v, lt)
end

function grade(R::MPolyRing, v::Array{GrpAbFinGenElem, 1})
  return MPolyRing_dec(R, v)
end

struct MPolyElem_dec{T} <: MPolyElem{T}
  f::MPolyElem{T}
  parent::MPolyRing_dec{T}
  function MPolyElem_dec(f::MPolyElem{T}, p) where {T}
    r = new{T}(f, p)
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

Nemo.symbols(R::MPolyRing_dec) = symbols(R.R)
Nemo.nvars(R::MPolyRing_dec) = nvars(R.R)

elem_type(::MPolyRing_dec{T}) where {T} = MPolyElem_dec{T}
elem_type(::Type{MPolyRing_dec{T}}) where {T} = MPolyElem_dec{T}
parent_type(::Type{MPolyElem_dec{T}}) where {T} = MPolyRing_dec{T}

(W::MPolyRing_dec)() = MPolyElem_dec(W.R(), W)
(W::MPolyRing_dec)(i::Int) = MPolyElem_dec(W.R(i), W)
(W::MPolyRing_dec)(i::RingElem) = MPolyElem_dec(W.R(i), W)
(W::MPolyRing_dec)(f::MPolyElem) = MPolyElem_dec(f, W)
(W::MPolyRing_dec)(g::MPolyElem_dec) = MPolyElem_dec(g.f, W)
one(W::MPolyRing_dec) = MPolyElem_dec(one(W.R), W)
zero(W::MPolyRing_dec) = MPolyElem_dec(zero(W.R), W)

+(a::MPolyElem_dec, b::MPolyElem_dec)   = MPolyElem_dec(a.f+b.f, a.parent)
-(a::MPolyElem_dec, b::MPolyElem_dec)   = MPolyElem_dec(a.f-b.f, a.parent)
-(a::MPolyElem_dec)   = MPolyElem_dec(-a.f, a.parent)
*(a::MPolyElem_dec, b::MPolyElem_dec)   = MPolyElem_dec(a.f*b.f, a.parent)
==(a::MPolyElem_dec, b::MPolyElem_dec)   = a.f == b.f
^(a::MPolyElem_dec, i::Int)    = MPolyElem_dec(a.f^i, a.parent)

function Oscar.mul!(a::MPolyElem_dec, b::MPolyElem_dec, c::MPolyElem_dec)
  return b*c
end

function Oscar.addeq!(a::MPolyElem_dec, b::MPolyElem_dec)
  return a+b
end

parent(a::MPolyElem_dec) = a.parent
length(a::MPolyElem_dec) = length(a.f)
monomial(a::MPolyElem_dec, i::Int) = parent(a)(monomial(a.f, i))
coeff(a::MPolyElem_dec, i::Int) = coeff(a.f, i)

function singular_ring(R::MPolyRing_dec; keep_ordering::Bool = false)
  return singular_ring(R.R, keep_ordering = keep_ordering)
end

MPolyCoeffs(f::MPolyElem_dec) = MPolyCoeffs(f.f)
MPolyExponentVectors(f::MPolyElem_dec) = MPolyExponentVectors(f.f)
function push_term!(M::MPolyBuildCtx{<:MPolyElem_dec{S}}, c::S, expv::Vector{Int}) where S <: RingElement  
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
    if isfiltrated(W)
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

function homogenous_components(a::MPolyElem_dec)
  D = parent(a).D
  d = parent(a).d
  h = Dict{elem_type(D), typeof(a)}()
  W = parent(a)
  for (c, m) = Base.Iterators.zip(MPolyCoeffs(a.f), Generic.MPolyMonomials(a.f))
    e = exponent_vector(m, 1)
    u = parent(a).D[0]
    for i=1:length(e)
      u += e[i]*d[i]
    end
    if haskey(h, u)
      h[u] += W(c*m)
    else
      h[u] = W(c*m)
    end
  end
  return h
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

isfiltrated(q::MPolyQuo) = isfiltrated(q.R)
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
    push_term!(a, R(1), [Int(x) for x = e])
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


