export weight, decorate, ishomogenous, homogenous_components, filtrate, grade, 
       homogenous_component

mutable struct MPolyRing_dec{T} <: AbstractAlgebra.Ring
  R::T
  D::GrpAbFinGen
  d::Array{GrpAbFinGenElem}
  lt
  Hecke.@declare_other
  function MPolyRing_dec(R::T, d::Array{GrpAbFinGenElem, 1}) where {T}
    r = new{T}()
    r.R = R
    r.D = parent(d[1])
    r.d = d
    return r
  end
  function MPolyRing_dec(R::T, d::Array{GrpAbFinGenElem, 1}, lt) where {T}
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
  A = AbelianGroup([0])
  return MPolyRing_dec(R, [0*A[1] for i = v], (x,y) -> x[1] < y[1])
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
  A = AbelianGroup([0])
  return MPolyRing_dec(R, [i*A[1] for i = v], (x,y) -> x[1] < y[1])
end

function grade(R::MPolyRing, v::Array{Int, 1})
  A = AbelianGroup([0])
  Hecke.set_special(A, :show_elem => show_special_elem_grad) 
  return MPolyRing_dec(R, [i*A[1] for i = v])
end

function filtrate(R::MPolyRing, v::Array{GrpAbFinGenElem, 1}, lt)
  return MPolyRing_dec(R, v, lt)
end

function grade(R::MPolyRing, v::Array{GrpAbFinGenElem, 1})
  return MPolyRing_dec(R, v)
end

struct MPolyElem_dec{T} <: AbstractAlgebra.RingElem
  f::T
  parent::AbstractAlgebra.Ring
  function MPolyElem_dec(f::T, p) where {T}
    r = new{T}(f, p)
    if isgraded(p) && length(r) > 1
      if !ishomogenous(r)
        error("element not homogenous")
      end
    end
    return r
  end
end

function show(io::IO, w::MPolyElem_dec)
  show(io, w.f)
end

elem_type(::MPolyRing_dec{T}) where {T} = MPolyElem_dec{elem_type(T)}
elem_type(::Type{MPolyRing_dec{T}}) where {T} = MPolyElem_dec{elem_type(T)}
parent_type(::Type{MPolyElem_dec{T}}) where {T} = MPolyRing_dec{parent_type(T)}

(W::MPolyRing_dec)(i::Int) = MPolyElem_dec(W.R(i), W)
(W::MPolyRing_dec)(i::RingElem) = MPolyElem_dec(W.R(i), W)
(W::MPolyRing_dec)(f::MPolyElem) = MPolyElem_dec(f, W)
(W::MPolyRing_dec)(g::MPolyElem_dec) = MPolyElem_dec(g.f, W)
one(W::MPolyRing_dec) = MPolyElem_dec(one(W.R), W)
zero(W::MPolyRing_dec) = MPolyElem_dec(zero(W.R), W)

+(a::MPolyElem_dec{T}, b::MPolyElem_dec{T}) where {T} = MPolyElem_dec(a.f+b.f, a.parent)
*(a::MPolyElem_dec{T}, b::MPolyElem_dec{T}) where {T} = MPolyElem_dec(a.f*b.f, a.parent)
==(a::MPolyElem_dec{T}, b::MPolyElem_dec{T}) where {T} = a.f == b.f
^(a::MPolyElem_dec{T}, i::Int)  where {T} = MPolyElem_dec(a.f^i, a.parent)

parent(a::MPolyElem_dec) = a.parent
length(a::MPolyElem_dec) = length(a.f)

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

function homogenous_component(a::MPolyElem_dec, d::GrpAbFinGenElem)
  R = parent(a).R
  r = R(0)
  for (c, m) = Base.Iterators.zip(MPolyCoeffs(a.f), Generic.MPolyMonomials(a.f))
    e = exponent_vector(m, 1)
    u = parent(a).D[0]
    for i=1:length(e)
      u += e[i]*d[i]
    end
    if u == d
      r += c*m
    end
  end
  return W(r)
end

base_ring(W::MPolyRing_dec) = base_ring(W.R)
Nemo.ngens(W::MPolyRing_dec) = Nemo.ngens(W.R)
Nemo.ngens(R::MPolyRing) = Nemo.nvars(R)
Nemo.gens(W::MPolyRing_dec) = map(W, gens(W.R))

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
  h = hom(FreeAbelianGroup(ngens(W)), W.d)
  fl, p = haspreimage(h, d)
  R = base_ring(W)
  @assert fl
  k, im = kernel(h)
  #need the positive elements in there...
  #Ax = b, Cx >= 0
  C = identity_matrix(FlintZZ, ngens(W))
  A = vcat([x.coeff for x = W.d])
  k = Hecke.solve_mixed(A', d.coeff', C)
  B = []
  for e = k
    a = MPolyBuildCtx(W.R)
    push_term!(a, R(1), [Int(x) for x = e])
    push!(B, W(finish(a)))
  end
  M = FreeModule(R, length(B), cached = false)
  g = function(y::MPolyElem_dec)
    x = elem_type(R)[R(0) for i=1:length(B)]
    for (c, m) = Base.Iterators.zip(MPolyCoeffs(y.f), Generic.MPolyMonomials(y.f))
      mm = W(m)
      p = Base.findfirst(i -> B[i] == mm, 1:length(B))
      if p === nothing
        error("element not in this component")
      end
      x[p] = c
    end
    return M(x)
  end
  h = MapFromFunc(x -> sum(x[i] * B[i] for i=1:length(B)), g, M, W)
  Hecke.set_special(M, :show => show_homo_comp, :data => (W, d))
  add_relshp(M, W, x -> sum(x[i] * B[i] for i=1:length(B)))
#  add_relshp(W, M, g)
  return M, h
end

basis(F::Generic.FreeModule) = gens(F) # needs to be in AA

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


