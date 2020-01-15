export weight, decorate, ishomogenous, homogenous_components

struct MPolyRing_dec{T} <: AbstractAlgebra.Ring
  R::T
  D::GrpAbFinGen
  d::Array{GrpAbFinGenElem}
  lt
end

function show(io::IO, W::MPolyRing_dec)
  print(io, "$(W.R) graded by ")
  println(IOContext(io, :compact => true), W.d)
end

function decorate(R::MPolyRing)
  A = AbelianGroup([0])
  return MPolyRing_dec{typeof(R)}(R, A, [0*A[1] for i = v], (x,y) -> x[1] < y[1])
end

function decorate(R::MPolyRing, v::Array{Int, 1})
  A = AbelianGroup([0])
  return MPolyRing_dec{typeof(R)}(R, A, [i*A[1] for i = v], (x,y) -> x[1] < y[1])
end

function decorate(R::MPolyRing, v::Array{GrpAbFinGenElem, 1}, lt)
  return MPolyRing_dec{typeof(R)}(R, parent(v[1]), v, lt)
end

struct MPolyElem_dec{T} <: AbstractAlgebra.RingElem
  f::T
  parent::AbstractAlgebra.Ring
end

function show(io::IO, w::MPolyElem_dec)
  show(io, w.f)
end

elem_type(::MPolyRing_dec{T}) where {T} = MPolyElem_dec{elem_type(T)}
parent_type(::MPolyElem_dec{T}) where {T} = MPolyRing_dec{parent_type(T)}

(W::MPolyRing_dec)(f::MPolyElem) = MPolyElem_dec(f, W)
(W::MPolyRing_dec)(g::MPolyElem_dec) = MPolyElem_dec(g.f, W)

+(a::MPolyElem_dec{T}, b::MPolyElem_dec{T}) where {T} = MPolyElem_dec{T}(a.f+b.f, a.parent)
*(a::MPolyElem_dec{T}, b::MPolyElem_dec{T}) where {T} = MPolyElem_dec{T}(a.f*b.f, a.parent)
^(a::MPolyElem_dec{T}, i::Int)  where {T} = MPolyElem_dec{T}(a.f^i, a.parent)

parent(a::MPolyElem_dec) = a.parent

function weight(a::MPolyElem_dec)
  w = parent(a).D[0]
  d = parent(a).d
  lt = parent(a).lt
  for c = MPolyExponentVectors(a.f)
    u = parent(a).D[0]
    for i=1:length(c)
      u += c[i]*d[i]
    end
    w = lt(w, u) ? u : w
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

