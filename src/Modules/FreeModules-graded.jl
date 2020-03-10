#TODO make d and S a function optionally - to support HUGE degree
export presentation
abstract type ModuleFP_dec{T} end
mutable struct FreeModule_dec{T} <: ModuleFP_dec{T}
  d::Array{GrpAbFinGenElem, 1}
  R::MPolyRing_dec{T}
  S::Array{Symbol, 1}
  AbstractAlgebra.@declare_other
  function FreeModule_dec(a,b,c)
    r = new{typeof(b.R)}()
    r.d = a
    r.R = b
    r.S = c
    return r
  end
end

function FreeModule(R::MPolyRing_dec, n::Int, name::String = "e"; cached::Bool = false)
  return FreeModule_dec([R.D[0] for i=1:n], R, [Symbol("$name[$i]") for i=1:n])
end

function FreeModule(R::MPolyRing_dec, d::Array{GrpAbFinGenElem, 1}, name::String = "e")
  return FreeModule_dec(d, R, [Symbol("$name[$i]") for i=1:length(d)])
end

^(R::MPolyRing_dec, n::Int) = FreeModule(R, n)

function AbstractAlgebra.extra_name(F::FreeModule_dec)
  t = Hecke.get_special(F, :twist)
  if t !== nothing
    n = Hecke.get_special(t[1], :name)
    if n !== nothing
      return "$n($(t[2]))"
    end
  end
  if length(Set(F.d)) == 1
    n = Hecke.get_special(F.R, :name)
    if n !== nothing
      return "$n^$(ngens(F))($(F.d[1]))"
    end
  end
  return nothing
end

function (F::FreeModule_dec)(a::GrpAbFinGenElem) 
  G = FreeModule(F.R, [x+a for x = F.d])
  Hecke.set_special(G, :twist => (F, a))
  return G
end

function show(io::IO, F::FreeModule_dec)
  @show_name(io, F)
  @show_special(io, F)

  print(io, "Free module of rank $(length(F.d)) over ")
  print(IOContext(io, :compact =>true), F.R)
  if isgraded(F.R)
    print(io, ", graded by\n")
  else
    print(io, ", filtrated by\n")
  end
  for i=1:dim(F)
    println(io, "\t$(F.S[i]) -> $(F.d[i])")
  end
end

dim(F::FreeModule_dec)  = length(F.d)
ngens(F::FreeModule_dec) = dim(F)

struct FreeModuleElem_dec{T}
  r::SRow{MPolyElem_dec{T}}
  parent::FreeModule_dec
end

elem_type(::Type{FreeModule_dec{T}}) where {T} = FreeModuleElem_dec{elem_type(T)}
parent_type(::Type{FreeModuleElem_dec{T}}) where {T} = FreeModule_dec{parent_type(T)}
elem_type(::FreeModule_dec{T}) where {T} = FreeModuleElem_dec{elem_type(T)}

function show(io::IO, e::FreeModuleElem_dec)
  if length(e.r) == 0
    print(io, 0)
    return
  end
  i = 1
  while i <= length(e.r)
    print(io, "(", e.r.values[i], ")*", e.parent.S[e.r.pos[i]])
    if i < length(e.r)
      print(io, " + ")
    end
    i += 1
  end
end

function basis(F::FreeModule_dec)
  bas = elem_type(F)[]
  for i=1:dim(F)
    s = Hecke.sparse_row(F.R, [(i, F.R(1))])
    push!(bas, FreeModuleElem_dec(s, F))
  end
  return bas
end
gens(F::FreeModule_dec) = basis(F)

function gen(F::FreeModule_dec, i::Int)
  @assert 0< i <= ngens(F)
  s = Hecke.sparse_row(F.R, [(i, F.R(1))])
  return FreeModuleElem_dec(s, F)
end

base_ring(F::FreeModule_dec) = F.R

function homogenous_component(F::FreeModule_dec, d::GrpAbFinGenElem)
  #TODO: lazy: ie. no enumeration of points
  #      aparently it is possible to get the number of points faster than the points
  W = base_ring(F)
  D = decoration(F)
  #have gens for W that can be combined
  #              F that can only be used
  #F ni f = sum c_i,j F[i]*w[j]
  #want deg(F[i]) + deg(w[j]) = d
  all = []
  for g = gens(F)
    Md, mMd = homogenous_component(W, d-degree(g))
    if dim(Md) > 0
      push!(all, (g, mMd))
    end
  end
  d = sum(x->dim(domain(x[1])), all)
  X = FreeModule(base_ring(W), d)
  B = []
  for (g, mMd) = all
    for x = gens(domain(mMd))
      push!(B, mMd(x)*g)
    end
  end
  function im(f)
    sum(f[i]*B[i] for i=1:dim(X))
  end
  function pr(g::FreeModuleElem_dec)
    z = F()
    for (p,v) = g.r
      i = findfirst(x->gen(F, p) == x[1], all)
      j = preimage(all[i][2], v)
      #put coeffs of j into proper spots in z
    end
    return z
  end
  return X, Hecke.MapFromFunc(im, pr, X, F)
end


-(a::FreeModuleElem_dec) = FreeModuleElem_dec(-a.r, a.parent)
-(a::FreeModuleElem_dec, b::FreeModuleElem_dec) = FreeModuleElem_dec(a.r-b.r, a.parent)
+(a::FreeModuleElem_dec, b::FreeModuleElem_dec) = FreeModuleElem_dec(a.r+b.r, a.parent)
*(a::MPolyElem_dec, b::FreeModuleElem_dec) = FreeModuleElem_dec(a*b.r, b.parent)
*(a::MPolyElem, b::FreeModuleElem_dec) = FreeModuleElem_dec(parent(b).R(a)*b.r, b.parent)
*(a::Int, b::FreeModuleElem_dec) = FreeModuleElem_dec(a*b.r, b.parent)
==(a::FreeModuleElem_dec, b::FreeModuleElem_dec) = a.r == b.r
zero(F::FreeModule_dec) = FreeModuleElem_dec(sparse_row(F.R, Tuple{Int, elem_type(F.R)}[]), F)
parent(a::FreeModuleElem_dec) = a.parent
iszero(a::FreeModuleElem_dec) = length(a.r) == 0

function degree(a::FreeModuleElem_dec)
  if iszero(a)
    error("zero has no degree")
  end
  first = true
  F = parent(a)
  W = base_ring(F)
  ww = W.D[0]
  local w
  for (p,v) = a.r
    w = degree(v)+F.d[p]
    if first
      ww = w
      first = false
    elseif isgraded(W)
      if ww != w
        error("elem not homogenous")
      end
    else
      if W.lt(ww, w) 
        ww = w
      end
    end
  end
  return w
end

function homogenous_components(a::FreeModuleElem_dec)
  res = Dict{GrpAbFinGenElem, FreeModuleElem_dec}()
  F = parent(a)
  for (p,v) = a.r
    c = homogenous_components(v)
    for (pp, vv) = c
      w = pp + F.d[p]
      if haskey(res, w)
        res[w] += vv*gen(F, p)
      else
        res[w] = vv*gen(F, p)
      end
    end
  end
  return res
end

function homogenous_component(a::FreeModuleElem_dec, g::GrpAbFinGenElem)
  F = parent(a)
  x = zero(F)
  for (p,v) = a.r
    x += homogenous_component(v, g-F.d[p])*gen(F, p)
  end
  return x
end

function ishomogenous(a::FreeModuleElem_dec)
  if iszero(a)
    return true
  end
  F = parent(a)
  first = true
  local d::GrpAbFinGenElem
  for (p,v) = a.r
    ishomogenous(v) || return false
    if first
      d = F.d[p] + degree(v)
      first = false
    else
      F.d[p] + degree(v) == d || return false
    end
  end
  return true
end

mutable struct BiModArray{T}
  O::Array{FreeModuleElem_dec{T}, 1}
  S::Singular.smodule
  F::FreeModule_dec
  SF::Singular.FreeMod
  function BiModArray(O::Array{<:FreeModuleElem_dec{T}, 1}) where {T}
    SF = singular_module(parent(O[1]))
    return BiModArray(O, SF)
  end
  function BiModArray(O::Array{<:FreeModuleElem_dec{T}, 1}, SF::Singular.FreeMod) where {T}
    r = new{T}()
    r.O = O
    r.SF = SF
    r.F = parent(O[1])
    return r
  end
  function BiModArray(F::FreeModule_dec{S}, s::Singular.smodule) where {S}
    r = new{elem_type(S)}()
    r.F = F
    r.SF = parent(s[1])
    r.S = s
    r.O = Array{FreeModuleElem_dec{elem_type(S)}, 1}(undef, Singular.ngens(s))
    return r
  end
end

function show(io::IO, F::BiModArray)
  println(io, "Array of length ", length(F))
  for i=1:length(F)
    if isassigned(F.O, i)
      println(io, i, " -> ", F.O[i])
    end
  end
  if isdefined(F, :S)
    println(io, "defined on the Singular side")
  end
end

length(F::BiModArray) = length(F.O)

function getindex(F::BiModArray, ::Val{:O}, i::Int)
  if !isassigned(F.O, i)
    F.O[i] = convert(F.F, F.S[i])
  end
  return F.O[i]
end

function getindex(F::BiModArray, ::Val{:S}, i::Int)
  if !isdefined(F, :S)
    F.S = Singular.smodule{elem_type(base_ring(F.SF))}(base_ring(F.SF), [convert(F.SF, x) for x = F.O]...)
  end
  return F.S[i]
end

function singular_assure(F::BiModArray)
  F[Val(:S), 1]
end

getindex(F::BiModArray, i::Int) = getindex(F, Val(:O), i)

function singular_module(F::FreeModule_dec)
  Sx = singular_ring(base_ring(F).R)
  return Singular.FreeModule(Sx, dim(F))
end

function convert(SF::Singular.FreeMod, m::FreeModuleElem_dec)
  g = Singular.gens(SF)
  e = 0*g[1]
  Sx = base_ring(SF)
  for (p,v) = m.r
    e += convert(Sx, v.f)*g[p]
  end
  return e
end

function convert(F::FreeModule_dec, s::Singular.svector)
  pv = Tuple{Int, elem_type(base_ring(F))}[]
  pos = Int[]
  values = []
  Rx = base_ring(F)
  R = base_ring(Rx)
  for (i, e, c) = s
    f = Base.findfirst(x->x==i, pos)
    if f === nothing
      push!(values, MPolyBuildCtx(base_ring(F).R))
      f = length(values)
      push!(pos, i)
    end
    push_term!(values[f], R(c), e)
  end
  pv = Tuple{Int, elem_type(Rx)}[(pos[i], base_ring(F)(finish(values[i]))) for i=1:length(pos)]
  return FreeModuleElem_dec(sparse_row(base_ring(F), pv), F)
end

isgraded(F::FreeModule_dec) = isgraded(F.R)
isfiltrated(F::FreeModule_dec) = isfiltrated(F.R)

abstract type ModuleFPHom_dec end
abstract type Map_dec{T1, T2} <: Map{T1, T2, Hecke.HeckeMap, ModuleFPHom_dec} end
  
mutable struct FreeModuleHom_dec{T1, T2} <: Map_dec{T1, T2} 
  header::MapHeader
  Hecke.@declare_other
  function FreeModuleHom_dec(F::FreeModule_dec{T}, G::S, a::Array{<:Any, 1}) where {T, S}
#    @assert isfiltrated(F) || all(ishomogenous, a) #neccessary and suffient according to Hans XXX
#same as non-homogenous elements are required, this too must not be enforced
    @assert all(x->parent(x) == G, a)
    @assert length(a) == ngens(F)
    #for filtrations, all is legal...
    r= new{typeof(F), typeof(G)}()
    function im_func(x::FreeModuleElem_dec)
      @assert parent(x) == F
      b = zero(G)
      for (i,v) = x.r
        b += v*a[i]
      end
      return b
    end
    function pr_func(x::FreeModuleElem_dec)
      @assert parent(x) == G
      c = coordinates(x, BiModArray(s))
      return FreeModuleElem_dec(c, F)
    end
    function pr_func(x::S)
      @assert parent(x) == G
      #assume S == SubQuoElem_dec which cannot be asserted here, the type if defined too late
      if isdefined(G, :quo)
        c = coordinates(x, BiModArray(vcat(s, gens(G.quo))))
      else
        c = coordinates(x, BiModArray(s))
      end
      return FreeModuleElem_dec(c[1:ngens(G)], G)
    end
    r.header = MapHeader{typeof(F), typeof(G)}(F, G, im_func, pr_func)
    return r
  end
end
(h::FreeModuleHom_dec)(a::FreeModuleElem_dec) = image(h, a)

hom(F::FreeModule_dec{T}, G, a) where {T} = FreeModuleHom_dec(F, G, a)

#TODO: should work with general Map_dec
#TODO: which should be renamed
function ishomogenous(h::T) where {T <: FreeModuleHom_dec}
  first = true
  local d::GrpAbFinGenElem
  for i = gens(domain(h))
    hi = h(i)
    iszero(hi) && continue

    if first
      d = degree(hi) - degree(i)
      first = false
    end
    d == degree(hi) - degree(i) || return false
  end
  return true
end

function degree(h::T) where {T <: FreeModuleHom_dec}
  first = true
  local d::GrpAbFinGenElem
  for i = gens(domain(h))
    hi = h(i)
    iszero(hi) && continue
    if first
      d = degree(hi) - degree(i)
      first = false
    end
    d == degree(hi) - degree(i) || error("hom is not homogenous")
  end
  if first
    error("hom is zero")
  end
  return d
end

function homogenous_components(h::T) where {T <: FreeModuleHom_dec}
  c = Dict{GrpAbFinGenElem, typeof(h)}()
  d = Dict{GrpAbFinGenElem, Array{Int, 1}}()
  F = domain(h)
  im = elem_type(codomain(h))[]
  for i = 1:ngens(F)
    hi = h(gen(F, i))
    push!(im, hi)
    iszero(hi) && continue
    x = degree(hi) - degree(gen(F, i))
    if haskey(d, x)
      push!(d[x], i)
    else
      d[x] = [i]
    end
  end
  for (k,v) = d
    c[k] = hom(F, codomain(h), [j in v ? im[j] : zero(codomain(h)) for j=1:ngens(F)])
  end
  return c
end

mutable struct SubQuo_dec{T} <: ModuleFP_dec{T}
  #meant to represent sub+ quo mod quo - as lazy as possible
  F::FreeModule_dec{T}
  sub::BiModArray
  quo::BiModArray
  sum::BiModArray
  std_sub::BiModArray
  std_quo::BiModArray
  AbstractAlgebra.@declare_other
  function SubQuo_dec(F::FreeModule_dec{R}, O::Array{<:FreeModuleElem_dec, 1}) where {R}
    r = new{R}()
    r.F = F
    r.sub = BiModArray(O)
    r.sum = r.sub
    return r
  end
  function SubQuo_dec(S::SubQuo_dec, O::Array{<:FreeModuleElem_dec{L}, 1}) where {L}
    r = new{parent_type(L)}()
    r.F = S.F
    r.sub = S.sub
    r.quo = BiModArray(O)
    r.sum = BiModArray(vcat(collect(r.sub), collect(r.quo)))
    return r
  end
  function SubQuo_dec(F::FreeModule_dec{R}, s::Singular.smodule) where {R}
    r = new{R}()
    r.F = F
    r.sub = BiModArray(F, s)
    r.sum = r.sub
    if s.isGB
      r.std_sub = r.sub
    end
    return r
  end
  function SubQuo_dec(F::FreeModule_dec{R}, s::Singular.smodule, t::Singular.smodule) where {R}
    r = new{R}()
    r.F = F
    r.sub = BiModArray(F, s)
    if s.isGB
      r.std_sub = r.sub
    end
    r.quo = BiModArray(F, t)
    if t.isGB
      r.std_quo = r.quo
    end
    r.sum = BiModArray(collect(r.sub), collect(r.quo))
    return r
  end
end

function show(io::IO, SQ::SubQuo_dec)
  @show_name(io, SQ)
  @show_special(io, SQ)

  if isdefined(SQ, :quo)
    println(io, "Subquotient of ", SQ.sub, " by ", SQ.quo)
  else
    println(io, "Subquotient by ", SQ.sub)
  end
end

@doc Markdown.doc"""
  A subquotien is (internally) given wia two submodules A and B of the same 
  FreeModule F. It represents $(A+B)/B$, so elements are given as elements
  in $A+B$
"""
struct SubQuoElem_dec{T}
 a::FreeModuleElem_dec{T}
 parent::SubQuo_dec
end

elem_type(::SubQuo_dec{T}) where {T} = SubQuoElem_dec{elem_type(T)}
parent_type(::SubQuoElem_dec{T}) where {T} = SubQuo_dec{parent_type(T)}
elem_type(::Type{SubQuo_dec{T}}) where {T} = SubQuoElem_dec{elem_type(T)}
parent_type(::Type{SubQuoElem_dec{T}}) where {T} = SubQuo_dec{parent_type(T)}

function sum_gb_assure(SQ::SubQuo_dec)
  singular_assure(SQ.sum)
  if SQ.sum.S.isGB
    return
  end
  SQ.sum = groebner_basis(SQ.sum)
end

function groebner_basis(F::BiModArray)
  singular_assure(F)
  if F.S.isGB
    return F
  end
  return BiModArray(F.F, Singular.std(F.S))
end

function show(io::IO, b::SubQuoElem_dec)
  print(io, b.a)
end

parent(b::SubQuoElem_dec) = b.parent

function (R::SubQuo_dec)(a::FreeModuleElem_dec; check::Bool = true)
  if check
    b = convert(R.sum.SF, a)
    sum_gb_assure(R)
    c = _reduce(b, R.sum.S)
    iszero(x) || error("not in the module")
  end
  return SubQuoElem_dec(a, R)
end

+(a::SubQuoElem_dec, b::SubQuoElem_dec) = SubQuoElem_dec(a.a+b.a, a.parent)
-(a::SubQuoElem_dec, b::SubQuoElem_dec) = SubQuoElem_dec(a.a-b.a, a.parent)
-(a::SubQuoElem_dec) = SubQuoElem_dec(-a.a, a.parent)
*(a::MPolyElem_dec, b::SubQuoElem_dec) = SubQuoElem_dec(a*b.a, b.parent)
*(a::MPolyElem, b::SubQuoElem_dec) = SubQuoElem_dec(a*b.a, b.parent)
*(a::Int, b::SubQuoElem_dec) = SubQuoElem_dec(a*b.a, b.parent)
==(a::SubQuoElem_dec, b::SubQuoElem_dec) = iszero(a-b)

function sub(F::FreeModule_dec, O::Array{<:FreeModuleElem_dec, 1})
  all(ishomogenous, O) || error("generators have to be homogenous")
  s = SubQuo_dec(F, O)
end

function sub(F::FreeModule_dec, O::Array{<:SubQuoElem_dec, 1})
  all(ishomogenous, O) || error("generators have to be homogenous")
  return SubQuo_dec(F, [x.a for x = O])
end

function sub(F::FreeModule_dec, s::SubQuo_dec)
  @assert !isdefined(s, :quo)
  return s
end

function quo(F::FreeModule_dec, O::Array{<:FreeModuleElem_dec, 1})
  S = SubQuo_dec(F, basis(F))
  return SubQuo_dec(S, O)
end

function quo(F::FreeModule_dec, O::Array{<:SubQuoElem_dec, 1})
  S = SubQuo_dec(F, basis(F))
  return SubQuo_dec(S, [x.a for x = O])
end

function quo(F::SubQuo_dec, O::Array{<:FreeModuleElem_dec, 1})
  all(ishomogenous, O) || error("generators have to be homogenous")
  if isdefined(F, :quo)
    t = BiModArray(F, O)
    t[Val{:S}, 1]
    F.quo[Val{:S}, 1]
    F.sub[Val{:S}, 1]
    s = t.S + F.quo.S
    return SubQuo_dec(F.F, F.sub.S, s)
  end
  return SubQuo_dec(F, O)
end

function quo(S::SubQuo_dec, O::Array{<:SubQuoElem_dec, 1})
  return SubQuo_dec(S, [x.a for x = O])
end

function quo(S::SubQuo_dec, T::SubQuo_dec)
  @assert !isdefined(T, :quo)
  return quo(S, gens(T))
end

function syzygie_module(F::BiModArray; sub = 0)
  F[Val(:S), 1]
  s = Singular.syz(F.S)
  if sub !== 0
    G = sub
  else
    G = FreeModule(base_ring(F.F), [degree(x) for x = F.O])
  end
  return SubQuo_dec(G, s)
end

function gens(F::SubQuo_dec)
  return map(x->SubQuoElem_dec(x, F), collect(F.sub))
end

function gen(F::SubQuo_dec, i::Int)
  return SubQuoElem_dec(F.sub[Val(:O), i], F)
end

ngens(F::SubQuo_dec) = length(F.sub)
base_ring(SQ::SubQuo_dec) = base_ring(SQ.F)

zero(SQ::SubQuo_dec) = SubQuoElem_dec(zero(SQ.F), SQ)

function Base.iterate(F::BiModArray, i::Int = 1)
  if i>length(F)
    return nothing
  else
    return F[i], i+1
  end
end
Base.eltype(::BiModArray{T}) where {T} = FreeModuleElem_dec{T} 

#??? A scalar product....
function *(a::FreeModuleElem_dec, b::Array{FreeModuleElem_dec, 1})
  @assert dim(parent(a)) == length(b)
  s = zero(parent(a))
  for (p,v) = a.r
    s += v*b[p]
  end
  return s
end

decoration(F::FreeModule_dec) = decoration(F.R)
decoration(R::MPolyRing_dec) = R.D
decoration(SQ::SubQuo_dec) = decoration(SQ.F)

function presentation(SQ::SubQuo_dec)
  #A+B/B is generated by A and B
  #the relations are A meet B? written wrt to A
  s = syzygie_module(SQ.sum)
  #TODO: wait for Hans to release Modulo(A, B) that does exactly this
  c = collect(s.sub)
  R = base_ring(SQ)
  F = FreeModule(R, [degree(x) for x = collect(SQ.sub)])
  q = elem_type(F)[]
  w = GrpAbFinGenElem[]

  for x = c
    b = sparse_row(R)
    e = zero(SQ.F)
    for (i,v) = x.r
      if i>ngens(F)
        break
      end
      e += v*gen(SQ.F, i)
      push!(b.pos, i)
      push!(b.values, v)
    end
    if iszero(e)
      continue
    end
    push!(q, FreeModuleElem_dec(b, F))
    push!(w, degree(q[end]))
  end
  #want R^a -> R^b -> SQ -> 0
  G = FreeModule(R, w)
  h_G_F = hom(G, F, q)
  @assert iszero(degree(h_G_F))
  h_F_SQ = hom(F, SQ, gens(SQ))
  @assert iszero(degree(h_F_SQ))
  Z = FreeModule(F.R, GrpAbFinGenElem[])
  Hecke.set_special(Z, :name => "Zero")
  h_SQ_Z = hom(SQ, Z, [zero(Z) for i=1:ngens(SQ)])
  return Hecke.ChainComplex(Oscar.ModuleFP_dec, Oscar.Map_dec[h_G_F, h_F_SQ, h_SQ_Z], check = false)
end

mutable struct SubQuoHom_dec{T1, T2} <: Map_dec{T1, T2}
  header::Hecke.MapHeader
  im::Array{<:Any, 1}
  function SubQuoHom_dec(D::SubQuo_dec, C::ModuleFP_dec, im::Array{<:Any, 1})
    first = true
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) == C, im)

    local deg::GrpAbFinGenElem
    b = gens(D)
    for i=1:length(im)
      if iszero(im[i])
        continue
      end
      if first
        deg = degree(b[i]) - degree(im[i])
        first = false
      else
        @assert deg == degree(b[i]) - degree(im[i])
      end
    end
    r = new{SubQuo_dec, typeof(C)}()
    r.header = Hecke.MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.im =  im
    return r
  end
end

function coordinates(a::FreeModuleElem_dec, SQ::SubQuo_dec)
  if iszero(a)
    return sparse_row(base_ring(parent(a)))
  end
  singular_assure(SQ.sum)
  b = BiModArray([a], SQ.sum.SF)
  singular_assure(b)
  s, r = Singular.lift(SQ.sum.S, b.S)
  if Singular.ngens(s) == 0 || iszero(s[1])
    error("elem not in module")
  end
  p = Int[]
  v = MPolyBuildCtx[]
  Sx = base_ring(SQ)
  Rx = Sx.R
  R = base_ring(Rx)
  for (i, e, c) = s[1]
    if i > ngens(SQ)
      break
    end
    if i in p
      f = findfirst(x->x==i, p)
    else
      push!(p, i)
      push!(v, MPolyBuildCtx(Rx))
      f = length(v)
    end
    push_term!(v[f], R(c), e)
  end
  pv = Tuple{Int, elem_type(Sx)}[(p[i], Sx(finish(v[i]))) for i=1:length(p)]
  return sparse_row(Sx, pv)
end

hom(D::SubQuo_dec, C::ModuleFP_dec, A::Array{<:Any, 1}) = SubQuoHom_dec(D, C, A)

function image(f::SubQuoHom_dec, a::SubQuoElem_dec)
  i = zero(codomain(f))
  D = domain(f)
  b = coordinates(a.a, D)
  for (p,v) = b
    i += v*f.im[p]
  end
  return i
end

function image(f::SubQuoHom_dec, a::FreeModuleElem_dec)
  i = zero(codomain(f))
  D = domain(f)
  b = coordinates(a, D)
  for (p,v) = b.r
    i += v*f.im[p]
  end
  return i
end
(f::SubQuoHom_dec)(a::FreeModuleElem_dec) = image(f, a)
(f::SubQuoHom_dec)(a::SubQuoElem_dec) = image(f, a)

function degree(h::SubQuoHom_dec)
  b = basis(domain(h).F)
  first = true
  for i = 1:length(b)
    if !iszero(h.im[i])
      return degree(h.im[i]) - degree(b[i])
    end
  end
end
#ishom, homcompo missing

function degree(a::FreeModuleElem_dec, C::SubQuo_dec)
  if !isdefined(C, :quo)
    return degree(a)
  end
  if !isdefined(C, :std_quo)
    C.quo[Val(:S), 1]
    C.std_quo = BiModArray(C.quo.F, Singular.std(C.quo.S))
  end
  x = _reduce(convert(C.quo.SF, a), C.std_quo.S)
  return degree(convert(C.F, x))
end

function ishomogenous(a::SubQuoElem_dec)
  C = parent(a)
  if !isdefined(C, :quo)
    return ishomogenous(a.a)
  end
  if !isdefined(C, :std_quo)
    singular_assure(C.quo)
    C.std_quo = BiModArray(C.quo.F, Singular.std(C.quo.S))
  end
  x = _reduce(convert(C.quo.SF, a.a), C.std_quo.S)
  return ishomogenous(convert(C.F, x))
end

function iszero(a::SubQuoElem_dec)
  C = parent(a)
  if !isdefined(C, :quo)
    return iszero(a.a)
  end
  if !isdefined(C, :std_quo)
    singular_assure(C.quo)
    C.std_quo = BiModArray(C.quo.F, Singular.std(C.quo.S))
  end
  x = _reduce(convert(C.quo.SF, a.a), C.std_quo.S)
  return iszero(x)
end

function degree(a::SubQuoElem_dec)
  return degree(a.a, parent(a))
end

function hom(F::FreeModule_dec, G::FreeModule_dec)
  @assert base_ring(F) == base_ring(G)
  GH = FreeModule(F.R, [y-x for x = F.d for y = G.d])
  GH.S = [Symbol("($i -> $j)") for i = F.S for j = G.S]
  Hecke.set_special(GH, :show => Hecke.show_hom, :hom => (F, G))

  #list is g1 - f1, g2-f1, g3-f1, ...
  X = Hecke.MapParent(F, G, "homomorphisms")
  n = ngens(F)
  m = ngens(G)
  R = base_ring(F)
  function im(x::FreeModuleElem_dec)
    return hom(F, G, [FreeModuleElem_dec(x.r[R, (i-1)*m+1:i*m], G) for i=1:n])
  end
  function pre(h::FreeModuleHom_dec)
    s = sparse_row(F.R)
    o = 0
    for i=1:n
      for (p,v) = h(gen(F, i)).r
        push!(s.pos, o+p)
        push!(s.values, v)
      end
      o += m
    end
    return FreeModuleElem_dec(s, GH)
  end
  return GH, Hecke.MapFromFunc(im, pre, GH, X)
end

function kernel(h::FreeModuleHom_dec)  #ONLY for free modules...
  G = domain(h)
  @assert isa(codomain(h), FreeModule_dec)
  #TODO allow sub-quo here as well
  b = BiModArray(map(h, basis(G)))
  k = syzygie_module(b, sub = G)
  @assert k.F == G
  c = collect(k.sub)
  return k, hom(k, parent(c[1]), c)
end

function image(h::FreeModuleHom_dec)
  s = sub(codomain(h), map(h, basis(domain(h))))
  return s, hom(s, codomain(h), map(h, basis(domain(h))))
end

function image(h::SubQuoHom_dec)
  s = sub(codomain(h), h.im)
  return s, hom(s, codomain(h), h.im)
end

function kernel(h::SubQuoHom_dec)
  error("not done yet")
end

function hom(M::SubQuo_dec, N::SubQuo_dec)
  p1 = presentation(M)
  p2 = presentation(N)
  k, mk = kernel(map(p2, 1))
  #Janko: have R^t1 -- g1 = map(p2, 0) -> R^t0 -> G
  #kernel g1: k -> R^t1
  #source: Janko's CA script: https://www.mathematik.uni-kl.de/~boehm/lehre/17_CA/ca.pdf
  D = decoration(M)
  F = FreeModule(base_ring(M), [iszero(x) ? D[0] : degree(x) for x = gens(k)])
  g2 = hom(F, codomain(mk), collect(k.sub)) #not clean - but maps not (yet) working
  #step 2
  H_s0_t0, mH_s0_t0 = hom(domain(map(p1, 2)), domain(map(p2, 2)))
  H_s1_t1, mH_s1_t1 = hom(domain(map(p1, 1)), domain(map(p2, 1)))
  D, pro, emb = direct_product(H_s0_t0, H_s1_t1, task = :both)

  H_s1_t0, mH_s1_t0 = hom(domain(map(p1, 1)), domain(map(p2, 2)))

  delta = hom(D, H_s1_t0, [preimage(mH_s1_t0, map(p1, 1)*mH_s0_t0(pro[1](g))-mH_s1_t1(pro[2](g))*map(p2, 1)) for g = gens(D)])

  H_s0_t1, mH_s0_t1 = hom(domain(map(p1, 2)), domain(map(p2, 1)))
  H_s1_t2, mH_s1_t2 = hom(domain(map(p1, 1)), F)

  E, pr = direct_product(H_s0_t1, H_s1_t2, task = :prod)

  rho = hom(E, D, [emb[1](preimage(mH_s0_t0, mH_s0_t1(pr[1](g))*map(p2, 1))) + 
                   emb[2](preimage(mH_s1_t1, map(p1, 1)*mH_s0_t1(pr[1](g)) - mH_s1_t2(pr[2](g))*g2)) for g = gens(E)])
  #need quo(image(delta), kern(rho))                 
  #return delta, rho
 
  H = quo(sub(D, kernel(delta)[1]), image(rho)[1])
  Hecke.set_special(H, :show => Hecke.show_hom, :hom => (M, N))

  #x in ker delta: mH_s0_t0(pro[1](x)) should be a hom from M to N
  function im(x::SubQuoElem_dec)
    @assert parent(x) == H
    return hom(M, N, [map(p2, 2)(mH_s0_t0(pro[1](x.a))(g)) for g = gens(domain(map(p1, 2)))])
  end

  function pr(f::SubQuoHom_dec)
    @assert domain(f) == M
    @assert codomain(f) == N
    Rs0 = domain(map(p1, 2))
    Rt0 = domain(map(p2, 2))
    g = hom(Rs0, Rt0, [FreeModuleElem_dec(preimage(map(p2, 2), f(map(p1, 2)(g)), N), Rt0) for g = gens(Rs0)])
    SubQuoElem_dec(emb[1](preimage(mH_s0_t0, g)), H)
  end
  return H, MapFromFunc(im, pr, H, Hecke.MapParent(M, N, "homomorphisms"))
end
#TODO
#  replace the +/- for the homs by proper constructors for homs and direct sums
#  relshp to store the maps elsewhere

function *(h::FreeModuleHom_dec, g::FreeModuleHom_dec) 
  @assert codomain(h) == domain(g)
  return hom(domain(h), codomain(g), [g(h(x)) for x = gens(domain(h))])
end
-(h::FreeModuleHom_dec, g::FreeModuleHom_dec) = hom(domain(h), codomain(h), [h(x) - g(x) for x = gens(domain(h))])
+(h::FreeModuleHom_dec, g::FreeModuleHom_dec) = hom(domain(h), codomain(h), [h(x) + g(x) for x = gens(domain(h))])

function direct_product(F::FreeModule_dec{T}...; task::Symbol = :sum) where {T}
  R = base_ring(F[1])
  G = FreeModule(R, vcat([f.d for f = F]...))
  G.S = []
  for i = 1:length(F)
    s = "("
    for j=1:i-1
      s *= "0, "
    end
    e = ""
    if i<length(F)
      e*=", "
    end
    for j=i+1:length(F)-1
      e *= "0, "
    end
    if i<length(F)
      e *= "0"
    end
    e*=")"

    for t = F[i].S
      push!(G.S, Symbol(s*string(t)*e))
    end
  end
  Hecke.set_special(G, :show => Hecke.show_direct_product, :direct_product => F)
  emb = []
  pro = []
  i = 0
  for f = F
    if task in [:sum, :both]
      push!(emb, hom(f, G, [gen(G, j+i) for j=1:ngens(f)]))
    end
    if task in [:prod, :both]
      push!(pro, hom(G, f, vcat(elem_type(f)[zero(f) for j=1:i], gens(f), elem_type(f)[zero(f) for j=i+ngens(f)+1:ngens(G)])))
    end
    i += ngens(f)
  end
  if task == :none
    return G
  elseif task == :sum
    return G, emb
  elseif task == :prod
    return G, pro
  elseif task == :both
    return G, pro, emb
  end
end

#############################
#TODO move to Hecke
function Base.getindex(r::Hecke.SRow, R::AbstractAlgebra.Ring, u::UnitRange)
  s = sparse_row(R)
  shift = 1-first(u)
  for (p,v) = r
    if p in u
      push!(s.pos, p+shift)
      push!(s.values, v)
    end
  end
  return s
end
##############################
#should be in Singular.jl
function Singular.intersection(a::Singular.smodule, b::Singular.smodule)
  c = base_ring(a)
  return Singular.Module(c, Singular.libSingular.id_Intersection(a.ptr, b.ptr, c.ptr))
end

function _reduce(a::Singular.smodule, b::Singular.smodule)
  @assert b.isGB
  p = Singular.libSingular.p_Reduce(a.ptr, b.ptr, base_ring(b).ptr)
  return Singular.Module(base_ring(b), p)
end

function _reduce(a::Singular.svector, b::Singular.smodule)
  @assert b.isGB
  p = _reduce(Singular.Module(base_ring(b), a), b)[1]
  return Singular.Module(base_ring(b), p)[1]
end


