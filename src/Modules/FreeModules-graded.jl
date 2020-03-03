#TODO make d and S a function optionally - to support HUGE degree
mutable struct FreeModule_dec
  d::Array{GrpAbFinGenElem, 1}
  R::MPolyRing_dec
  S::Array{Symbol, 1}
  Hecke.@declare_other
  function FreeModule_dec(a,b,c)
    r = new()
    r.d = a
    r.R = b
    r.S = c
    return r
  end
end

function FreeModule(R::MPolyRing_dec, n::Int; cached::Bool = false)
  return FreeModule_dec([R.D[0] for i=1:n], R, [Symbol("e[$i]") for i=1:n])
end

function FreeModule(R::MPolyRing_dec, d::Array{GrpAbFinGenElem, 1})
  return FreeModule_dec(d, R, [Symbol("e[$i]") for i=1:length(d)])
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

elem_type(::Type{FreeModule_dec}) = FreeModuleElem_dec
parent_type(::Type{FreeModuleElem_dec}) = FreeModule_dec

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
  bas = FreeModuleElem_dec[]
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

-(a::FreeModuleElem_dec) = FreeModuleElem_dec(-a.r, a.parent)
-(a::FreeModuleElem_dec, b::FreeModuleElem_dec) = FreeModuleElem_dec(a.r-b.r, a.parent)
+(a::FreeModuleElem_dec, b::FreeModuleElem_dec) = FreeModuleElem_dec(a.r+b.r, a.parent)
*(a::MPolyElem_dec, b::FreeModuleElem_dec) = FreeModuleElem_dec(a*b.r, b.parent)
==(a::FreeModuleElem_dec, b::FreeModuleElem_dec) = a.r == b.r
zero(F::FreeModule_dec) = FreeModuleElem_dec(sparse_row(F.R, Tuple{Int, elem_type(F.R)}[]), F)
parent(a::FreeModuleElem_dec) = a.parent

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

mutable struct BiModArray
  O::Array{FreeModuleElem_dec, 1}
  S::Singular.smodule
  F::FreeModule_dec
  SF::Singular.FreeMod
  function BiModArray(O::Array{<:FreeModuleElem_dec, 1})
    r = new()
    r.O = O
    r.SF = singular_module(parent(O[1]))
    r.F = parent(O[1])
    return r
  end
  function BiModArray(F::FreeModule_dec, s::Singular.smodule)
    r = new()
    r.F = F
    r.SF = parent(s[1])
    r.S = s
    r.O = Array{FreeModuleElem_dec, 1}(undef, Singular.ngens(s))
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

mutable struct SubQuo_dec 
  #meant to represent sub+ quo mod quo - as lazy as possible
  F::FreeModule_dec
  sub::BiModArray
  quo::BiModArray
  std_sub::BiModArray
  std_quo::BiModArray
  function SubQuo_dec(F::FreeModule_dec, O::Array{<:FreeModuleElem_dec, 1})
    r = new()
    r.F = F
    r.sub = BiModArray(O)
    return r
  end
  function SubQuo_dec(S::SubQuo_dec, O::Array{<:FreeModuleElem_dec, 1})
    r = new()
    r.F = S.F
    r.sub = S.sub
    r.quo = BiModArray(O)
    return r
  end
  function SubQuo_dec(F::FreeModule_dec, s::Singular.smodule)
    r = new()
    r.F = F
    r.sub = BiModArray(F, s)
    if s.isGB
      r.std_sub = r.sub
    end
    return r
  end
  function SubQuo_dec(F::FreeModule_dec, s::Singular.smodule, t::Singular.smodule)
    r = new()
    r.F = F
    r.sub = BiModArray(F, s)
    if s.isGB
      r.std_sub = r.sub
    end
    r.quo = BiModArray(F, t)
    if t.isGB
      r.std_quo = r.quo
    end
    return r
  end
end

function show(io::IO, SQ::SubQuo_dec)
  if isdefined(SQ, :quo)
    println(io, "Subquotient of ", SQ.sub, " by ", SQ.quo)
  else
    println(io, "Subquotient by ", SQ.sub)
  end
end

function sub(F::FreeModule_dec, O::Array{<:FreeModuleElem_dec, 1})
  all(ishomogenous, O) || error("generators have to be homogenous")
  return SubQuo_dec(F, O)
end

function quo(F::FreeModule_dec, O::Array{<:FreeModuleElem_dec, 1})
  S = SubQuo_dec(F, basis(F))
  return SubQuo_dec(S, O)
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

function syzygie_module(F::BiModArray)
  F[Val(:S), 1]
  s = Singular.syz(F.S)
  G = FreeModule(base_ring(F.F), [degree(x) for x = F.O])
  return SubQuo_dec(G, s)
end

function gens(F::SubQuo_dec)
  return F.sub
end

ngens(F::SubQuo_dec) = length(F.sub)

function iterate(F::BiModArray, i::Int = 1)
  if i>length(F)
    return nothing
  else
    return F[i], i+1
  end
end

function *(a::FreeModuleElem_dec, b::Array{FreeModuleElem_dec, 1})
  @assert dim(parent(a)) == length(b)
  s = zero(parent(a))
  for (p,v) = a.r
    s += v*b[p]
  end
  return s
end

iszero(a::FreeModuleElem_dec) = length(a.r) == 0

mutable struct SubQuoHom_dec{T1, T2} <: Map{T1, T2, Hecke.HeckeMap, SubQuoHom_dec}
  header::Hecke.MapHeader
  im::Array{FreeModuleElem_dec, 1}
  function SubQuoHom_dec(D::SubQuo_dec, C::SubQuo_dec, im::Array{<:FreeModuleElem_dec, 1})
    b = basis(D.F)
    first = true
    @assert length(im) == length(b)
    @assert parent(im[1]) == C.F #actually, membership should be tested..
    local deg::GrpAbFinGenElem
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
    r = new{SubQuo_dec, SubQuo_dec}()
    r.header = Hecke.MapHeader(D, C)
    r.im = im
    return r
  end
end

hom(D::SubQuo_dec, C::SubQuo_dec, A::Array{<:FreeModuleElem_dec, 1}) = SubQuoHom_dec(D, C, A)

function image(f::SubQuoHom_dec, a::FreeModuleElem_dec)
  i = zero(codomain(f).F)
  for (p,v) = a.r
    i += v*f.im[p]
  end
  return i
end
(f::SubQuoHom_dec)(a::FreeModuleElem_dec) = image(f, a)

function degree(h::SubQuoHom_dec)
  b = basis(domain(h).F)
  first = true
  for i = 1:length(b)
    if !iszero(h.im[i])
      return degree(h.im[i]) - degree(b[i])
    end
  end
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

function Singular.intersection(a::Singular.smodule, b::Singular.smodule)
  c = base_ring(a)
  return Singular.Module(c, Singular.libSingular.id_Intersection(a.ptr, b.ptr, c.ptr))
end

#=

  a = sum a_i e_i is homogenous <=> for all i: deg(a_i) + deg(e_i) is the same
  a -> b = sum a_i f_i  <=> for all i: deg(e_i) -deg(f_i) is constant  (or 0!!!)

    This only needs to work for the generators of the subquo
      and here, the degree is also dependent on the quo bit!
      What a mess

   Q: module elements have to be homogenous (if graded)?
       if not: hom can do different things on every homogenous component?
      how to check for filtration?
      cancellation might(?) render an element temporarily invalid(?)

hom(R^n, C/D) = (C/D)^n = oplus_1^n C/D = oplus(hom(R, C/D))      
  ni( f_1, ..., f_n), f_i in hom(R, C/D) (so f_i in C/D)

hom(A/B, C/D): A = <a_i | i>
  = oplus(C/D | i)/<all syz between a_i in A/B>
  ?????

=#
