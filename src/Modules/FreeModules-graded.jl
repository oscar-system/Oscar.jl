#TODO make d and S a function optionally - to support HUGE degree
export presentation, grading_group, ModuleOrdering

abstract type ModuleFP_dec{T} end

const Ring_dec = Union{MPolyRing_dec, MPolyQuo{<:Oscar.MPolyElem_dec}}
const RingElem_dec = Union{MPolyElem_dec, MPolyQuoElem{<:Oscar.MPolyElem_dec}}
#TODO: "fix" to allow QuoElem s as well...
# this requires
#  re-typeing of FreeModule
#  typing of BiModArray
# ... and all the rest.
# parametrization has to be by elem_type(coeff_ring) and not, like currently, the bottom coeff ring
# Also: qring is a Singular native. So it needs to be added to the ring creation

@attributes mutable struct FreeModule_dec{T} <: ModuleFP_dec{T}
  d::Vector{GrpAbFinGenElem}
  R::Ring_dec
  S::Vector{Symbol}

  function FreeModule_dec(a,b::Ring_dec,c)
    r = new{elem_type(b)}()
    r.d = a
    r.R = b
    r.S = c
    return r
  end
end

function FreeModule(R::Ring_dec, n::Int, name::String = "e"; cached::Bool = false) 
  return FreeModule_dec([grading_group(R)[0] for i=1:n], R, [Symbol("$name[$i]") for i=1:n])
end
free_module(R::Ring_dec, n::Int, name::String = "e"; cached::Bool = false) = FreeModule(R, n, name, cached = cached)

function FreeModule(R::Ring_dec, d::Vector{GrpAbFinGenElem}, name::String = "e"; cached::Bool = false)
  return FreeModule_dec(d, R, [Symbol("$name[$i]") for i=1:length(d)])
end
free_module(R::Ring_dec, d::Vector{GrpAbFinGenElem}, name::String = "e"; cached::Bool = false) = FreeModule(R, d, name, cached = cached)

#=XXX this cannot be as it is inherently ambiguous
  - FreeModule(R, n)
  - direct sum of rings, ie. a ring
  - set of n-th powers of R
thus the "category" needs to be set explicitly

^(R::Ring_dec, n::Int) = FreeModule(R, n)
=#

function AbstractAlgebra.extra_name(F::FreeModule_dec)
  t = get_attribute(F, :twist)
  if t !== nothing
    n = get_attribute(t[1], :name)
    if n !== nothing
      return "$n($(t[2]))"
    end
  end
  if length(Set(F.d)) == 1
    n = get_attribute(F.R, :name)
    if n !== nothing
      return "$n^$(ngens(F))($(-F.d[1]))"
    end
  end
  return nothing
end

function (F::FreeModule_dec)(a::GrpAbFinGenElem) 
  G = FreeModule(F.R, [x-a for x = F.d])
  set_attribute!(G, :twist => (F, a))
  return G
end

function (F::FreeModule_dec)()
  return FreeModuleElem_dec(sparse_row(base_ring(F)), F)
end

function show(io::IO, F::FreeModule_dec)
  @show_name(io, F)
  @show_special(io, F)

  print(io, "Free module of rank $(length(F.d)) over ")
  print(IOContext(io, :compact =>true), F.R)
  if is_graded(F.R)
    print(io, ", graded as ")
  else
    print(io, ", filtrated as ")
  end
  #TODO: on creation, in the resolution, sort the exponents
  md = MSet(F.d)
  first = true
  for (k,v) = md.dict
    if !first
      print(io, " + ")
    else
      first = false
    end
    print(IOContext(io, :compact => true), F.R, "^$v(", k, ")")
  end

#=
  i = 1
  while i < dim(F)
    d = F.d[i]
    j = 1
    while i+j <= dim(F) && d == F.d[i+j]
      j += 1
    end
    print(IOContext(io, :compact => true), F.R, "^$j")
    print(IOContext(io, :compact => true), "(", -d, ")")
    if i+j < dim(F)
      print(io, " + ")
    end
    i += j
  end
  =#
end

dim(F::FreeModule_dec)  = length(F.d)
ngens(F::FreeModule_dec) = dim(F)

struct FreeModuleElem_dec{T} <: AbstractAlgebra.ModuleElem{T}
  r::SRow{T}
  parent::FreeModule_dec{T}
end

elem_type(::Type{FreeModule_dec{T}}) where {T} = FreeModuleElem_dec{T}
parent_type(::Type{FreeModuleElem_dec{T}}) where {T} = FreeModule_dec{T}
elem_type(::FreeModule_dec{T}) where {T} = FreeModuleElem_dec{T}

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

function Base.getindex(F::FreeModule_dec, i::Int)
  i == 0 && return zero(F)
  return gen(F, i)
end

base_ring(F::FreeModule_dec) = F.R

#TODO: Parent - checks everywhere!!!

-(a::FreeModuleElem_dec) = FreeModuleElem_dec(-a.r, a.parent)
-(a::FreeModuleElem_dec, b::FreeModuleElem_dec) = FreeModuleElem_dec(a.r-b.r, a.parent)

function check_parent(a::FreeModuleElem_dec, b::FreeModuleElem_dec)
  if parent(a) !== parent(b)
    error("elements not compatible")
  end  
end

function +(a::FreeModuleElem_dec, b::FreeModuleElem_dec)
   check_parent(a, b)
   return FreeModuleElem_dec(a.r+b.r, a.parent)
end

*(a::MPolyElem_dec, b::FreeModuleElem_dec) = FreeModuleElem_dec(a*b.r, b.parent)
*(a::MPolyElem, b::FreeModuleElem_dec) = FreeModuleElem_dec(parent(b).R(a)*b.r, b.parent)
*(a::Int, b::FreeModuleElem_dec) = FreeModuleElem_dec(a*b.r, b.parent)
*(a::Integer, b::FreeModuleElem_dec) = FreeModuleElem_dec(b.parent.R(a)*b.r, b.parent)
*(a::fmpq, b::FreeModuleElem_dec) = FreeModuleElem_dec(b.parent.R(a)*b.r, b.parent)
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
    elseif is_graded(W)
      if ww != w
        error("elem not homogeneous")
      end
    else
      if W.lt(ww, w) 
        ww = w
      end
    end
  end
  return w
end

function homogeneous_components(a::FreeModuleElem_dec)
  res = Dict{GrpAbFinGenElem, FreeModuleElem_dec}()
  F = parent(a)
  for (p,v) = a.r
    c = homogeneous_components(v)
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

function homogeneous_component(a::FreeModuleElem_dec, g::GrpAbFinGenElem)
  F = parent(a)
  x = zero(F)
  for (p,v) = a.r
    x += homogeneous_component(v, g-F.d[p])*gen(F, p)
  end
  return x
end

function is_homogeneous(a::FreeModuleElem_dec)
  if iszero(a)
    return true
  end
  F = parent(a)
  first = true
  local d::GrpAbFinGenElem
  for (p,v) = a.r
    is_homogeneous(v) || return false
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
  O::Vector{FreeModuleElem_dec{T}}
  S::Singular.smodule
  F::FreeModule_dec
  SF::Singular.FreeMod

  function BiModArray(O::Vector{<:FreeModuleElem_dec{T}}) where {T}
    SF = singular_module(parent(O[1]))
    return BiModArray(O, SF)
  end

  function BiModArray(O::Vector{<:FreeModuleElem_dec{T}}, F::FreeModule_dec) where {T}
    SF = singular_module(F)
    return BiModArray(O, F, SF)
  end

  function BiModArray(O::Vector{<:FreeModuleElem_dec{T}}, SF::Singular.FreeMod) where {T}
    return BiModArray(O, parent(O[1]), SF)
  end

  function BiModArray(O::Vector{<:FreeModuleElem_dec{T}}, F::FreeModule_dec, SF::Singular.FreeMod) where {T}
    r = new{T}()
    r.O = O
    r.SF = SF
    r.F = F
    return r
  end

  function BiModArray(F::FreeModule_dec{S}, s::Singular.smodule) where {S}
    r = new{S}()
    r.F = F
    if Singular.ngens(s) == 0
      r.SF = Singular.FreeModule(base_ring(s), 0)
    else
      r.SF = parent(s[1])
    end
    r.S = s
    r.O = Vector{FreeModuleElem_dec{S}}(undef, Singular.ngens(s))
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
    F.S = Singular.smodule{elem_type(base_ring(F.SF))}(base_ring(F.SF), [convert(F.SF,x) for x = F.O]...)
  end
  return F.S[i]
end

function singular_assure(F::BiModArray)
  if length(F) == 0 && !isdefined(F, :S)
    F.S = Singular.smodule{elem_type(base_ring(F.SF))}(base_ring(F.SF), F.SF())
    return
  end
  F[Val(:S), 1]
end

getindex(F::BiModArray, i::Int) = getindex(F, Val(:O), i)

function singular_module(F::FreeModule_dec)
  Sx = singular_poly_ring(base_ring(F).R)
  return Singular.FreeModule(Sx, dim(F))
end

function convert(SF::Singular.FreeMod, m::FreeModuleElem_dec)
  g = Singular.gens(SF)
  e = SF()
  Sx = base_ring(SF)
  for (p,v) = m.r
    e += Sx(v.f)*g[p]
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

function (F::FreeModule_dec)(s::Singular.svector)
    convert(F, s)
end

is_graded(F::FreeModule_dec) = is_graded(F.R)
is_filtered(F::FreeModule_dec) = is_filtered(F.R)

abstract type ModuleFPHom_dec end
abstract type Map_dec{T1, T2} <: Map{T1, T2, Hecke.HeckeMap, ModuleFPHom_dec} end

mutable struct FreeModuleHom_dec{T1, T2} <: Map_dec{T1, T2} 
  header::MapHeader

  function FreeModuleHom_dec(F::FreeModule_dec{T}, G::S, a::Vector) where {T, S}
#    @assert is_filtered(F) || all(is_homogeneous, a) #necessary and sufficient according to Hans XXX
#same as non-homogeneous elements are required, this too must not be enforced
    @assert all(x->parent(x) == G, a)
    @assert length(a) == ngens(F)
    #for filtrations, all is legal...
    r = new{typeof(F), typeof(G)}()
    function im_func(x::FreeModuleElem_dec)
      b = zero(G)
      for (i,v) = x.r
        b += v*a[i]
      end
      return b
    end
    function pr_func(x::FreeModuleElem_dec)
      @assert parent(x) == G
      c = coordinates(x, sub(G, a))
      return FreeModuleElem_dec(c, F)
    end
    function pr_func(x)
      @assert parent(x) == G
      #assume S == SubQuoElem_dec which cannot be asserted here, the type if defined too late
      c = coordinates(x.a, sub(G, a))
      return FreeModuleElem_dec(c, F)
    end
    r.header = MapHeader{typeof(F), typeof(G)}(F, G, im_func, pr_func)
    return r
  end
end
(h::FreeModuleHom_dec)(a::FreeModuleElem_dec) = image(h, a)

hom(F::FreeModule_dec{T}, G, a) where {T} = FreeModuleHom_dec(F, G, a)

function identity_map(M::ModuleFP_dec)
  return hom(M, M, gens(M))
end

function is_homogeneous(h::T) where {T <: Map_dec}
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

function degree(h::T) where {T <: Map_dec}
  first = true
  local d::GrpAbFinGenElem
  D = domain(h)
  R = base_ring(D)
  for i = gens(domain(h))
    hi = h(i)
    iszero(hi) && continue
    if first
      d = degree(hi) - degree(i)
      first = false
    end
    dd = degree(hi) - degree(i)
    if is_filtered(R)
      if R.lt(d, dd)
        d = dd
      end
    else
      d == dd || error("hom is not homogeneous")
    end
  end
  if first
    error("hom is zero")
  end
  return d
end

function homogeneous_components(h::T) where {T <: Map_dec}
  c = Dict{GrpAbFinGenElem, typeof(h)}()
  d = Dict{GrpAbFinGenElem, Vector{Int}}()
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

@attributes mutable struct SubQuo_dec{T} <: ModuleFP_dec{T}
  #meant to represent sub+ quo mod quo - as lazy as possible
  F::FreeModule_dec{T}
  sub::BiModArray
  quo::BiModArray
  sum::BiModArray
  std_sub::BiModArray
  std_quo::BiModArray

  function SubQuo_dec(F::FreeModule_dec{R}, O::Vector{<:FreeModuleElem_dec}) where {R}
    r = new{R}()
    r.F = F
    r.sub = BiModArray(O, F, singular_module(F))
    r.sum = r.sub
    return r
  end
  function SubQuo_dec(S::SubQuo_dec, O::Vector{<:FreeModuleElem_dec{L}}) where {L}
    r = new{L}()
    r.F = S.F
    r.sub = S.sub
    r.quo = BiModArray(O, S.F, S.sub.SF)
    r.sum = BiModArray(vcat(collect(r.sub), collect(r.quo)), S.F, S.sub.SF)
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
    r.sum = BiModArray(vcat(collect(r.sub), collect(r.quo)))
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
  A subquotient is (internally) given wia two submodules A and B of the same 
  FreeModule F. It represents $(A+B)/B$, so elements are given as elements
  in $A+B$
"""
struct SubQuoElem_dec{T}
 a::FreeModuleElem_dec{T}
 parent::SubQuo_dec
end

elem_type(::SubQuo_dec{T}) where {T} = SubQuoElem_dec{T}
parent_type(::SubQuoElem_dec{T}) where {T} = SubQuo_dec{T}
elem_type(::Type{SubQuo_dec{T}}) where {T} = SubQuoElem_dec{T}
parent_type(::Type{SubQuoElem_dec{T}}) where {T} = SubQuo_dec{T}

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
    iszero(c) || error("not in the module")
  end
  return SubQuoElem_dec(a, R)
end

function (R::SubQuo_dec)(a::SubQuoElem_dec)
  if parent(a) == R
    return a
  end
  error("illegal coercion")
end

+(a::SubQuoElem_dec, b::SubQuoElem_dec) = SubQuoElem_dec(a.a+b.a, a.parent)
-(a::SubQuoElem_dec, b::SubQuoElem_dec) = SubQuoElem_dec(a.a-b.a, a.parent)
-(a::SubQuoElem_dec) = SubQuoElem_dec(-a.a, a.parent)
*(a::MPolyElem_dec, b::SubQuoElem_dec) = SubQuoElem_dec(a*b.a, b.parent)
*(a::MPolyElem, b::SubQuoElem_dec) = SubQuoElem_dec(a*b.a, b.parent)
*(a::Int, b::SubQuoElem_dec) = SubQuoElem_dec(a*b.a, b.parent)
*(a::Integer, b::SubQuoElem_dec) = SubQuoElem_dec(a*b.a, b.parent)
*(a::fmpq, b::SubQuoElem_dec) = SubQuoElem_dec(a*b.a, b.parent)
==(a::SubQuoElem_dec, b::SubQuoElem_dec) = iszero(a-b)

function sub(F::FreeModule_dec, O::Vector{<:FreeModuleElem_dec})
  all(is_homogeneous, O) || error("generators have to be homogeneous")
  s = SubQuo_dec(F, O)
end

function sub(F::FreeModule_dec, O::Vector{<:SubQuoElem_dec})
  all(is_homogeneous, O) || error("generators have to be homogeneous")
  return SubQuo_dec(F, [x.a for x = O])
end

function sub(F::FreeModule_dec, s::SubQuo_dec)
  @assert !isdefined(s, :quo)
  return s
end

function sub(S::SubQuo_dec, O::Vector{<:SubQuoElem_dec})
  t = sub(S.F, O)
  if isdefined(S, :quo)
    return quo(t, collect(S.quo))
  else
    return t
  end
end

function quo(F::FreeModule_dec, O::Vector{<:FreeModuleElem_dec})
  S = SubQuo_dec(F, basis(F))
  return SubQuo_dec(S, O)
end

function quo(F::FreeModule_dec, O::Vector{<:SubQuoElem_dec})
  S = SubQuo_dec(F, basis(F))
  return SubQuo_dec(S, [x.a for x = O])
end

function quo(F::SubQuo_dec, O::Vector{<:FreeModuleElem_dec})
  all(is_homogeneous, O) || error("generators have to be homogeneous")
  @assert parent(O[1]) == F.F
  if isdefined(F, :quo)
    F.sub[Val(:S), 1]
    [F.quo[Val(:O), i] for i = 1:length(F.quo.O)]
    s = Singular.smodule{elem_type(base_ring(F.quo.SF))}(base_ring(F.quo.SF), [convert(F.quo.SF, x) for x = [O; F.quo.O]]...)
    return SubQuo_dec(F.F, F.sub.S, s)
  end
  return SubQuo_dec(F, O)
end

function quo(S::SubQuo_dec, O::Vector{<:SubQuoElem_dec})
  return SubQuo_dec(S, [x.a for x = O])
end

function quo(S::SubQuo_dec, T::SubQuo_dec)
#  @assert !isdefined(T, :quo)
  return SubQuo_dec(S, T.sum.O)
end

function quo(F::FreeModule_dec, T::SubQuo_dec)
  @assert !isdefined(T, :quo)
  return quo(F, gens(T))
end

function syzygie_module(F::BiModArray; sub = 0)
  F[Val(:S), 1] #to force the existence of F.S
  s = Singular.syz(F.S)
  if sub !== 0
    G = sub
  else
  	[F[Val(:O), i] for i = 1:length(F.O)]
    z = grading_group(base_ring(F.F))[0]
    G = FreeModule(base_ring(F.F), [iszero(x) ? z : degree(x) for x = F.O])
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

function Base.iszero(F::SubQuo_dec)
  return all(iszero, gens(F))
end

function Base.getindex(F::SubQuo_dec, i::Int)
  i == 0 && return zero(F)
  return gen(F, i)
end

function Base.iterate(F::BiModArray, i::Int = 1)
  if i>length(F)
    return nothing
  else
    return F[i], i+1
  end
end
Base.eltype(::BiModArray{T}) where {T} = FreeModuleElem_dec{T} 

#??? A scalar product....
function *(a::FreeModuleElem_dec, b::Vector{FreeModuleElem_dec})
  @assert dim(parent(a)) == length(b)
  s = zero(parent(a))
  for (p,v) = a.r
    s += v*b[p]
  end
  return s
end

grading_group(F::FreeModule_dec) = grading_group(F.R)
grading_group(SQ::SubQuo_dec) = grading_group(SQ.F)

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
      if i>ngens(SQ)
        break
      end
      e += v*gen(SQ, i).a
      push!(b.pos, i)
      push!(b.values, v)
    end
    if length(b) == 0
      continue
    end
    push!(q, FreeModuleElem_dec(b, F))
    push!(w, degree(q[end]))
  end
  #want R^a -> R^b -> SQ -> 0
  #TODO sort grading_group and fix maps, same grading_group should be bundled (to match pretty printing)
  G = FreeModule(R, w)
  h_G_F = hom(G, F, q)
  @assert iszero(h_G_F) || iszero(degree(h_G_F))
  h_F_SQ = hom(F, SQ, gens(SQ))
  @assert iszero(h_F_SQ) || iszero(degree(h_F_SQ))
  Z = FreeModule(F.R, GrpAbFinGenElem[])
  set_attribute!(Z, :name => "Zero")
  h_SQ_Z = hom(SQ, Z, [zero(Z) for i=1:ngens(SQ)])
  return Hecke.ChainComplex(Oscar.ModuleFP_dec, Oscar.Map_dec[h_G_F, h_F_SQ, h_SQ_Z], check = false)
end

function presentation(F::FreeModule_dec)
  Z = FreeModule(F.R, GrpAbFinGenElem[])
  set_attribute!(Z, :name => "Zero")
  return Hecke.ChainComplex(ModuleFP_dec, Map_dec[hom(Z, F, FreeModuleElem_dec[]), hom(F, F, gens(F)), hom(F, Z, [zero(Z) for i=1:ngens(F)])], check = false)
end

mutable struct SubQuoHom_dec{T1, T2} <: Map_dec{T1, T2}
  header::Hecke.MapHeader
  im::Vector
  function SubQuoHom_dec(D::SubQuo_dec, C::ModuleFP_dec, im::Vector)
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
    r.header.preimage = x->preimage(r, x)
    r.im =  im
    return r
  end
end

function hom_tensor(G::ModuleFP_dec, H::ModuleFP_dec, A::Vector{ <: Map_dec})
  tG = get_attribute(G, :tensor_product)
  tG === nothing && error("both modules must be tensor products")
  tH = get_attribute(H, :tensor_product)
  tH === nothing && error("both modules must be tensor products")
  @assert length(tG) == length(tH) == length(A)
  @assert all(i-> domain(A[i]) == tG[i] && codomain(A[i]) == tH[i], 1:length(A))
  #gens of G are G[i][j] tensor G[h][l] for i != h and all j, l
  #such a pure tensor is mapped to A[i](G[i][j]) tensor A[h](G[j][l])
  #thus need the pure map - and re-create the careful ordering of the generators as in the 
  # constructor
  #store the maps? and possibly more data, like the ordeing
  error("not done yet")
  return hom(G, H)
end

function hom_product(G::ModuleFP_dec, H::ModuleFP_dec, A::Matrix{ <: Map_dec})
  tG = get_attribute(G, :tensor_product)
  tG === nothing && error("both modules must be direct products")
  tH = get_attribute(H, :tensor_product)
  tH === nothing && error("both modules must be direct products")
  @assert length(tG) == size(A, 1) && length(tH) == size(A, 2)
  @assert all((i,j)-> domain(A[i,j]) == tG[i] && codomain(A[i,j]) == tH[j], Base.Iterators.ProductIterator((1:size(A, 1), 1:size(A, 2))))
  #need the canonical maps..., maybe store them as well?
  error("not done yet")
end
# hom(prod -> X), hom(x -> prod)
# if too much time: improve the hom(A, B) in case of A and/or B are products - or maybe not...
# tensor and hom functors for chain complex
# dual: ambig: hom(M, R) or hom(M, Q(R))?

function coordinates(a::FreeModuleElem_dec, SQ::SubQuo_dec)
  if iszero(a)
    return sparse_row(base_ring(parent(a)))
  end
  [SQ.sub[Val(:O), i] for i = 1:length(SQ.sub.O)]
  if isdefined(SQ, :quo)
    [SQ.quo[Val(:O), i] for i = 1:length(SQ.quo.O)]
    generators = vcat(SQ.sub.O, SQ.quo.O)
  else
    generators = SQ.sub.O
  end
  S = Singular.smodule{elem_type(base_ring(SQ.sub.SF))}(base_ring(SQ.sub.SF), [convert(SQ.sub.SF, x) for x = generators]...)
  b = BiModArray([a], SQ.sum.SF)
  singular_assure(b)
  s, r = Singular.lift(S, b.S)
  if Singular.ngens(s) == 0 || iszero(s[1])
    error("elem not in module")
  end
  Sx = base_ring(SQ)
  Rx = Sx.R
  R = base_ring(Rx)
  return sparse_row(Sx, s[1], 1:ngens(SQ))
end


hom(D::SubQuo_dec, C::ModuleFP_dec, A::Vector) = SubQuoHom_dec(D, C, A)

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
  for (p,v) = b
    i += v*f.im[p]
  end
  return i
end

function preimage(f::SubQuoHom_dec, a::FreeModuleElem_dec)
  @assert parent(a) == codomain(f)
  D = domain(f)
  i = zero(D)
  b = coordinates(a, image(f)[1])
  for (p,v) = b
    i += v*gen(D, p)
  end
  return i
end

(f::SubQuoHom_dec)(a::FreeModuleElem_dec) = image(f, a)
(f::SubQuoHom_dec)(a::SubQuoElem_dec) = image(f, a)

function degree(h::SubQuoHom_dec)
  b = gens(domain(h))
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

function is_homogeneous(a::SubQuoElem_dec)
  C = parent(a)
  if !isdefined(C, :quo)
    return is_homogeneous(a.a)
  end
  if !isdefined(C, :std_quo)
    singular_assure(C.quo)
    C.std_quo = BiModArray(C.quo.F, Singular.std(C.quo.S))
  end
  x = _reduce(convert(C.quo.SF, a.a), C.std_quo.S)
  return is_homogeneous(convert(C.F, x))
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
  C = parent(a)
  if !isdefined(C, :quo)
    return degree(a.a)
  end
  if !isdefined(C, :std_quo)
    singular_assure(C.quo)
    C.std_quo = BiModArray(C.quo.F, Singular.std(C.quo.S))
  end
  x = _reduce(convert(C.quo.SF, a.a), C.std_quo.S)
  return degree(convert(parent(a.a), x))

end

function hom(F::FreeModule_dec, G::FreeModule_dec)
  @assert base_ring(F) == base_ring(G)
  GH = FreeModule(F.R, [y-x for x = F.d for y = G.d])
  GH.S = [Symbol("($i -> $j)") for i = F.S for j = G.S]
  set_attribute!(GH, :show => Hecke.show_hom, :hom => (F, G))

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
  if ngens(G) == 0
    s = sub(G, gens(G))
    return s, hom(s, G, gens(G))
  end
  g = map(h, basis(G))
  if isa(codomain(h), SubQuo_dec)
    g = [x.a for x = g]
    if isdefined(codomain(h), :quo)
      append!(g, collect(codomain(h).quo))
    end
  end
  #TODO allow sub-quo here as well
  b = BiModArray(g)
  k = syzygie_module(b)
  if isa(codomain(h), SubQuo_dec)
    s = collect(k.sub)
    k = sub(G, [FreeModuleElem_dec(x.r[1:dim(G)], G) for x = s])
  else
    #the syzygie_module creates a new free module to work in
    k = sub(G, [FreeModuleElem_dec(x.r, G) for x = collect(k.sub)])
  end
  @assert k.F == G
  c = collect(k.sub)
  return k, hom(k, parent(c[1]), c)
end

function image(h::FreeModuleHom_dec)
  si = [x for x = map(h, basis(domain(h))) if !iszero(x)]
  s = sub(codomain(h), si)
  return s, hom(s, codomain(h), si)
end

function image(h::SubQuoHom_dec)
  s = sub(codomain(h), h.im)
  return s, hom(s, codomain(h), h.im)
end

function kernel(h::SubQuoHom_dec)
  D = domain(h)
  R = base_ring(D)
  F = FreeModule(R, ngens(D))
  hh = hom(F, codomain(h), map(h, gens(D)))
  k = kernel(hh)
  @assert domain(k[2]) == k[1]
  @assert codomain(k[2]) == F
  hh = hom(F, D, gens(D))
  im = map(x->hh(k[2](x)), gens(k[1]))
  k = sub(D, im)
  return k, hom(k, D, im)
end

function free_resolution(F::FreeModule_dec)
  return presentation(F)
end

function free_resolution(S::SubQuo_dec, limit::Int = -1)
  p = presentation(S)
  mp = [map(p, j) for j=1:length(p)]
  D = grading_group(S)
  while true
    k, mk = kernel(mp[1])
    nz = findall(x->!iszero(x), gens(k))
    if length(nz) == 0 
      Z = FreeModule(base_ring(S), GrpAbFinGenElem[])
      set_attribute!(Z, :name => "Zero")
      h = hom(Z, domain(mp[1]), FreeModuleElem_dec[])
      insert!(mp, 1, h)
      break
    elseif limit != -1 && length(mp) > limit
      break
    end
    F = FreeModule(base_ring(S), [iszero(x) ? D[0] : degree(x) for x = gens(k)[nz]])
    g = hom(F, codomain(mk), collect(k.sub)[nz])
    insert!(mp, 1, g)
  end
  return Hecke.ChainComplex(ModuleFP_dec, mp, check = false, direction = :right)
end

function iszero(f::Map_dec)
  return all(iszero, map(f, gens(domain(f))))
end

function hom(M::ModuleFP_dec, N::ModuleFP_dec)
  p1 = presentation(M)
  p2 = presentation(N)
  k, mk = kernel(map(p2, 1))
  #Janko: have R^t1 -- g1 = map(p2, 0) -> R^t0 -> G
  #kernel g1: k -> R^t1
  #source: Janko's CA script: https://www.mathematik.uni-kl.de/~boehm/lehre/17_CA/ca.pdf
  D = grading_group(M)
  F = FreeModule(base_ring(M), GrpAbFinGenElem[iszero(x) ? D[0] : degree(x) for x = gens(k)])
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
  #need quo(kern(delta), image(rho))                 
 
  kDelta = kernel(delta)

  psi = kDelta[2]*pro[1]
  psi = hom(kDelta[1], H_s0_t0, [psi(g) for g = gens(kDelta[1])])

  H = quo(sub(D, kDelta[1]), image(rho)[1])
  set_attribute!(H, :show => Hecke.show_hom, :hom => (M, N))

  #x in ker delta: mH_s0_t0(pro[1](x)) should be a hom from M to N
  function im(x::SubQuoElem_dec)
    @assert parent(x) == H
    return hom(M, N, [map(p2, 2)(mH_s0_t0(pro[1](x.a))(preimage(map(p1, 2), g))) for g = gens(M)])
  end

  function pr(f::SubQuoHom_dec)
    @assert domain(f) == M
    @assert codomain(f) == N
    Rs0 = domain(map(p1, 2))
    Rt0 = domain(map(p2, 2))
    g = hom(Rs0, Rt0, [preimage(map(p2, 2), f(map(p1, 2)(g))) for g = gens(Rs0)])

    return H(preimage(psi, (preimage(mH_s0_t0, g))).a)
    return SubQuoElem_dec(emb[1](preimage(mH_s0_t0, g)), H)
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

##################################################
# direct product
##################################################
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
  set_attribute!(G, :show => Hecke.show_direct_product, :direct_product => F)
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

function direct_product(G::ModuleFP_dec...; task::Symbol = :none)
  F, pro, mF = direct_product([free_module(x) for x = G]..., task = :both)
  s = sub(F, vcat([[mF[i](y) for y = gens(G[i], free_module(G[i]))] for i=1:length(G)]...))
  q = vcat([[mF[i](y) for y = rels(G[i])] for i=1:length(G)]...)
  if length(q) != 0
    s = quo(s, q)
  end
  if task == :none
    return s
  elseif task == :prod
    return s, pro
  elseif task == :sum
    return s, mF
  else
    return s, pro, mF
  end
end
⊕(M::ModuleFP_dec...) = direct_product(M..., task = :none)


function Hecke.canonical_injection(G::ModuleFP_dec, i::Int)
  H = get_attribute(G, :direct_product)
  if H === nothing
    error("module not a direct product")
  end
  0<i<= length(H) || error("index out of bound")
  j = i == 1 ? 0 : sum(ngens(H[l]) for l=1:i-1) -1
  return hom(H[i], G, [G[l+j] for l = 1:ngens(H[i])])
end

function Hecke.canonical_projection(G::ModuleFP_dec, i::Int)
  H = get_attribute(G, :direct_product)
  if H === nothing
    error("module not a direct product")
  end
  0<i<= length(H) || error("index out of bound")
  j = i == 1 ? 0 : sum(ngens(H[l]) for l=1:i-1) 
  return hom(G, H[i], vcat([zero(H[i]) for l=1:j], gens(H[i]), [zero(H[i]) for l=1+j+ngens(H[i]):ngens(G)]))
end
    
##################################################
# Tensor
##################################################

function tensor_product(G::FreeModule_dec...; task::Symbol = :none)
  d = G[1].d
  s = G[1].S
  t = [[x] for x = 1:ngens(G[1])]
  for H = G[2:end]
    d = [x + y for x = d  for y = H.d]
    s = [Symbol("$x \\otimes $y") for x = s  for y = H.S]
    t = [push!(deepcopy(x), y) for x = t  for y = 1:ngens(H)]
  end

  F = FreeModule(G[1].R, d)
  F.S = s
  set_attribute!(F, :show => Hecke.show_tensor_product, :tensor_product => G)
  if task == :none
    return F
  end

  function pure(g::FreeModuleElem_dec...)
    @assert length(g) == length(G)
    @assert all(i-> parent(g[i]) == G[i], 1:length(G))
    z = [[x] for x = g[1].r.pos]
    zz = g[1].r.values
    for h = g[2:end]
      zzz = Vector{Int}[]
      zzzz = elem_type(F.R)[]
      for i = 1:length(z)
        for (p, v) = h.r
          push!(zzz, push!(deepcopy(z[i]), p))
          push!(zzzz, zz[i]*v)
        end
      end
      z = zzz
      zz = zzzz
    end
    return FreeModuleElem_dec(sparse_row(F.R, [findfirst(x->x == y, t) for y = z], zz), F)
  end
  function pure(T::Tuple)
    return pure(T...)
  end
  function inv_pure(e::FreeModuleElem_dec)
    if length(e.r.pos) == 0
      return Tuple(zero(g) for g = G)
    end
    @assert length(e.r.pos) == 1
    @assert isone(e.r.values[1])
    return Tuple(gen(G[i], t[e.r.pos[1]][i]) for i = 1:length(G))
  end

  return F, MapFromFunc(pure, inv_pure, Hecke.TupleParent(Tuple([g[0] for g = G])), F)
end

⊗(G::ModuleFP_dec...) = tensor_product(G..., task = :none)

function free_module(F::FreeModule_dec)
  return F
end

function free_module(F::SubQuo_dec)
  return F.F
end

gens(F::FreeModule_dec, G::FreeModule_dec) = gens(F)
gens(F::SubQuo_dec, G::FreeModule_dec) = [x.a for x = gens(F)]
rels(F::FreeModule_dec) = elem_type(F)[]
rels(F::SubQuo_dec) = isdefined(F, :quo) ? collect(F.quo) : elem_type(F.F)[]

@doc Markdown.doc"""
    tensor_product(G::ModuleFP_dec...; task::Symbol = :map) -> SubQuo_dec, Map

Given modules $G_i$ compute the tensor product $G_1\otimes \cdots \otimes G_n$.
If `task` is set to ":map", a map $\phi$ is returned that
maps tuples in $G_1 \times \cdots \times G_n$ to pure tensors
$g_1 \otimes \cdots \otimes g_n$. The map admits a preimage as well.
"""
function tensor_product(G::ModuleFP_dec...; task::Symbol = :none)
  F, mF = tensor_product([free_module(x) for x = G]..., task = :map)
  s = sub(F, vec([mF(x) for x = Base.Iterators.ProductIterator(Tuple(gens(x, free_module(x)) for x = G))]))
  q = vcat([vec([mF(x) for x = Base.Iterators.ProductIterator(Tuple(i == j ? rels(G[i]) : gens(free_module(G[i])) for i=1:length(G)))]) for j=1:length(G)]...) 
  if length(q) != 0
    s = quo(s, q)
  end
  if task == :none
    return s
  else
    return s, mF
  end
end

#TODO
#  (hom lift) => hom and tensor functor
#  filtrations
#  more constructors
#################################################
#
#################################################
function homogeneous_component(F::T, d::GrpAbFinGenElem) where {T <: Union{FreeModule_dec, SubQuo_dec, MPolyIdeal{<:MPolyElem_dec}}}

  #TODO: lazy: ie. no enumeration of points
  #      aparently it is possible to get the number of points faster than the points
  W = base_ring(F)
  D = grading_group(W)
  #have gens for W that can be combined
  #              F that can only be used
  #F ni f = sum c_i,j F[i]*w[j]
  #want deg(F[i]) + deg(w[j]) = d
  all = []
  for g = gens(F)
    if iszero(g)
      continue
    end
    Md, mMd = homogeneous_component(W, d - degree(g))
    #TODO careful <0> is 0-dim but non empty...
    if dim(Md) > 0
      push!(all, (g, mMd))
    end
  end

  B = elem_type(F)[]
  for (g, mMd) = all
    for x = gens(domain(mMd))
      y = mMd(x) * g
      iszero(y) && continue
      @assert !iszero(y)
      push!(B, y)
    end
  end

#TODO: vector_space(QQ, module elems)
  if isa(F, MPolyIdeal)
    X, h = vector_space(base_ring(base_ring(F)), B)
    set_attribute!(X, :show => show_homo_comp, :data => (F, d))
    return X, h
  end

  function im(f)
    sum(f[i]*B[i] for i=1:dim(X))
  end

  deg = length(B)
  X = FreeModule(base_ring(W), deg)
  set_attribute!(X, :show => show_homo_comp, :data => (F, d))

  function pr(g::S) where {S <: Union{FreeModuleElem_dec, SubQuoElem_dec}}
    #TODO: add X() and some sane setting of coeffs in FreeElems
    @assert elem_type(F) == typeof(g)
    @assert parent(g) == F
    z = zero(X)
    for (p,v) = g.r
      i = findfirst(x->gen(F, p) == x[1], all)
      j = preimage(all[i][2], v)
      o = i == 1 ?  0 : sum(x->dim(domain(x[2])), all[1:i-1]) 
      for k = 1:dim(parent(j))
        z.v[1,o+k] = j[k]
      end
    end
    return z
  end
  return X, Hecke.MapFromFunc(im, pr, X, F)
end
