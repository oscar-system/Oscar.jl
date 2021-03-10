#TODO make d and S a function optionally - to support HUGE degree
export presentation

abstract type ModuleFP{T} end

const CRing = Union{MPolyRing, MPolyQuo{<:Oscar.MPolyElem}, MPolyRing_dec, MPolyQuo{<:Oscar.MPolyElem_dec}}
const CRingElem = Union{MPolyElem, MPolyQuoElem{<:Oscar.MPolyElem}, MPolyElem_dec, MPolyQuoElem{<:Oscar.MPolyElem_dec}}
#TODO: "fix" to allow QuoElem s as well...
# this requires
#  re-typeing of FreeModule
#  typing of BiModArray
# ... and all the rest.
# parametrization has to be by elem_type(coeff_ring) and not, like currently, the bottom coeff ring
# Also: qring is a Singular native. So it needs to be added to the ring creation

mutable struct FreeMod{T} <: ModuleFP{T}
  R::CRing
  n::Int
  S::Array{Symbol, 1}
  AbstractAlgebra.@declare_other

  function FreeMod(n,b::CRing,c)
    r = new{elem_type(b)}()
    r.n = n
    r.R = b
    r.S = c
    return r
  end
end

function FreeMod(R::CRing, n::Int, name::String = "e"; cached::Bool = false) 
  return FreeMod(n, R, [Symbol("$name[$i]") for i=1:n])
end
free_module(R::CRing, n::Int, name::String = "e"; cached::Bool = false) = FreeMod(R, n, name, cached = cached)

#=XXX this cannot be as it is inherently ambigous
  - FreeModule(R, n)
  - direct sum of rings, ie. a ring
  - set of n-th powers of R
thus the "category" needs to be set explicitly

^(R::Ring_dec, n::Int) = FreeModule(R, n)
=#

function AbstractAlgebra.extra_name(F::FreeMod)
  if F.n == 1
    n = Hecke.get_special(F.R, :name)
    if n !== nothing
      return "$n^$(ngens(F))($(-F.d[1]))"
    end
  end
  return nothing
end

function (F::FreeMod)()
  return FreeModuleElem(sparse_row(base_ring(F)), F)
end

function show(io::IO, F::FreeMod)
  @show_name(io, F)
  @show_special(io, F)

  print(io, "Free module of rank $(F.n) over ")
  print(IOContext(io, :compact =>true), F.R)
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

dim(F::FreeMod) = F.n
ngens(F::FreeMod) = dim(F)

struct FreeModuleElem{T}
  r::SRow{T}
  parent::FreeMod{T}
end

elem_type(::Type{FreeMod{T}}) where {T} = FreeModuleElem{T}
parent_type(::Type{FreeModuleElem{T}}) where {T} = FreeMod{T}
elem_type(::FreeMod{T}) where {T} = FreeModuleElem{T}

function show(io::IO, e::FreeModuleElem)
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

function basis(F::FreeMod)
  bas = elem_type(F)[]
  for i=1:dim(F)
    s = Hecke.sparse_row(F.R, [(i, F.R(1))])
    push!(bas, FreeModuleElem(s, F))
  end
  return bas
end
gens(F::FreeMod) = basis(F)

function gen(F::FreeMod, i::Int)
  @assert 0 < i <= ngens(F)
  s = Hecke.sparse_row(F.R, [(i, F.R(1))])
  return FreeModuleElem(s, F)
end

function Base.getindex(F::FreeMod, i::Int)
  i == 0 && return zero(F)
  return gen(F, i)
end

base_ring(F::FreeMod) = F.R

#TODO: Parent - checks everywhere!!!

-(a::FreeModuleElem) = FreeModuleElem(-a.r, a.parent)

function check_parent(a::FreeModuleElem, b::FreeModuleElem)
  if parent(a) !== parent(b)
    error("elements not compatible")
  end  
end

function +(a::FreeModuleElem, b::FreeModuleElem)
   check_parent(a, b)
   return FreeModuleElem(a.r+b.r, a.parent)
end

function -(a::FreeModuleElem, b::FreeModuleElem)
    check_parent(a,b)
    return FreeModuleElem(a.r+b.r, a.parent)
end

function ==(a::FreeModuleElem, b::FreeModuleElem) 
    check_parent(a,b)
    return a.r == b.r
end

*(a::MPolyElem_dec, b::FreeModuleElem) = FreeModuleElem(a*b.r, b.parent)
*(a::MPolyElem, b::FreeModuleElem) = FreeModuleElem(a*b.r, b.parent)
*(a::Int, b::FreeModuleElem) = FreeModuleElem(a*b.r, b.parent)
*(a::Integer, b::FreeModuleElem) = FreeModuleElem(b.parent.R(a)*b.r, b.parent)
*(a::fmpq, b::FreeModuleElem) = FreeModuleElem(b.parent.R(a)*b.r, b.parent)
zero(F::FreeMod) = FreeModuleElem(sparse_row(F.R, Tuple{Int, elem_type(F.R)}[]), F)
parent(a::FreeModuleElem) = a.parent
iszero(a::FreeModuleElem) = Hecke.iszero(a.r)

mutable struct ModuleGens{T}
  O::Array{FreeModuleElem{T}, 1}
  S::Singular.smodule
  F::FreeMod
  SF::Singular.FreeMod

  function ModuleGens(O::Array{<:FreeModuleElem{T}, 1}) where {T}
    SF = singular_module(parent(O[1]))
    return ModuleGens(O, SF)
  end

  function ModuleGens(O::Array{<:FreeModuleElem{T}, 1}, F::FreeMod) where {T}
    SF = singular_module(F)
    return ModuleGens(O, F, SF)
  end

  function ModuleGens(O::Array{<:FreeModuleElem{T}, 1}, SF::Singular.FreeMod) where {T}
    return ModuleGens(O, parent(O[1]), SF)
  end

  function ModuleGens(O::Array{<:FreeModuleElem{T}, 1}, F::FreeMod, SF::Singular.FreeMod) where {T}
    r = new{T}()
    r.O = O
    r.SF = SF
    r.F = F
    return r
  end

  function ModuleGens(F::FreeMod{S}, s::Singular.smodule) where {S}
    r = new{S}()
    r.F = F
    if Singular.ngens(s) == 0
      r.SF = Singular.FreeModule(base_ring(s), 0)
    else
      r.SF = parent(s[1])
    end
    r.S = s
    #r.O = Array{FreeModuleElem_dec{S}, 1}(undef, Singular.ngens(s))
    r.O = [convert(F, s[i]) for i=1:Singular.ngens(s)]
    return r
  end
end

function show(io::IO, F::ModuleGens)
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

length(F::ModuleGens) = length(F.O)

function getindex(F::ModuleGens, ::Val{:O}, i::Int)
  if !isassigned(F.O, i)
    F.O[i] = convert(F.F, F.S[i])
  end
  return F.O[i]
end

function getindex(F::ModuleGens, ::Val{:S}, i::Int)
  if !isdefined(F, :S)
    F.S = Singular.smodule{elem_type(base_ring(F.SF))}(base_ring(F.SF), [convert(F.SF,x) for x = F.O]...)
  end
  return F.S[i]
end

function singular_assure(F::ModuleGens)
  if length(F) == 0 && !isdefined(F, :S)
    F.S = Singular.smodule{elem_type(base_ring(F.SF))}(base_ring(F.SF), F.SF())
    return
  end
  F[Val(:S), 1]
end

getindex(F::ModuleGens, i::Int) = getindex(F, Val(:O), i)

function singular_module(F::FreeMod)
  Sx = singular_ring(base_ring(F).R)
  return Singular.FreeModule(Sx, dim(F))
end

function convert(SF::Singular.FreeMod, m::FreeModuleElem)
  g = Singular.gens(SF)
  e = SF()
  Sx = base_ring(SF)
  for (p,v) = m.r
    e += Sx(v.f)*g[p]
  end
  return e
end

function convert(F::FreeMod, s::Singular.svector)
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
  return FreeModuleElem(sparse_row(base_ring(F), pv), F)
end

abstract type ModuleFPHom end
abstract type ModuleMap{T1, T2} <: Map{T1, T2, Hecke.HeckeMap, ModuleFPHom} end

mutable struct FreeModuleHom{T1, T2} <: ModuleMap{T1, T2} 
  header::MapHeader
  Hecke.@declare_other

  function FreeModuleHom(F::FreeMod{T}, G::S, a::Array{<:Any, 1}) where {T, S}
#    @assert isfiltrated(F) || all(ishomogenous, a) #neccessary and suffient according to Hans XXX
#same as non-homogenous elements are required, this too must not be enforced
    @assert all(x->parent(x) == G, a)
    @assert length(a) == ngens(F)
    #for filtrations, all is legal...
    r = new{typeof(F), typeof(G)}()
    function im_func(x::FreeModuleElem)
      b = zero(G)
      for (i,v) = x.r
        b += v*a[i]
      end
      return b
    end
    function pr_func(x::FreeModuleElem)
      @assert parent(x) == G
      c = coordinates(x, sub(G, a)) #??
      return FreeModuleElem(c, F)
    end
    function pr_func(x)
      @assert parent(x) == G
      #assume S == SubQuoElem_dec which cannot be asserted here, the type if defined too late
      c = coordinates(x.a, sub(G, a)) #??
      return FreeModuleElem(c, F)
    end
    r.header = MapHeader{typeof(F), typeof(G)}(F, G, im_func, pr_func)
    return r
  end
end
(h::FreeModuleHom)(a::FreeModuleElem) = image(h, a)

hom(F::FreeMod{T}, G, a) where {T} = FreeModuleHom(F, G, a)

function identity_map(M::ModuleFP)
  return hom(M, M, gens(M))
end

mutable struct SubQuo{T} <: ModuleFP{T}
  #meant to represent sub+ quo mod quo - as lazy as possible
  F::FreeMod{T}
  sub::ModuleGens
  quo::ModuleGens
  sum::ModuleGens
  std_sub::ModuleGens
  std_quo::ModuleGens
  AbstractAlgebra.@declare_other

  function SubQuo(F::FreeMod{R}, O::Array{<:FreeModuleElem, 1}) where {R}
    r = new{R}()
    r.F = F
    r.sub = ModuleGens(O, F, singular_module(F))
    r.sum = r.sub
    return r
  end
  function SubQuo(S::SubQuo, O::Array{<:FreeModuleElem{L}, 1}) where {L}
    r = new{L}()
    r.F = S.F
    r.sub = S.sub
    r.quo = ModuleGens(O, S.F, S.sub.SF)
    r.sum = ModuleGens(vcat(collect(r.sub), collect(r.quo)), S.F, S.sub.SF)
    return r
  end
  function SubQuo(F::FreeMod{R}, s::Singular.smodule) where {R}
    r = new{R}()
    r.F = F
    r.sub = ModuleGens(F, s)
    r.sum = r.sub
    if s.isGB
      r.std_sub = r.sub
    end
    return r
  end
  function SubQuo(F::FreeMod{R}, s::Singular.smodule, t::Singular.smodule) where {R}
    r = new{R}()
    r.F = F
    r.sub = ModuleGens(F, s)
    if s.isGB
      r.std_sub = r.sub
    end
    r.quo = ModuleGens(F, t)
    if t.isGB
      r.std_quo = r.quo
    end
    r.sum = ModuleGens(vcat(collect(r.sub), collect(r.quo)))
    return r
  end
end

function show(io::IO, SQ::SubQuo)
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
struct SubQuoElem{T} # this needs to be redone TODO
 a::FreeModuleElem{T}
 parent::SubQuo
end

elem_type(::SubQuo{T}) where {T} = SubQuoElem{T}
parent_type(::SubQuoElem{T}) where {T} = SubQuo{T}
elem_type(::Type{SubQuo{T}}) where {T} = SubQuoElem{T}
parent_type(::Type{SubQuoElem{T}}) where {T} = SubQuo{T}

function sum_gb_assure(SQ::SubQuo)
  singular_assure(SQ.sum)
  if SQ.sum.S.isGB
    return
  end
  SQ.sum = groebner_basis(SQ.sum)
end

function groebner_basis(F::ModuleGens)
  singular_assure(F)
  if F.S.isGB
    return F
  end
  return ModuleGens(F.F, Singular.std(F.S))
end

function show(io::IO, b::SubQuoElem)
  print(io, b.a)
end

parent(b::SubQuoElem) = b.parent

function (R::SubQuo)(a::FreeModuleElem; check::Bool = true)
  if check
    b = convert(R.sum.SF, a)
    sum_gb_assure(R)
    c = _reduce(b, R.sum.S)
    iszero(c) || error("not in the module")
  end
  return SubQuoElem(a, R)
end

function (R::SubQuo)(a::SubQuoElem)
  if parent(a) == R
    return a
  end
  error("illegal coercion")
end

+(a::SubQuoElem, b::SubQuoElem) = SubQuoElem(a.a+b.a, a.parent)
-(a::SubQuoElem, b::SubQuoElem) = SubQuoElem(a.a-b.a, a.parent)
-(a::SubQuoElem) = SubQuoElem(-a.a, a.parent)
*(a::MPolyElem_dec, b::SubQuoElem) = SubQuoElem(a*b.a, b.parent)
*(a::MPolyElem, b::SubQuoElem) = SubQuoElem(a*b.a, b.parent)
*(a::Int, b::SubQuoElem) = SubQuoElem(a*b.a, b.parent)
*(a::Integer, b::SubQuoElem) = SubQuoElem(a*b.a, b.parent)
*(a::fmpq, b::SubQuoElem) = SubQuoElem(a*b.a, b.parent)
==(a::SubQuoElem, b::SubQuoElem) = iszero(a-b)

function sub(F::FreeMod, O::Array{<:FreeModuleElem, 1})
  s = SubQuo(F, O)
end

function sub(F::FreeMod, O::Array{<:SubQuoElem, 1})
  return SubQuo(F, [x.a for x = O])
end

function sub(F::FreeMod, s::SubQuo)
  @assert !isdefined(s, :quo)
  return s
end

function sub(S::SubQuo, O::Array{<:SubQuoElem, 1})
  t = sub(S.F, O)
  if isdefined(S, :quo)
    return quo(t, collect(S.quo))
  else
    return t
  end
end

function quo(F::FreeMod, O::Array{<:FreeModuleElem, 1})
  S = SubQuo(F, basis(F))
  return SubQuo(S, O)
end

function quo(F::FreeMod, O::Array{<:SubQuoElem, 1})
  S = SubQuo(F, basis(F))
  return SubQuo(S, [x.a for x = O])
end

function quo(F::SubQuo, O::Array{<:FreeModuleElem, 1})
  @assert parent(O[1]) == F.F
  if isdefined(F, :quo)
    F.sub[Val(:S), 1]
    [F.quo[Val(:O), i] for i = 1:length(F.quo.O)]
    s = Singular.smodule{elem_type(base_ring(F.quo.SF))}(base_ring(F.quo.SF), [convert(F.quo.SF, x) for x = [O; F.quo.O]]...)
    return SubQuo(F.F, F.sub.S, s)
  end
  return SubQuo(F, O)
end

function quo(S::SubQuo, O::Array{<:SubQuoElem, 1})
  return SubQuo(S, [x.a for x = O])
end

function quo(S::SubQuo, T::SubQuo)
#  @assert !isdefined(T, :quo)
  return SubQuo(S, T.sum.O)
end

function quo(F::FreeMod, T::SubQuo)
  @assert !isdefined(T, :quo)
  return quo(F, gens(T))
end

function syzygie_module(F::ModuleGens; sub = 0)
  F[Val(:S), 1] #to force the existence of F.S
  s = Singular.syz(F.S)
  if sub !== 0
    G = sub
  else
    G = FreeMod(base_ring(F.F), length(F.O))
  end
  return SubQuo(G, s)
end

function gens(F::SubQuo)
  return map(x->SubQuoElem(x, F), collect(F.sub))
end

function gen(F::SubQuo, i::Int)
  return SubQuoElem(F.sub[Val(:O), i], F)
end

ngens(F::SubQuo) = length(F.sub)
base_ring(SQ::SubQuo) = base_ring(SQ.F)

zero(SQ::SubQuo) = SubQuoElem(zero(SQ.F), SQ)

function Base.iszero(F::SubQuo)
  return all(iszero, gens(F))
end

function Base.getindex(F::SubQuo, i::Int)
  i == 0 && return zero(F)
  return gen(F, i)
end

function Base.iterate(F::ModuleGens, i::Int = 1)
  if i>length(F)
    return nothing
  else
    return F[i], i+1
  end
end
Base.eltype(::ModuleGens{T}) where {T} = FreeModuleElem{T} 

#??? A scalar product....
function *(a::FreeModuleElem, b::Array{FreeModuleElem, 1})
  @assert dim(parent(a)) == length(b)
  s = zero(parent(a))
  for (p,v) = a.r
    s += v*b[p]
  end
  return s
end

function presentation(SQ::SubQuo)
  #A+B/B is generated by A and B
  #the relations are A meet B? written wrt to A
  s = syzygie_module(SQ.sum)
  #TODO: wait for Hans to release Modulo(A, B) that does exactly this
  c = collect(s.sub)
  R = base_ring(SQ)
  F = FreeMod(R, length(SQ.sub))
  q = elem_type(F)[]

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
    push!(q, FreeModuleElem(b, F))
  end
  #want R^a -> R^b -> SQ -> 0
  #TODO sort decoration and fix maps, same decoration should be bundled (to match pretty printing)
  G = FreeMod(R, length(s.sub))
  h_G_F = hom(G, F, q)
  #@assert iszero(h_G_F) || iszero(degree(h_G_F)) #???
  h_F_SQ = hom(F, SQ, gens(SQ))
  #@assert iszero(h_F_SQ) || iszero(degree(h_F_SQ)) #???
  Z = FreeMod(F.R, 0)
  Hecke.set_special(Z, :name => "Zero")
  h_SQ_Z = hom(SQ, Z, [zero(Z) for i=1:ngens(SQ)])
  return Hecke.ChainComplex(Oscar.ModuleFP, Oscar.ModuleMap[h_G_F, h_F_SQ, h_SQ_Z], check = false)
end

function presentation(F::FreeMod)
  Z = FreeMod(F.R, 0)
  Hecke.set_special(Z, :name => "Zero")
  return Hecke.ChainComplex(ModuleFP, ModuleMap[hom(Z, F, FreeModuleElem[]), hom(F, F, gens(F)), hom(F, Z, [zero(Z) for i=1:ngens(F)])], check = false)
end

mutable struct SubQuoHom{T1, T2} <: ModuleMap{T1, T2}
  header::Hecke.MapHeader
  im::Array{<:Any, 1}
  function SubQuoHom(D::SubQuo, C::ModuleFP, im::Array{<:Any, 1})
    first = true
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) == C, im)

    r = new{SubQuo, typeof(C)}()
    r.header = Hecke.MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im =  im
    return r
  end
end

function hom_tensor(G::ModuleFP, H::ModuleFP, A::Array{ <: ModuleMap, 1})
  tG = get_special(G, :tensor_product)
  tG === nothing && error("both modules must be tensor products")
  tH = get_special(H, :tensor_product)
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

function hom_prod_prod(G::ModuleFP, H::ModuleFP, A::Array{ <: ModuleMap, 2})
  tG = get_special(G, :tensor_product)
  tG === nothing && error("both modules must be direct products")
  tH = get_special(H, :tensor_product)
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

function coordinates(a::FreeModuleElem, SQ::SubQuo)
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
  b = ModuleGens([a], SQ.sum.SF)
  singular_assure(b)
  s, r = Singular.lift(S, b.S)
  if Singular.ngens(s) == 0 || iszero(s[1])
    error("elem not in module")
  end
  Sx = base_ring(SQ)
  Rx = Sx.R
  R = base_ring(Rx)
  return sparse_row(Rx, s[1], 1:ngens(SQ))
end


hom(D::SubQuo, C::ModuleFP, A::Array{<:Any, 1}) = SubQuoHom(D, C, A)

function image(f::SubQuoHom, a::SubQuoElem)
  i = zero(codomain(f))
  D = domain(f)
  b = coordinates(a.a, D)
  for (p,v) = b
    i += v*f.im[p]
  end
  return i
end

function image(f::SubQuoHom, a::FreeModuleElem)
  i = zero(codomain(f))
  D = domain(f)
  b = coordinates(a, D)
  for (p,v) = b
    i += v*f.im[p]
  end
  return i
end

function preimage(f::SubQuoHom, a::FreeModuleElem)
  @assert parent(a) == codomain(f)
  D = domain(f)
  i = zero(D)
  b = coordinates(a, image(f)[1])
  for (p,v) = b
    i += v*gen(D, p)
  end
  return i
end

(f::SubQuoHom)(a::FreeModuleElem) = image(f, a)
(f::SubQuoHom)(a::SubQuoElem) = image(f, a)

#ishom, homcompo missing

function iszero(a::SubQuoElem)
  C = parent(a)
  if !isdefined(C, :quo)
    return iszero(a.a)
  end
  if !isdefined(C, :std_quo)
    singular_assure(C.quo)
    C.std_quo = ModuleGens(C.quo.F, Singular.std(C.quo.S))
  end
  x = _reduce(convert(C.quo.SF, a.a), C.std_quo.S)
  return iszero(x)
end

function hom(F::FreeMod, G::FreeMod)
  @assert base_ring(F) == base_ring(G)
  GH = FreeMod(F.R, F.n * G.n)
  GH.S = [Symbol("($i -> $j)") for i = F.S for j = G.S]
  Hecke.set_special(GH, :show => Hecke.show_hom, :hom => (F, G))

  #list is g1 - f1, g2-f1, g3-f1, ...
  X = Hecke.MapParent(F, G, "homomorphisms")
  n = ngens(F)
  m = ngens(G)
  R = base_ring(F)
  function im(x::FreeModuleElem)
    return hom(F, G, [FreeModuleElem(x.r[R, (i-1)*m+1:i*m], G) for i=1:n])
  end
  function pre(h::FreeModuleHom)
    s = sparse_row(F.R)
    o = 0
    for i=1:n
      for (p,v) = h(gen(F, i)).r
        push!(s.pos, o+p)
        push!(s.values, v)
      end
      o += m
    end
    return FreeModuleElem(s, GH)
  end
  return GH, Hecke.MapFromFunc(im, pre, GH, X)
end

function kernel(h::FreeModuleHom)  #ONLY for free modules...
  G = domain(h)
  if ngens(G) == 0
    s = sub(G, gens(G))
    return s, hom(s, G, gens(G))
  end
  g = map(h, basis(G))
  if isa(codomain(h), SubQuo)
    g = [x.a for x = g]
    if isdefined(codomain(h), :quo)
      append!(g, collect(codomain(h).quo))
    end
  end
  #TODO allow sub-quo here as well
  b = ModuleGens(g)
  k = syzygie_module(b)
  if isa(codomain(h), SubQuo)
    s = collect(k.sub)
    k = sub(G, [FreeModuleElem(x.r[1:dim(G)], G) for x = s])
  else
    #the syzygie_module creates a new free module to work in
    k = sub(G, [FreeModuleElem(x.r, G) for x = collect(k.sub)])
  end
  @assert k.F == G
  c = collect(k.sub)
  return k, hom(k, parent(c[1]), c)
end

function image(h::FreeModuleHom)
  si = [x for x = map(h, basis(domain(h))) if !iszero(x)]
  s = sub(codomain(h), si)
  return s, hom(s, codomain(h), si)
end

function image(h::SubQuoHom)
  s = sub(codomain(h), h.im)
  return s, hom(s, codomain(h), h.im)
end

function kernel(h::SubQuoHom)
  D = domain(h)
  R = base_ring(D)
  F = FreeMod(R, ngens(D))
  hh = hom(F, codomain(h), map(h, gens(D)))
  k = kernel(hh)
  @assert domain(k[2]) == k[1]
  @assert codomain(k[2]) == F
  hh = hom(F, D, gens(D))
  im = map(x->hh(k[2](x)), gens(k[1]))
  k = sub(D, im)
  return k, hom(k, D, im)
end

function free_resolution(F::FreeMod)
  return presentation(F)
end

function free_resolution(S::SubQuo, limit::Int = -1)
  p = presentation(S)
  mp = [map(p, j) for j=1:length(p)]
  D = decoration(S)
  while true
    k, mk = kernel(mp[1])
    nz = findall(x->!iszero(x), gens(k))
    if length(nz) == 0 
      Z = FreeMod(base_ring(S), 0)
      Hecke.set_special(Z, :name => "Zero")
      h = hom(Z, domain(mp[1]), FreeModuleElem[])
      insert!(mp, 1, h)
      break
    elseif limit != -1 && length(mp) > limit
      break
    end
    F = FreeMod(base_ring(S), length(nz))
    g = hom(F, codomain(mk), collect(k.sub)[nz])
    insert!(mp, 1, g)
  end
  return Hecke.ChainComplex(ModuleFP, mp, check = false, direction = :right)
end

function Hecke.ring(I::MPolyIdeal)
  return parent(gen(I, 1))
end

function free_resolution(I::MPolyIdeal)
  F = free_module(Hecke.ring(I), 1)
  S = sub(F, [x * gen(F, 1) for x = gens(I)])
  n = Hecke.find_name(I)
  if n !== nothing
    AbstractAlgebra.set_name!(S, string(n))
  end
  return free_resolution(S)
end

function free_resolution(Q::MPolyQuo)
  F = free_module(Q.R, 1)
  q = quo(F, [x * gen(F, 1) for x = gens(Q.I)])
  n = Hecke.find_name(Q)
  if n !== nothing
    AbstractAlgebra.set_name!(q, String(n))
  end
  return free_resolution(q)
end

function iszero(f::ModuleMap)
  return all(iszero, map(f, gens(domain(f))))
end

function hom(M::ModuleFP, N::ModuleFP)
  p1 = presentation(M)
  p2 = presentation(N)
  k, mk = kernel(map(p2, 1))
  #Janko: have R^t1 -- g1 = map(p2, 0) -> R^t0 -> G
  #kernel g1: k -> R^t1
  #source: Janko's CA script: https://www.mathematik.uni-kl.de/~boehm/lehre/17_CA/ca.pdf
  F = FreeMod(base_ring(M), ngens(k))
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
  Hecke.set_special(H, :show => Hecke.show_hom, :hom => (M, N))

  #x in ker delta: mH_s0_t0(pro[1](x)) should be a hom from M to N
  function im(x::SubQuoElem)
    @assert parent(x) == H
    return hom(M, N, [map(p2, 2)(mH_s0_t0(pro[1](x.a))(preimage(map(p1, 2), g))) for g = gens(M)])
  end

  function pr(f::SubQuoHom)
    @assert domain(f) == M
    @assert codomain(f) == N
    Rs0 = domain(map(p1, 2))
    Rt0 = domain(map(p2, 2))
    g = hom(Rs0, Rt0, [preimage(map(p2, 2), f(map(p1, 2)(g))) for g = gens(Rs0)])

    return H(preimage(psi, (preimage(mH_s0_t0, g))).a)
    return SubQuoElem(emb[1](preimage(mH_s0_t0, g)), H) #???
  end
  return H, MapFromFunc(im, pr, H, Hecke.MapParent(M, N, "homomorphisms"))
end

#TODO
#  replace the +/- for the homs by proper constructors for homs and direct sums
#  relshp to store the maps elsewhere

function *(h::FreeModuleHom, g::FreeModuleHom) 
  @assert codomain(h) == domain(g)
  return hom(domain(h), codomain(g), [g(h(x)) for x = gens(domain(h))])
end
-(h::FreeModuleHom, g::FreeModuleHom) = hom(domain(h), codomain(h), [h(x) - g(x) for x = gens(domain(h))])
+(h::FreeModuleHom, g::FreeModuleHom) = hom(domain(h), codomain(h), [h(x) + g(x) for x = gens(domain(h))])

##################################################
# direct product
##################################################
function direct_product(F::FreeMod{T}...; task::Symbol = :sum) where {T}
  R = base_ring(F[1])
  G = FreeMod(R, sum([f.n for f = F]...))
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

function direct_product(G::ModuleFP...; task::Symbol = :none)
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
⊕(M::ModuleFP...) = direct_product(M..., task = :none)


function Hecke.canonical_injection(G::ModuleFP, i::Int)
  H = Hecke.get_special(G, :direct_product)
  if H === nothing
    error("module not a direct product")
  end
  0<i<= length(H) || error("index out of bound")
  j = i == 1 ? 0 : sum(ngens(H[l]) for l=1:i-1) -1
  return hom(H[i], G, [G[l+j] for l = 1:ngens(H[i])])
end

function Hecke.canonical_projection(G::ModuleFP, i::Int)
  H = Hecke.get_special(G, :direct_product)
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

function tensor_product(G::FreeMod...; task::Symbol = :none)
  s = G[1].S
  t = [[x] for x = 1:ngens(G[1])]
  for H = G[2:end]
    s = [Symbol("$x \\otimes $y") for x = s  for y = H.S]
    t = [push!(deepcopy(x), y) for x = t  for y = 1:ngens(H)]
  end

  F = FreeMod(G[1].R, prod([g.n for g in G]...))
  F.S = s
  Hecke.set_special(F, :show => Hecke.show_tensor_product, :tensor_product => G)
  if task == :none
    return F
  end

  function pure(g::FreeModuleElem...)
    @assert length(g) == length(G)
    @assert all(i-> parent(g[i]) == G[i], 1:length(G))
    z = [[x] for x = g[1].r.pos]
    zz = g[1].r.values
    for h = g[2:end]
      zzz = Array{Int, 1}[]
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
    return FreeModuleElem(sparse_row(F.R, [findfirst(x->x == y, t) for y = z], zz), F)
  end
  function pure(T::Tuple)
    return pure(T...)
  end
  function inv_pure(e::FreeModuleElem)
    if length(e.r.pos) == 0
      return Tuple(zero(g) for g = G)
    end
    @assert length(e.r.pos) == 1
    @assert isone(e.r.values[1])
    return Tuple(gen(G[i], t[e.r.pos[1]][i]) for i = 1:length(G))
  end

  return F, MapFromFunc(pure, inv_pure, Hecke.TupleParent(Tuple([g[0] for g = G])), F)
end

⊗(G::ModuleFP...) = tensor_product(G..., task = :none)

function free_module(F::FreeMod)
  return F
end

function free_module(F::SubQuo)
  return F.F
end

gens(F::FreeMod, G::FreeMod) = gens(F)
gens(F::SubQuo, G::FreeMod) = [x.a for x = gens(F)]
rels(F::FreeMod) = elem_type(F)[]
rels(F::SubQuo) = isdefined(F, :quo) ? collect(F.quo) : elem_type(F.F)[]

@doc Markdown.doc"""
    tensor_product(G::ModuleFP...; task::Symbol = :map) -> SubQuo_dec, Map

Given modules $G_i$ compute the tensor product $G_1\otimes \cdots \otimes G_n$.
If `task` is set to ":map", a map $\phi$ is returned that
maps tuples in $G_1 \times \cdots \times G_n$ to pure tensors
$g_1 \otimes \cdots \otimes g_n$. The map admits a preimage as well.
"""
function tensor_product(G::ModuleFP...; task::Symbol = :none)
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


#############################
function homology(C::Hecke.ChainComplex{ModuleFP})
  H = SubQuo[]
  for i=1:length(C)-1
    push!(H, quo(kernel(C.maps[i+1])[1], image(C.maps[i])[1])[1])
  end
  return H
end
#############################
#TODO move to Hecke
#  re-evaluate and use or not
function Base.getindex(r::Hecke.SRow, u::UnitRange)
  R = base_ring(r)
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

function getindex(a::Hecke.SRow, b::AbstractArray{Int, 1})
  if length(a.pos) == 0
    return a
  end
  m = minimum(b)
  b = sparse_row(parent(a.values[1]))
  for (k,v) = a
    if k in b
      push!(b.pos, k-b+1)
      push!(b.values, v)
    end
  end
  return b
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

#TODO: tensor_product from Raul's H is broken
