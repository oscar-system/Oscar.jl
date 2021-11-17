export FreeMod_dec, FreeModElem_dec, decoration

abstract type ModuleFP_dec{T} <: ModuleFP{T} end
abstract type AbstractFreeMod_dec{T} <: AbstractFreeMod{T} end
abstract type AbstractSubQuo_dec{T} <: AbstractSubQuo{T} end

abstract type AbstractFreeModElem_dec{T} <: AbstractFreeModElem{T} end
abstract type AbstractSubQuoElem_dec{T} <: AbstractSubQuoElem{T} end

const CRing_dec = Union{MPolyRing_dec, MPolyQuo{<:Oscar.MPolyElem_dec}}
const CRingElem_dec = Union{MPolyElem_dec, MPolyQuoElem{<:Oscar.MPolyElem_dec}}
#TODO: other name for CRing_dec -> which?


mutable struct FreeMod_dec{T <: CRingElem_dec} <: AbstractFreeMod_dec{T}
  F::FreeMod{T}
  d::Vector{GrpAbFinGenElem}

  AbstractAlgebra.@declare_other

  function FreeMod_dec{T}(R::CRing_dec,S::Vector{Symbol},d::Vector{GrpAbFinGenElem}) where T <: CRingElem_dec
    r = new{elem_type(R)}()
    r.F = FreeMod{T}(length(d),R,S)
    r.d = d
    return r
  end
end

function FreeMod_dec(R::CRing_dec, n::Int, name::String = "e"; cached::Bool = false) 
  return FreeMod_dec{elem_type(R)}(R, [Symbol("$name[$i]") for i=1:n], [decoration(R)[0] for i=1:n])
end
free_module_dec(R::CRing_dec, n::Int, name::String = "e"; cached::Bool = false) = FreeMod_dec(R, n, name, cached = cached)

# if one des not provide names for the generators, the standard names e_i are used for the unit vectors
function FreeMod_dec(R::CRing_dec, d::Vector{GrpAbFinGenElem}, name::String = "e"; cached::Bool = false) 
  return FreeMod_dec{elem_type(R)}(R, [Symbol("$name[$i]") for i=1:length(d)],d)
end
free_module_dec(R::CRing_dec, d::Vector{GrpAbFinGenElem}, name::String = "e"; cached::Bool = false) = FreeMod_dec(R, d, name, cached = cached)


function AbstractAlgebra.extra_name(F::FreeMod_dec)
  t = Hecke.get_special(F, :twist)
  if t !== nothing
    n = Hecke.get_special(t[1], :name)
    if n !== nothing
      return "$n($(t[2]))"
    end
  end
  if length(Set(F.d)) == 1
    n = Hecke.get_special(forget_decoration(F).R, :name)
    if n !== nothing
      return "$n^$(ngens(F))($(-F.d[1]))"
    end
  end
  return nothing
end

function show(io::IO, F::FreeMod_dec)
  @show_name(io, F)
  @show_special(io, F)

  print(io, "Decorated free module of rank $(rank(F)) over ")
  print(IOContext(io, :compact =>true), forget_decoration(F).R)

  i = 1
  while i < dim(F)
    d = F.d[i]
    j = 1
    while i+j <= dim(F) && d == F.d[i+j]
      j += 1
    end
    print(IOContext(io, :compact => true), forget_decoration(F).R, "^$j")
    print(IOContext(io, :compact => true), "(", -d, ")")
    if i+j < dim(F)
      print(io, " + ")
    end
    i += j
  end
end

function forget_decoration(F::FreeMod_dec)
  return F.F
end

@doc Markdown.doc"""
    base_ring(F::AbstractFreeMod)

Return the underlying ring of `F`.
"""
base_ring(F::FreeMod_dec) = forget_decoration(F).R

rank(F::FreeMod_dec) = rank(forget_decoration(F))

decoration(F::FreeMod_dec) = F.d
decoration(R::MPolyRing_dec) = R.D

# two free modules are equal if the rank and the ring are
function Base.:(==)(F::FreeMod_dec, G::FreeMod_dec)
  return forget_decoration(F) == forget_decoration(G) && F.d == G.d
end

# elements of free modules are encoded in SRows
struct FreeModElem_dec{T} <: AbstractFreeModElem_dec{T}
  coords::SRow{T} # also usable via coeffs()
  parent::FreeMod_dec{T}

  function FreeModElem_dec{T}(coords::SRow{T}, parent::FreeMod_dec{T}) where T
    r = new{T}(coords,parent)
    return r
  end
end

FreeModElem_dec(c::SRow{T}, parent::FreeMod_dec{T}) where T = FreeModElem_dec{T}(c, parent)

function FreeModElem_dec(c::Vector{T}, parent::FreeMod_dec{T}) where T
  @assert length(c) == rank(parent)
  sparse_coords = sparse_row(base_ring(parent), collect(1:rank(parent)), c)
  return FreeModElem_dec{T}(sparse_coords,parent)
end

function (F::FreeMod_dec{T})(c::SRow{T}) where T
  return FreeModElem_dec(c, F)
end

function (F::FreeMod_dec{T})(c::Vector{T}) where T 
  return FreeModElem_dec(c, F)
end

function (F::FreeMod_dec)()
  return FreeModElem_dec(sparse_row(base_ring(F)), F)
end

function FreeModElem(coords::SRow{T}, parent::FreeMod_dec{T}) where T <: CRingElem_dec
  return FreeModElem_dec{T}(coords, parent)
end


elem_type(::Type{FreeMod_dec{T}}) where {T} = FreeModElem_dec{T}
parent_type(::Type{FreeModElem_dec{T}}) where {T} = FreeMod_dec{T}
elem_type(::FreeMod_dec{T}) where {T} = FreeModElem_dec{T}
parent_type(::FreeModElem_dec{T}) where {T} = FreeMod_dec{T}


function generator_symbols(F::FreeMod_dec)
  return generator_symbols(forget_decoration(F))
end
@enable_all_show_via_expressify FreeModElem_dec


function degree_homogeneous_helper(u::FreeModElem_dec)
  if iszero(u)
    return nothing, true
  end
  first = true
  homogeneous_when_filtrated = true #this variable is only changed in the filtrated case
  F = parent(u)
  W = base_ring(F)
  ww = W.D[0]
  local w
  for (p,v) in coords(u)
    w = degree(v)+F.d[p]
    if first
      ww = w
      first = false
    elseif isgraded(W)
      if ww != w
        return nothing, false
      end
    else
      if ww != w
        homogeneous_when_filtrated = false
      end
      if W.lt(ww, w) 
        ww = w
      end
    end
  end
  return w, homogeneous_when_filtrated
end

function degree(a::FreeModElem_dec)
  if iszero(a)
    error("zero has no degree")
  end
  first = true
  F = parent(a)
  W = base_ring(F)
  ww = W.D[0]
  local w
  for (p,v) in coords(a)
    w = degree(v)+F.d[p]
    if first
      ww = w
      first = false
    elseif isgraded(W)
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

function homogeneous_components(a::FreeModElem_dec)
  res = Dict{GrpAbFinGenElem, FreeModElem_dec}()
  F = parent(a)
  for (p,v) in coords(a)
    c = homogeneous_components(v)
    for (pp, vv) in c
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

function homogeneous_component(a::FreeModElem_dec, g::GrpAbFinGenElem)
  F = parent(a)
  x = zero(F)
  for (p,v) in coords(a)
    x += homogeneous_component(v, g-F.d[p])*gen(F, p)
  end
  return x
end

function ishomogeneous(a::FreeModElem_dec)
  if iszero(a)
    return true
  end
  F = parent(a)
  first = true
  local d::GrpAbFinGenElem
  for (p,v) in coords(a)
    ishomogeneous(v) || return false
    if first
      d = F.d[p] + degree(v)
      first = false
    else
      F.d[p] + degree(v) == d || return false
    end
  end
  return true
end

# Weight vector or function?
# Should we already grade ModuleGens?
# Should it be possible to construct ungraded SubQuo with graded elements? (I.e. should the constructors
# accept AbstractFreeMod and AbstractFreeModElem instead of FreeMod and FreeModElem?)
# proceed with FreeModHom_dec?