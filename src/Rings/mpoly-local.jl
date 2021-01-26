export Localization, MPolyElemLoc, ideal, singular_assure, numerator,
       denominator, groebner_basis, minimal_generators

###############################################################################
# Constructors for localized polynomial ring and its elements                 #
###############################################################################

function _sort_helper(p::MPolyElem)
  return var_index(lt(p))
end

mutable struct MPolyRingLoc{T} <: AbstractAlgebra.Ring where T <: AbstractAlgebra.FieldElem
  base_ring::MPolyRing{T}
  max_ideal::Oscar.MPolyIdeal
  nvars::Int64

  function MPolyRingLoc(R::MPolyRing{S}, m::Oscar.MPolyIdeal) where {S}
    r = new{S}()
    r.base_ring = R
    r.max_ideal = ideal(R, sort(m.gens.O, by=_sort_helper)) # sorts maxideal to be of shape x_1-a_1, x_2-a_2, ..., x_n-a_n
    r.nvars = nvars(R)
    return r
  end
end

function Oscar.Localization(R::MPolyRing{S}, m::Oscar.MPolyIdeal) where S
  return MPolyRingLoc(R, m)
end

struct MPolyElemLoc{T} <: AbstractAlgebra.RingElem where {T}
  frac::AbstractAlgebra.Generic.Frac
  parent::MPolyRingLoc{T}
  function MPolyElemLoc(f::MPolyElem{T}, m::Oscar.MPolyIdeal) where {T}
    R = parent(f)
    (R != base_ring(m)) && error("Parent rings do not match!")
    r = new{T}(f//R(1), Localization(R, m))
    return r
  end
  function MPolyElemLoc(f::AbstractAlgebra.Generic.Frac, m::Oscar.MPolyIdeal)
    R = parent(numerator(f))
    B = base_ring(R)
    (R != base_ring(m)) && error("Parent rings do not match!")
    pt = lc.([gen(R, i)-m.gens.Ox[i] for i in 1:nvars(R)]) # This should be easier, somehow ...
    (evaluate(denominator(f), pt) == base_ring(R)(0)) && error("Element does not belong to the localization.")
    r = new{elem_type(B)}(f, Localization(R, m))
  end
end

###############################################################################
# Basic functions                                                             #
###############################################################################

function Base.show(io::IO, W::MPolyRingLoc)
  print("Localization of the ", W.base_ring, " at the maximal ", W.max_ideal)
end

function Base.show(io::IO, w::MPolyElemLoc)
  show(io, w.frac)
end

Nemo.symbols(R::MPolyRingLoc) = symbols(R.base_ring)
Nemo.nvars(R::MPolyRingLoc) = nvars(R.base_ring)
Nemo.parent(f::MPolyElemLoc) = f.parent
Nemo.numerator(f::MPolyElemLoc) = numerator(f.frac)
Nemo.denominator(f::MPolyElemLoc) = denominator(f.frac)

elem_type(::MPolyRingLoc{T}) where {T} = MPolyElemLoc{T}
elem_type(::Type{MPolyRingLoc{T}}) where {T} = MPolyElemLoc{T}
parent_type(::Type{MPolyElemLoc{T}}) where {T} = MPolyRingLoc{T}

###############################################################################
# Arithmetics                                                                 #
###############################################################################

(W::MPolyRingLoc)() = MPolyElemLoc(W.base_ring(), W.max_ideal)
(W::MPolyRingLoc)(i::Int) = MPolyElemLoc(W.base_ring(i), W.max_ideal)
(W::MPolyRingLoc)(i::RingElem) = MPolyElemLoc(W.base_ring(i), W.max_ideal)
(W::MPolyRingLoc)(f::MPolyElem) = MPolyElemLoc(f, W.max_ideal)
(W::MPolyRingLoc)(g::AbstractAlgebra.Generic.Frac) = MPolyElemLoc(g, W.max_ideal)
(W::MPolyRingLoc)(g::MPolyElemLoc) = W(g.frac)
Base.one(W::MPolyRingLoc) = MPolyElemLoc(one(W.base_ring), W.max_ideal)
Base.zero(W::MPolyRingLoc) = MPolyElemLoc(zero(W.base_ring), W.max_ideal)

+(a::MPolyElemLoc, b::MPolyElemLoc)   = MPolyElemLoc(a.frac+b.frac, a.parent.max_ideal)
-(a::MPolyElemLoc, b::MPolyElemLoc)   = MPolyElemLoc(a.frac-b.frac, a.parent.max_ideal)
-(a::MPolyElemLoc)   = MPolyElemLoc(-a.frac, a.parent.max_ideal)
*(a::MPolyElemLoc, b::MPolyElemLoc)   = MPolyElemLoc(a.frac*b.frac, a.parent.max_ideal)
==(a::MPolyElemLoc, b::MPolyElemLoc)   = a.frac == b.frac
^(a::MPolyElemLoc, i::Int)    = MPolyElemLoc(a.frac^i, a.parent.max_ideal)

function Oscar.mul!(a::MPolyElemLoc, b::MPolyElemLoc, c::MPolyElemLoc)
  return b*c
end

function Oscar.addeq!(a::MPolyElemLoc, b::MPolyElemLoc)
  return a+b
end

###############################################################################
# Constructors for ideals                                                     #
###############################################################################

function singular_ring_loc(R::MPolyRingLoc{T}; ord::Symbol = :negdegrevlex) where T
  return Singular.PolynomialRing(Oscar.singular_ring(base_ring(R.base_ring)),
              [string(x) for x = Nemo.symbols(R)],
              ordering = ord,
              cached = false)[1]
end

mutable struct BiPolyArrayLoc{S}
  O::Array{S, 1}
  S::Singular.sideal
  Ox::MPolyRingLoc
  Sx::Singular.PolyRing

  function BiPolyArrayLoc(Ox::MPolyRingLoc{T}, b::Singular.sideal) where {T}
    r = new{elem_type(Ox)}()
    r.S = b
    r.Ox = Ox
    r.Sx = base_ring(b)
    R = Ox.base_ring
    m = Ox.max_ideal
    phi = hom(R, R, m.gens.O)
    r.O = Ox.(phi.([R(x) for x = gens(b)]))
    return r
  end
  function BiPolyArrayLoc(a::Array{T, 1}; ord::Symbol = :negdegrevlex) where T <: MPolyElemLoc
    r = new{T}()
    r.O = a
    r.Ox = parent(a[1])
    r.Sx = singular_ring_loc(r.Ox, ord = ord)
    return r
  end
end

mutable struct MPolyIdealLoc{S} <: Ideal{S}
  gens::BiPolyArrayLoc{S}
  min_gens::BiPolyArrayLoc{S}
  gb::BiPolyArrayLoc{S}
  dim::Int

  function MPolyIdealLoc(Ox::T, s::Singular.sideal) where {T <: MPolyRingLoc}
    r = new{elem_type(T)}()
    r.dim = -1 # not known
    r.gens = BiPolyArrayLoc(Ox, s)
    if s.isGB
      r.gb = gens
    end
    return r
  end
  function MPolyIdealLoc(B::BiPolyArrayLoc{T}) where T
    r = new{T}()
    r.dim = -1
    r.gens = B
    return r
  end
  function MPolyIdealLoc(g::Array{T, 1}) where {T <: MPolyElemLoc}
    r = new{T}()
    r.dim = -1 # not known
    r.gens = BiPolyArrayLoc(g)
    return r
  end
end

###############################################################################
# Basic ideal functions                                                       #
###############################################################################

function Base.getindex(A::BiPolyArrayLoc, ::Val{:S}, i::Int)
  if !isdefined(A, :S)
    A.S = Singular.Ideal(A.Sx, [A.Sx(x) for x = A.O])
  end
  return A.S[i]
end

function Base.getindex(A::BiPolyArrayLoc, ::Val{:O}, i::Int)
  if !isassigned(A.O, i)
    A.O[i] = A.Ox(A.S[i])
  end
  return A.O[i]
end

function Base.length(A::BiPolyArrayLoc)
  if isdefined(A, :S)
    return Singular.ngens(A.S)
  else
    return length(A.O)
  end
end

function Base.iterate(A::BiPolyArrayLoc, s::Int = 1)
  if s > length(A)
    return nothing
  end
  return A[Val(:O), s], s + 1
end

Base.eltype(::BiPolyArrayLoc{S}) where S = S

function Base.show(io::IO, I::MPolyIdealLoc)
  print(io, "ideal generated by: ")
  g = collect(I.gens)
  first = true
  for i = g
    if first
      print(io, i)
      first = false
    else
      print(io, ", ", i)
    end
  end
  print(io, "")
end

function Base.show(io::IO, ::IJuliaMime, I::MPolyIdealLoc)
  print(io, "\$")
  math_html(io, I)
  print(io, "\$")
end

function math_html(io::IO, I::MPolyIdealLoc)
  print(io, "\\text{ideal generated by: }")
  g = collect(I.gens)
  first = true
  for i = g
    if first
      math_html(io, i)
      first = false
    else
      print(io, ", ")
      math_html(io, i)
    end
  end
  print(io, "")
end

###############################################################################
# Ideal constructor functions                                                 #
###############################################################################

function ideal(g::Array{T, 1}) where T <: MPolyElemLoc
  return MPolyIdealLoc(g)
end

function ideal(Rx::MPolyRingLoc, g::Array{<:Any, 1})
  f = elem_type(Rx)[Rx(f) for f = g]
  return ideal(f)
end

function ideal(Rx::MPolyRingLoc, g::Singular.sideal)
  return MPolyIdealLoc(Rx, g)
end

# Computes the Singular.jl data of an MPolyIdealLoc if it is not defined yet.
function singular_assure(I::MPolyIdealLoc)
  singular_assure(I.gens)
end

function singular_assure(I::BiPolyArrayLoc)
  if !isdefined(I, :S)
    R = I.Ox.base_ring
    m = I.Ox.max_ideal
    Q = I.Ox
    phi = hom(R, R, [2*gen(R, i)-m.gens.O[i] for i in 1:nvars(R)])
    I.S = Singular.Ideal(I.Sx, [I.Sx(phi(numerator(x))) for x = I.O])
  end
end

###############################################################################
# Ideal arithmetic                                                            #
###############################################################################

function Base.:*(I::MPolyIdealLoc, J::MPolyIdealLoc)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdealLoc(I.gens.Ox, I.gens.S * J.gens.S)
end

function Base.:+(I::MPolyIdealLoc, J::MPolyIdealLoc)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdealLoc(I.gens.Ox, I.gens.S + J.gens.S)
end
Base.:-(I::MPolyIdealLoc, J::MPolyIdealLoc) = I+J

function Base.:^(I::MPolyIdealLoc, j::Int)
  singular_assure(I)
  return MPolyIdealLoc(I.gens.Ox, I.gens.S^j)
end

###############################################################################
# Groebner bases                                                              #
###############################################################################

function base_ring(I::MPolyIdealLoc)
  return I.gens.Ox
end

function groebner_assure(I::MPolyIdealLoc)
  if !isdefined(I, :gb)
    if !isdefined(I.gens, :S)
      singular_assure(I)
    end
    R = I.gens.Sx
    i = Singular.std(I.gens.S)
    I.gb = BiPolyArrayLoc(I.gens.Ox, i)
  end
end

function groebner_basis(I::MPolyIdealLoc; ord::Symbol = :negdegrevlex)
  if ord != :negdegrevlex
    B = BiPolyArrayLoc(I.gens.O, ord = ord)
    singular_assure(B)
    R = B.Sx
    !Oscar.Singular.has_local_ordering(R) && error("The ordering has to be a local ordering.")
    i = Singular.std(B.S)
    I.gb = BiPolyArrayLoc(I.gens.Ox, i)
  else
    groebner_assure(I)
  end
  return I.gb.O
end

###############################################################################
# Ideal functions                                                             #
###############################################################################

function Base.:(==)(I::MPolyIdealLoc, J::MPolyIdealLoc)
  singular_assure(I)
  singular_assure(J)
  return Singular.equal(I.gens.S, J.gens.S)
end

function dim(I::MPolyIdealLoc)
  if I.dim > -1
    return I.dim
  end
  groebner_assure(I)
  I.dim = Singular.dimension(I.gb.S)
  return I.dim
end

function minimal_generators(I::MPolyIdealLoc)
  if !isdefined(I.gens, :S)
    singular_assure(I)
  end
  if !isdefined(I, :min_gens)
    if isdefined(I, :gb)
      sid = Singular.Ideal(I.gb.Sx, Singular.libSingular.idMinBase(I.gb.S.ptr, I.gb.Sx.ptr))
      I.min_gens = BiPolyArrayLoc(I.gb.Ox, sid)
    else
      sid = Singular.Ideal(I.gens.Sx, Singular.libSingular.idMinBase(I.gens.S.ptr, I.gens.Sx.ptr))
      I.min_gens = BiPolyArrayLoc(I.gens.Ox, sid)
    end
  end
  return I.min_gens.O
end
