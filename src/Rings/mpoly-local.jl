export Localization, MPolyElem_loc, ideal, singular_assure, oscar_assure,
       numerator, denominator, groebner_basis, minimal_generators

###############################################################################
# Constructors for localized polynomial ring and its elements                 #
###############################################################################

function _sort_helper(p::MPolyElem)
  return var_index(lt(p))
end

mutable struct MPolyRing_loc{T} <: AbstractAlgebra.Ring
  base_ring::MPolyRing{T}
  max_ideal::Oscar.MPolyIdeal
  nvars::Int64

  function MPolyRing_loc(R :: MPolyRing{S}, m::Oscar.MPolyIdeal) where {S}
    r = new{S}()
    r.base_ring = R
    r.max_ideal = ideal(R, sort(m.gens.O, by=_sort_helper)) # sorts maxideal to be of shape x_1-a_1, x_2-a_2, ..., x_n-a_n
    r.nvars = nvars(R)
    return r
  end
end

function Oscar.Localization(R::MPolyRing{S}, m::Oscar.MPolyIdeal) where S
  return MPolyRing_loc(R, m)
end

struct MPolyElem_loc{T} <: AbstractAlgebra.RingElem where {T}
  frac::AbstractAlgebra.Generic.Frac
  parent::MPolyRing_loc{T}
  function MPolyElem_loc(f::MPolyElem{T}, m::Oscar.MPolyIdeal) where {T}
    R = parent(f)
    (R != base_ring(m)) && error("Parent rings do not match!")
    r = new{T}(f//R(1), Localization(R, m))
    return r
  end
  function MPolyElem_loc(f::AbstractAlgebra.Generic.Frac, m::Oscar.MPolyIdeal)
    R = parent(numerator(f))
    B = base_ring(R)
    (R != base_ring(m)) && error("Parent rings do not match!")
    pt = lc.([gen(R, i)-m.gens.Ox[i] for i in 1:nvars(R)]) # This should be easier, somehow ...
    (evaluate(denominator(f), pt) == base_ring(R)(0)) && error("Element does not belong to the localization.")i
    r = new{elem_type(B)}(f, Localization(R, m))
  end
end

###############################################################################
# Basic functions                                                             #
###############################################################################

function Base.show(io::IO, W::MPolyRing_loc)
  print("Localization of the ", W.base_ring, " at the maximal ", W.max_ideal)
end

function Base.show(io::IO, w::MPolyElem_loc)
  show(io, w.frac)
end

Nemo.symbols(R::MPolyRing_loc) = symbols(R.base_ring)
Nemo.nvars(R::MPolyRing_loc) = nvars(R.base_ring)
Nemo.parent(f::MPolyElem_loc) = f.parent
Nemo.numerator(f::MPolyElem_loc) = numerator(f.frac)
Nemo.denominator(f::MPolyElem_loc) = denominator(f.frac)

elem_type(::MPolyRing_loc{T}) where {T} = MPolyElem_loc{T}
elem_type(::Type{MPolyRing_loc{T}}) where {T} = MPolyElem_loc{T}
parent_type(::Type{MPolyElem_loc{T}}) where {T} = MPolyRing_loc{T}

###############################################################################
# Arithmetics                                                                 #
###############################################################################

(W::MPolyRing_loc)() = MPolyElem_loc(W.base_ring(), W.max_ideal)
(W::MPolyRing_loc)(i::Int) = MPolyElem_loc(W.base_ring(i), W.max_ideal)
(W::MPolyRing_loc)(i::RingElem) = MPolyElem_loc(W.base_ring(i), W.max_ideal)
(W::MPolyRing_loc)(f::MPolyElem) = MPolyElem_loc(f, W.max_ideal)
(W::MPolyRing_loc)(g::AbstractAlgebra.Generic.Frac) = MPolyElem_loc(g, W.max_ideal)
(W::MPolyRing_loc)(g::MPolyElem_loc) = W(g.frac)
Base.one(W::MPolyRing_loc) = MPolyElem_loc(one(W.base_ring), W.max_ideal)
Base.zero(W::MPolyRing_loc) = MPolyElem_loc(zero(W.base_ring), W.max_ideal)

+(a::MPolyElem_loc, b::MPolyElem_loc)   = MPolyElem_loc(a.frac+b.frac, a.parent.max_ideal)
-(a::MPolyElem_loc, b::MPolyElem_loc)   = MPolyElem_loc(a.frac-b.frac, a.parent.max_ideal)
-(a::MPolyElem_loc)   = MPolyElem_loc(-a.frac, a.parent.max_ideal)
*(a::MPolyElem_loc, b::MPolyElem_loc)   = MPolyElem_loc(a.frac*b.frac, a.parent.max_ideal)
==(a::MPolyElem_loc, b::MPolyElem_loc)   = a.frac == b.frac
^(a::MPolyElem_loc, i::Int)    = MPolyElem_loc(a.frac^i, a.parent.max_ideal)

function Oscar.mul!(a::MPolyElem_loc, b::MPolyElem_loc, c::MPolyElem_loc)
  return b*c
end

function Oscar.addeq!(a::MPolyElem_loc, b::MPolyElem_loc)
  return a+b
end

###############################################################################
# Constructors for ideals                                                     #
###############################################################################

function singular_ring_loc(R::MPolyRing_loc{T}; ord::Symbol = :negdegrevlex) where T
  return Singular.PolynomialRing(Oscar.singular_ring(base_ring(R.base_ring)),
              [string(x) for x = Nemo.symbols(R)],
              ordering = ord,
              cached = false)[1]
end

mutable struct BiPolyArray_loc{S}
  O::Array{S, 1}
  S::Singular.sideal
  Ox
  Sx

  function BiPolyArray_loc(Ox::T, b::Singular.sideal) where {T <: MPolyRing_loc}
    r = new{elem_type(T)}()
    r.S = b
    r.Ox = Ox
    r.Sx = base_ring(b)
    R = Ox.base_ring
    m = Ox.max_ideal
    phi = hom(R, R, m.gens.O)
    r.O = Ox.(phi.([convert(R, x) for x = gens(b)]))
    return r
  end
  function BiPolyArray_loc(a::Array{T, 1}; ord::Symbol = :negdegrevlex) where T <: MPolyElem_loc
    r = new{T}()
    r.O = a
    r.Ox = parent(a[1])
    r.Sx = singular_ring_loc(r.Ox, ord = ord)
    return r
  end
end

mutable struct MPolyIdeal_loc{S} <: Ideal{S}
  gens::BiPolyArray_loc{S}
  min_gens::BiPolyArray_loc{S}
  gb::BiPolyArray_loc{S}
  dim::Int

  function MPolyIdeal_loc(Ox::T, s::Singular.sideal) where {T <: MPolyRing_loc}
    r = new{elem_type(T)}()
    r.dim = -1 #not known
    r.gens = BiPolyArray_loc(Ox, s)
    if s.isGB
      r.gb = gens
    end
    return r
  end
  function MPolyIdeal_loc(B::BiPolyArray_loc{T}) where T
    r = new{T}()
    r.dim = -1
    r.gens = B
    return r
  end
  function MPolyIdeal_loc(g::Array{T, 1}) where {T <: MPolyElem_loc}
    r = new{T}()
    r.dim = -1 #not known
    r.gens = BiPolyArray_loc(g)
    return r
  end
end

###############################################################################
# Basic ideal functions                                                       #
###############################################################################

function Base.getindex(A::BiPolyArray_loc, ::Val{:S}, i::Int)
  if !isdefined(A, :S)
    A.S = Singular.Ideal(A.Sx, [convert(A.Sx, x) for x = A.O])
  end
  return A.S[i]
end

function Base.getindex(A::BiPolyArray_loc, ::Val{:O}, i::Int)
  if !isassigned(A.O, i)
    A.O[i] = convert(A.Ox, A.S[i])
  end
  return A.O[i]
end

function Base.length(A::BiPolyArray_loc)
  if isdefined(A, :S)
    return Singular.ngens(A.S)
  else
    return length(A.O)
  end
end

function Base.iterate(A::BiPolyArray_loc, s::Int = 1)
  if s > length(A)
    return nothing
  end
  return A[Val(:O), s], s+1
end

Base.eltype(::BiPolyArray_loc{S}) where S = S

function Base.show(io::IO, I::MPolyIdeal_loc)
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

function Base.show(io::IO, ::IJuliaMime, I::MPolyIdeal_loc)
  print(io, "\$")
  math_html(io, I)
  print(io, "\$")
end

function math_html(io::IO, I::MPolyIdeal_loc)
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

function ideal(g::Array{T, 1}) where T <: MPolyElem_loc
  return MPolyIdeal_loc(g)
end

function ideal(Rx::MPolyRing_loc, g::Array{<:Any, 1})
  f = elem_type(Rx)[Rx(f) for f = g]
  return ideal(f)
end

function ideal(Rx::MPolyRing_loc, g::Singular.sideal)
  return MPolyIdeal_loc(Rx, g)
end

function singular_assure(I::MPolyIdeal_loc)
  singular_assure(I.gens)
end

function singular_assure(I::BiPolyArray_loc)
  if !isdefined(I, :S)
    R = I.Ox.base_ring
    m = I.Ox.max_ideal
    Q = I.Ox
    phi = hom(R, R, [2*gen(R, i)-m.gens.O[i] for i in 1:nvars(R)])
    I.S = Singular.Ideal(I.Sx, [convert(I.Sx, phi(numerator(x))) for x = I.O])
  end
end

#function oscar_assure(I::MPolyIdeal_loc)
#  if !isdefined(I.gens, :O)
#    R = I.gens.Ox.base_ring
#    m = I.gens.Ox.max_ideal
#    Q = I.gens.Ox
#    phi = hom(R, R, m.gens.O)
#    I.gens.O = Q.(phi.([convert(I.gens.Ox.base_ring, x) for x = gens(I.gens.S)]))
#  end
#end

###############################################################################
# Ideal arithmetic                                                            #
###############################################################################

function Base.:*(I::MPolyIdeal_loc, J::MPolyIdeal_loc)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal_loc(I.gens.Ox, I.gens.S * J.gens.S)
end

function Base.:+(I::MPolyIdeal_loc, J::MPolyIdeal_loc)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal_loc(I.gens.Ox, I.gens.S + J.gens.S)
end
Base.:-(I::MPolyIdeal_loc, J::MPolyIdeal_loc) = I+J

function Base.:^(I::MPolyIdeal_loc, j::Int)
  singular_assure(I)
  return MPolyIdeal_loc(I.gens.Ox, I.gens.S^j)
end

###############################################################################
# Groebner bases                                                              #
###############################################################################

function base_ring(I::MPolyIdeal_loc)
  return I.gens.Ox
end

function groebner_assure(I::MPolyIdeal_loc)
  if !isdefined(I, :gb)
    if !isdefined(I.gens, :S)
      singular_assure(I)
    end
    R = I.gens.Sx
    i = Singular.std(I.gens.S)
    I.gb = BiPolyArray_loc(I.gens.Ox, i)
  end
end

function groebner_basis(I::MPolyIdeal_loc; ord::Symbol = :negdegrevlex)
  if ord != :negdegrevlex
    B = BiPolyArray_loc(I.gens.O, ord = ord)
    singular_assure(B)
    R = B.Sx
    !Oscar.Singular.has_local_ordering(R) && error("The ordering has to be a local ordering.")
    i = Singular.std(B.S)
    I.gb = BiPolyArray_loc(I.gens.Ox, i)
  else
    groebner_assure(I)
  end
  return I.gb.O
end

###############################################################################
# Ideal functions                                                             #
###############################################################################

function Base.:(==)(I::MPolyIdeal_loc, J::MPolyIdeal_loc)
  singular_assure(I)
  singular_assure(J)
  return Singular.equal(I.gens.S, J.gens.S)
end

function dim(I::MPolyIdeal_loc)
  if I.dim > -1
    return I.dim
  end
  groebner_assure(I)
  I.dim = Singular.dimension(I.gb.S)
  return I.dim
end

function minimal_generators(I::MPolyIdeal_loc)
  if !isdefined(I.gens, :S)
    singular_assure(I)
  end
  if !isdefined(I, :min_gens)
    if isdefined(I, :gb)
      sid = Singular.Ideal(I.gb.Sx, Singular.libSingular.idMinBase(I.gb.S.ptr, I.gb.Sx.ptr))
      I.min_gens = BiPolyArray_loc(I.gb.Ox, sid)
    else
      sid = Singular.Ideal(I.gens.Sx, Singular.libSingular.idMinBase(I.gens.S.ptr, I.gens.Sx.ptr))
      I.min_gens = BiPolyArray_loc(I.gens.Ox, sid)
    end
  end
  return I.min_gens.O
end
