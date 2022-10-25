export Localization, MPolyElemLoc, ideal, singular_assure, numerator,
       denominator, groebner_basis, minimal_generators

###############################################################################
# Constructors for localized polynomial ring and its elements                 #
###############################################################################

function _sort_helper(p::MPolyElem)
  return var_index(leading_term(p))
end

@Markdown.doc """
MPolyRingLoc{T} <: AbstractAlgebra.Ring where T <: AbstractAlgebra.FieldElem

The localization Pₘ = { f/g : f, g ∈ R, g ∉ 𝔪 } of a polynomial ring P = 𝕜[x₁,…,xₙ] over a 
base ring 𝕜 at a maximal ideal 𝔪 = ⟨x₁-a₁,…,xₙ-aₙ⟩ with coefficients aᵢ∈ 𝕜. 

The data being stored consists of 
  * a multivariate polynomial ring P = 𝕜[x₁,…,xₙ] with 𝕜 of type T;
  * the maximal ideal 𝔪 = ⟨x₁-a₁,…,xₙ-aₙ⟩⊂ P;
  * the number of variables n.
"""
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

@Markdown.doc """
MPolyElemLoc{T} <: AbstractAlgebra.RingElem where {T}

An element f/g in an instance of `MPolyRingLoc{T}`.

The data being stored consists of 
  * the fraction f/g as an instance of `AbstractAlgebra.Generic.Frac`;
  * the parent instance of `MPolyRingLoc{T}`.
"""
struct MPolyElemLoc{T} <: AbstractAlgebra.RingElem where {T}
  frac::AbstractAlgebra.Generic.Frac
  parent::MPolyRingLoc{T}

  # pass with checked = false to skip the non-trivial denominator check
  function MPolyElemLoc{T}(f::AbstractAlgebra.Generic.Frac,
                           p::MPolyRingLoc{T}, checked = true) where {T}
    R = p.base_ring
    B = base_ring(R)
    R != parent(numerator(f)) && error("Parent rings do not match")
    T != elem_type(B) && error("Type mismatch")
    if checked
      # this code seems to assume m.gens is of the form [xi - ai]_i
      m = p.max_ideal
      # This should be easier, somehow ...
      pt = leading_coefficient.([gen(R, i)-m.gens.O[i] for i in 1:nvars(R)])
      if evaluate(denominator(f), pt) == base_ring(R)(0)
        error("Element does not belong to the localization.")
      end
    end
    return new(f, p)
  end
end

function MPolyElemLoc(f::MPolyElem{T}, m::Oscar.MPolyIdeal) where {T}
  R = parent(f)
  return MPolyElemLoc{T}(f//R(1), Localization(R, m), false)
end

function MPolyElemLoc(f::AbstractAlgebra.Generic.Frac, m::Oscar.MPolyIdeal)
  R = parent(numerator(f))
  B = base_ring(R)
  return MPolyElemLoc{elem_type(B)}(f, Localization(R, m))
end

###############################################################################
# Basic functions                                                             #
###############################################################################

function Base.deepcopy_internal(a::MPolyElemLoc{T}, dict::IdDict) where T
  return MPolyElemLoc{T}(Base.deepcopy_internal(a.frac, dict), a.parent, false)
end

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

function check_parent(a::MPolyElemLoc{T}, b::MPolyElemLoc{T}, thr::Bool = true) where {T}
  if parent(a) == parent(b)
    return true
  elseif thr
     error("Parent rings do not match")
  else
    return false
  end
end


###############################################################################
# Arithmetics                                                                 #
###############################################################################

(W::MPolyRingLoc)() = W(W.base_ring())
(W::MPolyRingLoc)(i::Int) = W(W.base_ring(i))
(W::MPolyRingLoc)(i::RingElem) = W(W.base_ring(i))

function (W::MPolyRingLoc{T})(f::MPolyElem) where {T}
  return MPolyElemLoc{T}(f//one(parent(f)), W)
end

function (W::MPolyRingLoc{T})(g::AbstractAlgebra.Generic.Frac) where {T}
  return MPolyElemLoc{T}(g, W)
end

(W::MPolyRingLoc)(g::MPolyElemLoc) = W(g.frac)
Base.one(W::MPolyRingLoc) = W(1)
Base.zero(W::MPolyRingLoc) = W(0)

# Since a.parent.max_ideal is maximal and AA's frac arithmetic is reasonable,
# none of these ring operations should generate bad denominators.
# If this turns out to be a problem, remove the last false argument.
function +(a::MPolyElemLoc{T}, b::MPolyElemLoc{T}) where {T}
  check_parent(a, b)
  return MPolyElemLoc{T}(a.frac + b.frac, a.parent, false)
end

function -(a::MPolyElemLoc{T}, b::MPolyElemLoc{T}) where {T}
  check_parent(a, b)
  return MPolyElemLoc{T}(a.frac - b.frac, a.parent, false)
end

function -(a::MPolyElemLoc{T}) where {T}
  return MPolyElemLoc{T}(-a.frac, a.parent, false)
end

function *(a::MPolyElemLoc{T}, b::MPolyElemLoc{T}) where {T}
  check_parent(a, b)
  return MPolyElemLoc{T}(a.frac*b.frac, a.parent, false)
end

function ==(a::MPolyElemLoc{T}, b::MPolyElemLoc{T}) where {T}
  return check_parent(a, b, false) && a.frac == b.frac
end

function ^(a::MPolyElemLoc{T}, i::Int) where {T}
  return MPolyElemLoc{T}(a.frac^i, a.parent, false)
end

function Oscar.mul!(a::MPolyElemLoc, b::MPolyElemLoc, c::MPolyElemLoc)
  return b*c
end

function Oscar.addeq!(a::MPolyElemLoc, b::MPolyElemLoc)
  return a+b
end

function Base.:(//)(a::MPolyElemLoc{T}, b::MPolyElemLoc{T}) where {T}
  check_parent(a, b)
  return MPolyElemLoc{T}(a.frac//b.frac, a.parent)
end

###############################################################################
# Constructors for ideals                                                     #
###############################################################################

@Markdown.doc """
singular_ring_loc(R::MPolyRingLoc{T}; ord::Symbol = :negdegrevlex) where T

Sets up the singular ring in the backend to perform, for instance, standard basis computations
in `R`.
"""
function singular_ring_loc(R::MPolyRingLoc{T}; ord::Symbol = :negdegrevlex) where T
  return Singular.PolynomialRing(Oscar.singular_coeff_ring(base_ring(R.base_ring)),
              [string(x) for x = Nemo.symbols(R)],
              ordering = ord,
              cached = false)[1]
end

@Markdown.doc """
IdealGensLoc{S}

The main workhorse for translation of localizations of polynomial algebras at 
maximal ideals (i.e. instances of `MPolyRingLoc`) to the singular ring with 
local orderings in the backend.

This struct stores the following data:
  * an `MPolyRingLoc`, the parent ring Pₘ on the OSCAR side;
  * a `Singular.PolyRing`, the parent ring on the Singular side;
  * a `Vector{S}` of elements hᵢ = fᵢ/gᵢ ∈ Pₘ ;
  * a `Singular.sideal` storing only the numerator fᵢ for each one of the elements hᵢ above (after applying a geometric shift taking the maximal ideal 𝔪 to the ideal of the origin).
See also the various `*_assure` methods associated with this struct.
"""
mutable struct IdealGensLoc{S}
  O::Vector{S}
  S::Singular.sideal
  Ox::MPolyRingLoc
  Sx::Singular.PolyRing
  isGB::Bool

  function IdealGensLoc(Ox::MPolyRingLoc{T}, b::Singular.sideal) where {T}
    r = new{elem_type(Ox)}()
    r.S     = b
    r.Ox    = Ox
    r.Sx    = base_ring(b)
    R       = Ox.base_ring
    m       = Ox.max_ideal
    phi     = hom(R, R, m.gens.O)
    r.O     = Ox.(phi.([R(x) for x = gens(b)]))
    r.isGB  = false
    return r
  end
  function IdealGensLoc(a::Vector{T}; ord::Symbol = :negdegrevlex) where T <: MPolyElemLoc
    r = new{T}()
    r.O     = a
    r.Ox    = parent(a[1])
    r.Sx    = singular_ring_loc(r.Ox, ord = ord)
    r.isGB  = false
    return r
  end
end

@Markdown.doc """
MPolyIdealLoc{S} <: Ideal{S}

An ideal I in an instance of `MPolyRingLoc{S}`.

The data being stored consists of 
  * an instance of `IdealGensLoc{S}` for the set of generators of I
and some further fields used for caching.
"""
mutable struct MPolyIdealLoc{S} <: Ideal{S}
  gens::IdealGensLoc{S}
  min_gens::IdealGensLoc{S}
  gb::Dict{MonomialOrdering, IdealGensLoc{S}}
  dim::Int

  function MPolyIdealLoc(Ox::T, s::Singular.sideal) where {T <: MPolyRingLoc}
    r = new{elem_type(T)}()
    r.dim = -1 # not known
    r.gens = IdealGensLoc(Ox, s)
    r.gb = Dict()
    return r
  end
  function MPolyIdealLoc(B::IdealGensLoc{T}) where T
    r = new{T}()
    r.dim = -1
    r.gens = B
    r.gb = Dict()
    return r
  end
  function MPolyIdealLoc(g::Vector{T}) where {T <: MPolyElemLoc}
    r = new{T}()
    r.dim = -1 # not known
    r.gens = IdealGensLoc(g)
    r.gb = Dict()
    return r
  end
end

@enable_all_show_via_expressify MPolyIdealLoc

function AbstractAlgebra.expressify(a::MPolyIdealLoc; context = nothing)
  return Expr(:call, :ideal, [expressify(g, context = context) for g in collect(a.gens)]...)
end

###############################################################################
# Basic ideal functions                                                       #
###############################################################################

function Base.getindex(A::IdealGensLoc, ::Val{:S}, i::Int)
  if !isdefined(A, :S)
    A.S = Singular.Ideal(A.Sx, [A.Sx(x) for x = A.O])
  end
  return A.S[i]
end

function Base.getindex(A::IdealGensLoc, ::Val{:O}, i::Int)
  if !isassigned(A.O, i)
    A.O[i] = A.Ox(A.S[i])
  end
  return A.O[i]
end

function Base.length(A::IdealGensLoc)
  if isdefined(A, :S)
    return Singular.ngens(A.S)
  else
    return length(A.O)
  end
end

function Base.iterate(A::IdealGensLoc, s::Int = 1)
  if s > length(A)
    return nothing
  end
  return A[Val(:O), s], s + 1
end

Base.eltype(::IdealGensLoc{S}) where S = S

###############################################################################
# Ideal constructor functions                                                 #
###############################################################################

@Markdown.doc """
ideal(h::Vector{T}) where T <: MPolyElemLoc

Construct the ideal I = ⟨h₁,…,hᵣ⟩ generated by the elements hᵢ in `h`.
"""
function ideal(g::Vector{T}) where T <: MPolyElemLoc
  return MPolyIdealLoc(g)
end

function ideal(Rx::MPolyRingLoc, g::Vector)
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

function singular_assure(I::IdealGensLoc)
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

#= function groebner_assure(I::MPolyIdealLoc)
 =   if !isdefined(I, :gb)
 =     if !isdefined(I.gens, :S)
 =       singular_assure(I)
 =     end
 =     R = I.gens.Sx
 =     i = Singular.std(I.gens.S)
 =     I.gb = IdealGensLoc(I.gens.Ox, i)
 =   end
 = end
 =
 = function groebner_basis(I::MPolyIdealLoc; ordering::Symbol = :negdegrevlex)
 =   if ordering != :negdegrevlex
 =     B = IdealGensLoc(I.gens.O, ordering = ordering)
 =     singular_assure(B)
 =     R = B.Sx
 =     !Oscar.Singular.has_local_ordering(R) && error("The ordering has to be a local ordering.")
 =     i = Singular.std(B.S)
 =     I.gb = IdealGensLoc(I.gens.Ox, i)
 =   else
 =     groebner_assure(I)
 =   end
 =   return I.gb.O
 = end =#

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
      I.min_gens = IdealGensLoc(I.gb.Ox, sid)
    else
      sid = Singular.Ideal(I.gens.Sx, Singular.libSingular.idMinBase(I.gens.S.ptr, I.gens.Sx.ptr))
      I.min_gens = IdealGensLoc(I.gens.Ox, sid)
    end
  end
  return I.min_gens.O
end

function groebner_assure(I::MPolyIdealLoc, ordering::MonomialOrdering=negdegrevlex(gens(base_ring(I).base_ring)), complete_reduction::Bool = false)
    if get(I.gb, ordering, -1) == -1
        I.gb[ordering]  = groebner_basis(I.gens, ordering, complete_reduction)
    end

    return I.gb[ordering]
end

function groebner_basis(B::IdealGensLoc, ordering::MonomialOrdering, complete_reduction::Bool = false)
   singular_assure(B)
   R = B.Sx
   !Oscar.Singular.has_local_ordering(R) && error("The ordering has to be a local ordering.")
   I  = Singular.Ideal(R, gens(B.S)...)
   i  = Singular.std(I)
   BA = IdealGensLoc(B.Ox, i)
   BA.isGB  = true
   if isdefined(BA, :S)
       BA.S.isGB  = true
   end

   return BA
end

function groebner_basis(I::MPolyIdealLoc; ordering::MonomialOrdering = negdegrevlex(gens(base_ring(I).base_ring)), complete_reduction::Bool=false)
    groebner_assure(I, ordering, complete_reduction)
    return collect(I.gb[ordering])
end
