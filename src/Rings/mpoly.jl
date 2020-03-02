#module MPolyModule

import AbstractAlgebra: PolyRing, PolynomialRing, total_degree, degree, Ideal,
                        MPolyElem, Generic.MPolyExponentVectors, Generic.MPolyCoeffs,
                        Generic.MPolyBuildCtx, Generic.push_term!, Generic.finish, MPolyRing,
                        base_ring, ngens, gens, dim, ordering, SetMap, Map
import Nemo                        
import Nemo: fmpz, fmpq

import Singular
import Base: +, *, ==, ^, -

import Hecke
import Hecke: MapHeader

export PolynomialRing, total_degree, degree, MPolyElem, ordering, ideal, groebner_basis

#TODO/ to think
#  default in Nemo is     :lex
#             Singular is :degrevlex -> better for std
#by default, use different orders???
#make BiPolyArray use different orders for both? make the type depend on it?
#
#for std: abstraction to allow Christian to be used
#
#type for orderings, use this...
#in general: all algos here needs revision: do they benefit from gb or not?

mutable struct BiPolyArray
  O::Array{MPolyElem, 1} 
  S::Singular.sideal
  Ox #Oscar Poly Ring
  Sx # Singular Poly Ring, poss. with different ordering
  function BiPolyArray(a::Array{T, 1}; keep_ordering::Bool = true) where {T <: MPolyElem}
    r = new()
    r.O = a
    r.Ox = parent(a[1])
    r.Sx = singular_ring(r.Ox, keep_ordering = keep_ordering)
    return r
  end
  function BiPolyArray(Ox, b::Singular.sideal)
    r =new()
    r.S = b
    r.O = Array{MPolyElem}(undef, Singular.ngens(b))
    r.Ox = Ox
    r.Sx = base_ring(r.S)
    return r
  end
end

function Base.getindex(A::BiPolyArray, ::Val{:S}, i::Int)
  if !isdefined(A, :S)
    A.S = Singular.Ideal(A.Sx, [convert(A.Sx, x) for x = A.O])
  end
  return A.S[i]
end

function Base.getindex(A::BiPolyArray, ::Val{:O}, i::Int)
  if !isassigned(A.O, i)
    A.O[i] = convert(A.Ox, A.S[i])
  end
  return A.O[i]
end

function Base.length(A::BiPolyArray)
  if isdefined(A, :S)
    return Singular.ngens(A.S)
  else
    return length(A.O)
  end
end

function Base.iterate(A::BiPolyArray, s::Int = 1)
  if s > length(A)
    return nothing
  end
  return A[Val(:O), s], s+1
end

function Base.convert(Ox::MPolyRing, f::MPolyElem) 
  O = base_ring(Ox)
  g = MPolyBuildCtx(Ox)
  for (c, e) = Base.Iterators.zip(MPolyCoeffs(f), MPolyExponentVectors(f))
    push_term!(g, O(c), e)
  end
  return finish(g)
end

function Base.convert(::Type{fmpz}, a::Singular.n_Z)
  return fmpz(BigInt(a))
end

function Base.convert(::Type{fmpq}, a::Singular.n_Q)
  return fmpq(Base.Rational{BigInt}(a))
end

function (::Nemo.FlintRationalField)(a::Singular.n_Q)
  return convert(fmpq, a)
end

function (S::Singular.Rationals)(a::fmpq)
  b = Base.Rational{BigInt}(a)
  return S(b)
end
(F::Singular.N_ZpField)(a::Nemo.gfp_elem) = F(lift(a))
(F::Nemo.GaloisField)(a::Singular.n_Zp) = F(Int(a))

singular_ring(::Nemo.FlintRationalField) = Singular.Rationals()
singular_ring(F::Nemo.GaloisField) = Singular.Fp(Int(characteristic(F)))

function singular_ring(Rx::Nemo.FmpqMPolyRing; keep_ordering::Bool = true)
  if keep_ordering
    return Singular.PolynomialRing(singular_ring(base_ring(Rx)), 
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = ordering(Rx),
              cached = false)[1]
  else
    return Singular.PolynomialRing(singular_ring(base_ring(Rx)), 
              [string(x) for x = Nemo.symbols(Rx)],
              cached = false)[1]
  end          
end

function singular_ring(Rx::Generic.MPolyRing{T}; keep_ordering::Bool = true) where {T <: RingElem}
  if keep_ordering
    return Singular.PolynomialRing(singular_ring(base_ring(Rx)), 
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = ordering(Rx),
              cached = false)[1]
  else
    return Singular.PolynomialRing(singular_ring(base_ring(Rx)), 
              [string(x) for x = Nemo.symbols(Rx)],
              cached = false)[1]
  end          
end

function singular_ring(Rx::Nemo.FmpqMPolyRing, ord::Symbol)
  return Singular.PolynomialRing(singular_ring(base_ring(Rx)), 
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = ord,
              cached = false)[1]
end

function singular_ring(Rx::Generic.MPolyRing{T}, ord::Symbol) where {T <: RingElem}
  return Singular.PolynomialRing(singular_ring(base_ring(Rx)), 
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = ord,
              cached = false)[1]
end

mutable struct MPolyIdeal <: Ideal{MPolyElem}
  gens::BiPolyArray
  gb::BiPolyArray
  dim::Int

  function MPolyIdeal(g::Array{T, 1}) where {T <: MPolyElem}
    r = new()
    r.dim = -1 #not known
    r.gens = BiPolyArray(g, keep_ordering = false)
    return r
  end
  function MPolyIdeal(Ox::MPolyRing, s::Singular.sideal)
    r = new()
    r.dim = -1 #not known
    r.gens = BiPolyArray(Ox, s)
    if s.isGB
      r.gb = gens
    end
    return r
  end
end

function Base.show(io::IO, I::MPolyIdeal)
  println(io, I.gens)
end

function ideal(g::Array{T, 1}) where {T <: MPolyElem}
  return MPolyIdeal(g)
end

function ideal(g::Array{Any, 1})
  return ideal(typeof(g[1])[x for x = g])
end

function ideal(Rx::MPolyRing, g::Array{<:Any, 1})
  f = elem_type(Rx)[Rx(f) for f = g]
  return ideal(f)
end

function singular_assure(I::MPolyIdeal)
  if !isdefined(I.gens, :S)
    I.gens.S = Singular.Ideal(I.gens.Sx, [convert(I.gens.Sx, x) for x = I.gens.O])
  end
end

function *(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal(I.gens.Ox, I.gens.S * J.gens.S)
end

function +(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal(I.gens.Ox, I.gens.S + J.gens.S)
end
-(I::MPolyIdeal, J::MPolyIdeal) = I+J

function ==(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return Singular.equal(I.gens.S, J.gens.S)
end

function ^(I::MPolyIdeal, j::Int)
  singular_assure(I)
  return MPolyIdeal(I.gens.Ox, I.gens.S^j)
end

function Base.intersect(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal(I.gens.Ox, Singular.intersection(I.gens.S, J.gens.S))
end

function ngens(I::MPolyIdeal)
  return length(I.gens)
end

function Base.issubset(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return Singular.contains(I.gens.S, J.gens.S)
end

function gens(I::MPolyIdeal)
  return [I.gens[Val(:O), i] for i=1:ngens(I)]
end

function groebner_assure(I::MPolyIdeal)
  if !isdefined(I, :gb)
    singular_assure(I)
    I.gb = BiPolyArray(I.gens.Ox, Singular.std(I.gens.S))
  end
end

function dim(I::MPolyIdeal)
  if I.dim > -1
    return I.dim
  end
  groebner_assure(I)
  I.dim = Singular.dimension(I.gb.S)
  return I.dim
end


function Base.in(f::MPolyElem, I::MPolyIdeal)
  groebner_assure(I)
  Sx = base_ring(I.gb.S)
  return Singular.iszero(reduce(convert(Sx, f), I.gb.S))
end

function base_ring(I::MPolyIdeal)
  return I.gens.Ox
end

function im_func(f::MPolyElem, S::MPolyRing, i::Array{Int, 1})
  O = base_ring(S)
  g = MPolyBuildCtx(S)
  for (c, e) = Base.Iterators.zip(MPolyCoeffs(f), MPolyExponentVectors(f))
    f = zeros(Int, nvars(S))
    for j=1:length(e)
      if i[j] == 0
        e[j] != 0 && error("illegal map: var $(j) is used")
      else
        f[i[j]] = e[j]
      end
    end
    push_term!(g, O(c), f)
  end
  return finish(g)
end


abstract type OscarMap <: SetMap end
mutable struct MPolyHom_vars{T1, T2}  <: Map{T1, T2, Hecke.HeckeMap, MPolyHom_vars}
  header::Hecke.MapHeader
  Hecke.@declare_other
  i::Array{Int, 1}

  function MPolyHom_vars{T1, T2}(R::T1, S::T2, i::Array{Int, 1}) where {T1 <: MPolyRing, T2 <: MPolyRing}
    r = new()
    p = sortperm(i)
    j = Int[]
    for h = 1:length(p)
      if i[p[h]] != 0
        j = p[h:length(p)]
        break
      end
    end
    r.header = MapHeader{T1, T2}(R, S, x -> im_func(x, S, i), y-> im_func(y, R, j))
    r.i = i
    return r
  end

  function MPolyHom_vars{T1, T2}(R::T1, S::T2; type::Symbol = :none) where {T1 <: MPolyRing, T2 <: MPolyRing}

    if type == :names
      i = Int[]
      for h = symbols(R)
        push!(i, findfirst(x -> x == h, symbols(S)))
      end
      return MPolyHom_vars(R, S, i)
    end
    error("type not supported")
  end
end

(f::MPolyHom_vars)(g::MPolyElem) = image(f, g)

function Hecke.hom(R::MPolyRing, S::MPolyRing, i::Array{Int, 1})
  return MPolyHom_vars{typeof(R), typeof(S)}(R, S, i)
end

function eliminate(I::MPolyIdeal, l::Array{<:MPolyElem, 1})
  #wrong!!! ideal intersection is not eliminate!!! 
  J = ideal(base_ring(I), [x for x in gens(base_ring(I)) if !(x in l)])
  return intersect(I, J)
end

function groebner_basis(I::MPolyIdeal)
  groebner_assure(I)
  return collect(I.gb)
end

function groebner_basis(I::MPolyIdeal, ord::Symbol)
  R = singular_ring(base_ring(I), ord)
  i = Singular.Ideal(R, [convert(R, x) for x = gens(I)])
  return collect(BiPolyArray(base_ring(I), i))
end

#end #MPolyModule
