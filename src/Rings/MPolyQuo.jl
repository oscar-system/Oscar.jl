##############################################################################
#
# quotient rings
#
##############################################################################

#TODO: add to singular_ring natively as this is potentially one
mutable struct MPolyQuo{S} <: AbstractAlgebra.Ring
  R::MPolyRing
  I::MPolyIdeal{S}
  AbstractAlgebra.@declare_other

  function MPolyQuo(R, I) where S
    @assert base_ring(I) == R
    r = new{elem_type(R)}()
    r.R = R
    r.I = I
    return r
  end
end

function show(io::IO, Q::MPolyQuo)
  Hecke.@show_name(io, Q)
  Hecke.@show_special(io, Q)
  io = IOContext(io, :compact => true)
  print(io, "Quotient of $(Q.R) by $(Q.I)")
end

gens(Q::MPolyQuo) = [Q(x) for x = gens(Q.R)]
ngens(Q::MPolyQuo) = ngens(Q.R)
gen(Q::MPolyQuo, i::Int) = Q(gen(Q.R, i))
Base.getindex(Q::MPolyQuo, i::Int) = Q(Q.R[i])

#TODO: think: do we want/ need to keep f on the Singular side to avoid conversions?
#      or use Bill's divrem to speed things up?
mutable struct MPolyQuoElem{S} <: RingElem
  f::S
  P::MPolyQuo{S}
end

AbstractAlgebra.expressify(a::MPolyQuoElem; context = nothing) = expressify(a.f, context = context)

function show(io::IO, a::MPolyQuoElem)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function singular_ring(Rx::MPolyQuo; keep_ordering::Bool = true)
  Sx = singular_ring(Rx.R, keep_ordering = keep_ordering)
  groebner_assure(Rx.I)
  Q = Sx(Singular.libSingular.rQuotientRing(Rx.I.gb.S.ptr, Sx.ptr))
  return Q
end

parent_type(::MPolyQuoElem{S}) where S = MPolyQuo{S}
parent_type(::Type{MPolyQuoElem{S}}) where S = MPolyQuo{S}
elem_type(::MPolyQuo{S})  where S= MPolyQuoElem{S}
elem_type(::Type{MPolyQuo{S}})  where S= MPolyQuoElem{S}

canonical_unit(a::MPolyQuoElem) = one(parent(a))

parent(a::MPolyQuoElem) = a.P

+(a::MPolyQuoElem, b::MPolyQuoElem) = MPolyQuoElem(a.f+b.f, a.P)
-(a::MPolyQuoElem, b::MPolyQuoElem) = MPolyQuoElem(a.f-b.f, a.P)
-(a::MPolyQuoElem) = MPolyQuoElem(-a.f, a.P)
*(a::MPolyQuoElem, b::MPolyQuoElem) = MPolyQuoElem(a.f*b.f, a.P)
^(a::MPolyQuoElem, b::Integer) = MPolyQuoElem(Base.power_by_squaring(a.f, b), a.P)

function Oscar.mul!(a::MPolyQuoElem, b::MPolyQuoElem, c::MPolyQuoElem)
  a.f = b.f*c.f
  return a
end

function Oscar.addeq!(a::MPolyQuoElem, b::MPolyQuoElem)
  a.f += b.f
  return a
end

function simplify!(a::MPolyQuoElem)
  R = parent(a)
  I = R.I
  groebner_assure(I)
  singular_assure(I.gb)
  Sx = base_ring(I.gb.S)
  I.gb.S.isGB = true
  f = a.f
  a.f = convert(I.gens.Ox, reduce(convert(Sx, f), I.gb.S))
  return a
end

function ==(a::MPolyQuoElem, b::MPolyQuoElem)
  simplify!(a)
  simplify!(b)
  return a.f == b.f
end

function quo(R::MPolyRing, I::MPolyIdeal) 
  q = MPolyQuo(R, I)
  function im(a::MPolyElem)
    return MPolyQuoElem(a, q)
  end
  function pr(a::MPolyQuoElem)
    return a.f
  end
  return q, MapFromFunc(im, pr, R, q)
end

lift(a::MPolyQuoElem) = a.f

(Q::MPolyQuo)() = MPolyQuoElem(Q.R(), Q)
(Q::MPolyQuo)(a::MPolyQuoElem) = a
(Q::MPolyQuo)(a) = MPolyQuoElem(Q.R(a), Q)

zero(Q::MPolyQuo) = Q(0)

one(Q::MPolyQuo) = Q(1)

#TODO: find a more descriptive, meaningful name
function _kbase(Q::MPolyQuo)
  I = Q.I
  groebner_assure(I)
  s = Singular.kbase(I.gb.S)
  if iszero(s)
    error("ideal was no zero-dimensional")
  end
  return [convert(Q.R, x) for x = gens(s)]
end

#TODO: the reverse map...
# problem: the "canonical" reps are not the monomials.
function vector_space(K::AbstractAlgebra.Field, Q::MPolyQuo)
  R = Q.R
  @assert K == base_ring(R)
  l = _kbase(Q)
  V = free_module(K, length(l))
  function im(a::Generic.FreeModuleElem)
    @assert parent(a) == V
    b = R(0)
    for i=1:length(l)
      c = a[i]
      if !iszero(c)
        b += c*l[i]
      end
    end
    return Q(b)
  end
  return MapFromFunc(im, V, Q)
end

################################################################################
#
#  To fix printing of fraction fields of MPolyQuo
#
################################################################################

function AbstractAlgebra.expressify(a::AbstractAlgebra.Generic.Frac{T};
                                    context = nothing) where {T <: MPolyQuoElem}
  n = numerator(a, false)
  d = denominator(a, false)
  if isone(d)
    return expressify(n, context = context)
  else
    return Expr(:call, ://, expressify(n, context = context), expressify(d, context = context))
  end
end
