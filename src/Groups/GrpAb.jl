# additional functions for abelian groups (type `GrpAbFinGen`)

# Restrict a morphism f : G -> H to the surjective morphism g : G -> im(G)
function restrict_codomain(f::GrpAbFinGenMap)
  G = domain(f)
  H, Htocd = image(f, false)
  imgs = elem_type(H)[]
  for g in gens(G)
    fl, h = haspreimage(Htocd, f(g))
    @assert fl
    push!(imgs, h)
  end
  return hom(G, H, imgs)
end


# Compute whether `x` is contained in `G`, modulo natural embeddings.
# (This is analogous to `issubset`.)
is_element(x::GrpAbFinGenElem, G::GrpAbFinGen) = issubset(sub([x])[1], G)

function is_finiteorder(a::GrpAbFinGenElem)
  G, m = snf(a.parent)
  b = preimage(m, a)
  return any(i -> iszero(G.snf[i]) && !iszero(b[i]), 1:ngens(G))
end

function order(::Type{T}, a::GrpAbFinGenElem) where T <: IntegerUnion
   return T(order(a))
end

# order(x::GrpAbFinGenElem) = order(ZZRingElem, x) # see Hecke.jl/src/GrpAb/Elem.jl

#TODO: do we need `one`? (Hecke defines `is_one(x::GrpAbFinGen)`.)
# one(x::GrpAbFinGen) == zero(x)
# one(x::GrpAbFinGenElem) == zero(x)
# Base.inv(x::GrpAbFinGenElem) = -x
# x^n for GrpAbFinGenElem -> n*x

function comm(x::GrpAbFinGenElem, y::GrpAbFinGenElem)
  x.parent == y.parent || error("Elements must belong to the same group")
  zero(parent(x))
end

has_gens(G::GrpAbFinGen) = true

#TODO: improve?
small_generating_set(G::GrpAbFinGen) = gens(G)


################################################################################
#
#   Conjugacy Classes in GrpAbFinGen groups: just wrap elements
#
################################################################################

struct GrpAbFinGenConjClass{T<:GrpAbFinGen, S<:Union{GrpAbFinGenElem,GrpAbFinGen}} <: GroupConjClass{T, S}
   X::T
   repr::S
end

function Base.hash(C::GrpAbFinGenConjClass, h::UInt)
  return Base.hash(Representative(C), h)
end

function Base.show(io::IO, C::GrpAbFinGenConjClass)
  print(io, string(representative(C)), " ^ ", string(C.X))
end

==(a::GrpAbFinGenConjClass, b::GrpAbFinGenConjClass) = representative(a) == representative(b)

Base.length(::Type{T}, C::GrpAbFinGenConjClass) where T <: IntegerUnion = T(1)

Base.length(C::GrpAbFinGenConjClass) = ZZRingElem(1)


################################################################################
#
#   Conjugacy classes of elements
#
################################################################################

conjugacy_class(G::GrpAbFinGen, g::GrpAbFinGenElem) = GrpAbFinGenConjClass(G, g)

Base.eltype(::Type{GrpAbFinGenConjClass{T,S}}) where {T,S} = S

Base.rand(C::GrpAbFinGenConjClass) = representative(C)

Base.rand(rng::Random.AbstractRNG, C::GrpAbFinGenConjClass) = representative(C)

number_conjugacy_classes(G::GrpAbFinGen) = ZZRingElem(order(G))

number_conjugacy_classes(::Type{T}, G::GrpAbFinGen) where T <: IntegerUnion = T(order(G))

conjugacy_classes(G::GrpAbFinGen) = [GrpAbFinGenConjClass(G, x) for x in G]

is_conjugate(G::GrpAbFinGen, x::GrpAbFinGenElem, y::GrpAbFinGenElem) = (x == y)

function representative_action(G::GrpAbFinGen, x::GrpAbFinGenElem, y::GrpAbFinGenElem)
   x == y ? (true, zero(G)) : (false, nothing)
end


################################################################################
#
#   Conjugacy classes of subgroups
#
################################################################################

conjugacy_class(G::T, H::T) where T <: GrpAbFinGen = GrpAbFinGenConjClass(G, H)

# function Base.rand(C::GroupConjClass{S,T}) where S <: GrpAbFinGen where T <: GrpAbFinGen
#    return C.repr
# end

# function Base.rand(rng::Random.AbstractRNG, C::GroupConjClass{S,T}) where S <: GrpAbFinGen where T <: GrpAbFinGen
#    return C.repr
# end

function conjugacy_classes_subgroups(G::GrpAbFinGen)
   return [conjugacy_class(G, H) for (H, mp) in subgroups(G)]
end

#TODO: Does Hecke provide the following?
# subgroup_reps(G::GrpAbFinGen; order::ZZRingElem = ZZRingElem(-1))
# conjugacy_classes_maximal_subgroups(G::GrpAbFinGen)
# maximal_subgroup_reps(G::GrpAbFinGen)
# low_index_subgroup_reps(G::GrpAbFinGen, n::Int)

conjugate_group(G::GrpAbFinGen, x::GrpAbFinGenElem) = G

Base.:^(G::GrpAbFinGen, x::GrpAbFinGenElem) = G

function is_conjugate(G::GrpAbFinGen, H::GrpAbFinGen, K::GrpAbFinGen)
   if is_subgroup(H, K)[1] && is_subgroup(K, H)[1]
     return true
   else
     return false
   end
end

function representative_action(G::GrpAbFinGen, H::GrpAbFinGen, K::GrpAbFinGen)
   if is_subgroup(H, K)[1] && is_subgroup(K, H)[1]
     return true, zero(G)
   else
     return false, nothing
   end
end

is_conjugate_subgroup(G::T, U::T, V::T) where T <: GrpAbFinGen = is_subgroup(V, U)[1]


Base.IteratorSize(::Type{<:GrpAbFinGenConjClass}) = Base.HasLength()

Base.iterate(C::GrpAbFinGenConjClass) = iterate(C, 0)

function Base.iterate(C::GrpAbFinGenConjClass, state::Int)
  if state == 0
    return representative(C), 1
  else
    return nothing
  end
end


################################################################################
#
# Normal Structure
#
################################################################################

function core(G::GrpAbFinGen, H::GrpAbFinGen)
  flag, emb = is_subgroup(H, G)
  flag || error("H  must be a subgroup of G")
  return (H, emb)
end

function normalizer(G::GrpAbFinGen, H::GrpAbFinGen)
  is_subgroup(H, G)[1] || error("H  must be a subgroup of G")
  return (G, identity_map(G))
end

function normalizer(G::GrpAbFinGen, x::GrpAbFinGenElem)
  is_element(x, G) || throw(ArgumentError("x does not lie in G"))
  return (G, identity_map(G))
end

function normal_closure(G::GrpAbFinGen, H::GrpAbFinGen)
  flag, emb = is_subgroup(H, G)
  flag || error("H  must be a subgroup of G")
  return (H, emb)
end

pcore(G::GrpAbFinGen, p::IntegerUnion) = sylow_subgroup(G, p)


################################################################################
#
# Specific Subgroups
#
################################################################################

fitting_subgroup(G::GrpAbFinGen) = (G, identity_map(G))

function frattini_subgroup(G::GrpAbFinGen)
   is_finite(G) || throw(ArgumentError("G is not finite"))
   subgens = GrpAbFinGenElem[]
   for x in gens(G)
     for (p, e) in collect(factor(order(x)))
       x = p*x
     end
     if !is_zero(x)
       push!(subgens, x)
     end
   end
   return sub(G, subgens)
end

function socle(G::GrpAbFinGen)
   is_finite(G) || throw(ArgumentError("G is not finite"))
   subgens = GrpAbFinGenElem[]
   for x in gens(G)
     n = 1
     ord = order(x)
     for (p, e) in collect(factor(ord))
       n = p*n
     end
     x = divexact(ord, n)*x
     if !is_zero(x)
       push!(subgens, x)
     end
   end
   return sub(G, subgens)
end

solvable_radical(G::GrpAbFinGen) = (G, identity_map(G))


################################################################################
#
# Sylow & Hall Subgroups
#
################################################################################

#TODO: how to compute complements?
# complement_class_reps(G::T, N::T) where T <: GrpAbFinGen
# complement_system(G::GrpAbFinGen)

function sylow_subgroup(G::GrpAbFinGen, p::IntegerUnion)
   is_finite(G) || throw(ArgumentError("G is not finite"))
   is_prime(p) || throw(ArgumentError("p is not a prime"))
   subgens = GrpAbFinGenElem[]
   for x in gens(G)
     ord = order(x)
     while mod(ord, p) == 0
       ord = divexact(ord, p)
     end
     x = ord*x
     if !is_zero(x)
       push!(subgens, x)
     end
   end
   return sub(G, subgens)
end

function sylow_system(G::GrpAbFinGen)
   is_finite(G) || throw(ArgumentError("G is not finite"))
   result = GrpAbFinGen[]
   for (p, e) in collect(factor(order(G)))
     push!(result, sylow_subgroup(G, p)[1])
   end
   return result
end

# no longer documented, better use `hall_subgroup_reps`
hall_subgroup(G::GrpAbFinGen, P::AbstractVector{<:IntegerUnion}) = hall_subgroup_reps(G, P)[1]

function hall_subgroup_reps(G::GrpAbFinGen, P::AbstractVector{<:IntegerUnion})
   is_finite(G) || throw(ArgumentError("G is not finite"))
   P = unique(P)
   all(is_prime, P) || throw(ArgumentError("The integers must be prime"))
   subgens = GrpAbFinGenElem[]
   for x in gens(G)
     ord = order(x)
     for p in P
       ordp = ord
       while mod(ordp, p) == 0
         ordp = divexact(ordp, p)
       end
       px = ordp*x
       if !is_zero(px)
         push!(subgens, px)
       end
     end
   end
   return [sub(G, subgens)[1]]
end

function hall_system(G::GrpAbFinGen)
   is_solvable(G) || error("G must be solvable")
   primes = [p for (p, e) in factor(order(G))]
   result = GrpAbFinGen[]
   for P in subsets(Set(primes))
     push!(result, hall_subgroup_reps(G, collect(P))[1])
   end
   return result
end

################################################################################
#
#   Properties
#
################################################################################

is_almostsimple(G::GrpAbFinGen) = false

is_finitelygenerated(G::GrpAbFinGen) = true

is_perfect(G::GrpAbFinGen) = is_trivial(G)

is_pgroup(G::GrpAbFinGen) = is_pgroup_with_prime(G)[1]

is_quasisimple(G::GrpAbFinGen) = false

is_simple(G::GrpAbFinGen) = is_finite(G) && is_prime(order(G))

is_solvable(G::GrpAbFinGen) = true

is_sporadic_simple(G::GrpAbFinGen) = false

function is_pgroup_with_prime(::Type{T}, G::GrpAbFinGen) where T <: IntegerUnion
  is_trivial(G) && return true, nothing
  is_finite(G) || return false, nothing
  flag, p, e = is_prime_power_with_data(order(G))
  flag && return true, T(p)
  return false, nothing
end

is_pgroup_with_prime(G::GrpAbFinGen) = is_pgroup_with_prime(ZZRingElem, G)

nilpotency_class(G::GrpAbFinGen) = (order(G) == 1 ? 0 : 1)

# helper for prime_of_pgroup: this helper is efficient thanks to
# @gapattribute, but we also want prime_of_pgroup to accept an optional
# type argument; so we cannot use @gapattribute directly. Instead we set
# up this auxiliary _prime_of_pgroup which then is called by the real
# prime_of_pgroup.
# TODO: enhance @gapattribute so this is not necessary
function _prime_of_pgroup(G::GrpAbFinGen)
  flag, p, e = is_prime_power_with_data(order(G))
  if !flag
    error("only supported for non-trivial p-groups")
  end
  return p
end

function prime_of_pgroup(::Type{T}, G::GrpAbFinGen) where T <: IntegerUnion
  return ZZRingElem(_prime_of_pgroup(G))
end

prime_of_pgroup(G::GrpAbFinGen) = prime_of_pgroup(ZZRingElem, G)

#TODO There is no underlying GAP attribute, we do not store the prime.
#   has_prime_of_pgroup(G::GrpAbFinGen) = false
#   function set_prime_of_pgroup(G::GrpAbFinGen, p::IntegerUnion)
#     # do nothing
#   end
# The same holds for `fitting_subgroup`, `frattini_subgroup`, `socle`, ...
