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

function order(::Type{T}, G::GrpAbFinGen) where T <: IntegerUnion
   return T(order(G))
end

#TODO: do we need `one`? (Hecke defines `is_one(x::GrpAbFinGen)`.)
# one(x::GrpAbFinGen) == zero(x)
# one(x::GrpAbFinGenElem) == zero(x)
# Base.inv(x::GrpAbFinGenElem) = -x
# x^n for GrpAbFinGenElem -> n*x

Base.:^(g::GrpAbFinGenElem, x::GrpAbFinGenElem) = g

function comm(x::GrpAbFinGenElem, y::GrpAbFinGenElem)
  @req (x.parent == y.parent) "elements must belong to the same group"
  return zero(parent(x))
end

has_gens(G::GrpAbFinGen) = true

function small_generating_set(G::GrpAbFinGen)
   is_snf(G) && return gens(G)
   S, mp = snf(G)
   return [mp(x) for x in gens(S)]
end

function index(::Type{I}, G::T, H::T) where I <: IntegerUnion where T <: GrpAbFinGen
   @req is_subgroup(H, G)[1] "H must be a subgroup of G"
   f = count(x -> x == 0, snf(G)[1].snf) - count(x -> x == 0, snf(H)[1].snf)
   @req f == 0 "index is supported only for subgroups of finite index"
   return I(divexact(order(torsion_subgroup(G)[1]), order(torsion_subgroup(H)[1]), check = false))
end


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

number_conjugacy_classes(G::GrpAbFinGen) = order(ZZRingElem, G)

number_conjugacy_classes(::Type{T}, G::GrpAbFinGen) where T <: IntegerUnion = order(T, G)

conjugacy_classes(G::GrpAbFinGen) = [GrpAbFinGenConjClass(G, x) for x in G]

is_conjugate(G::GrpAbFinGen, x::GrpAbFinGenElem, y::GrpAbFinGenElem) = (x == y)

function is_conjugate_with_data(G::GrpAbFinGen, x::GrpAbFinGenElem, y::GrpAbFinGenElem)
   x == y ? (true, zero(G)) : (false, nothing)
end

################################################################################
#
#   Conjugacy classes of subgroups
#
################################################################################

conjugacy_class(G::T, H::T) where T <: GrpAbFinGen = GrpAbFinGenConjClass(G, H)

function conjugacy_classes_subgroups(G::GrpAbFinGen)
   @req is_finite(G) "G is not finite"
   return [conjugacy_class(G, H) for (H, mp) in subgroups(G)]
end

function subgroup_reps(G::GrpAbFinGen; order::ZZRingElem = ZZRingElem(-1))
   if order > 0 && mod(Hecke.order(G), order) != 0
     # `subgroups` would throw an error
     return GrpAbFinGen[]
   end
   return [H for (H, mp) in subgroups(G, order = order)]
end

function low_index_subgroup_reps(G::GrpAbFinGen, n::Int)
   @req (n > 0) "index must be positive"
   res = [G]
   ord = order(G)
   for i in 2:n
     if mod(ord, i) == 0
       append!(res, [H for (H, mp) in subgroups(G, index = i)])
     end
   end
   return res
end

function maximal_subgroup_reps(G::GrpAbFinGen)
   @req is_finite(G) "G is not finite"
   primes = [p for (p, e) in factor(order(G))]
   res = typeof(G)[]
   for p in primes
     append!(res, [H for (H, mp) in subgroups(G, index = p)])
   end
   return res
end

function conjugacy_classes_maximal_subgroups(G::GrpAbFinGen)
   return [conjugacy_class(G, H) for H in maximal_subgroup_reps(G)]
end

conjugate_group(G::GrpAbFinGen, x::GrpAbFinGenElem) = G

Base.:^(G::GrpAbFinGen, x::GrpAbFinGenElem) = G

function is_conjugate(G::GrpAbFinGen, H::GrpAbFinGen, K::GrpAbFinGen)
   return is_subgroup(H, G)[1] && is_subgroup(K, G)[1] &&
          is_subgroup(H, K)[1] && is_subgroup(K, H)[1]
end

function is_conjugate_with_data(G::GrpAbFinGen, H::GrpAbFinGen, K::GrpAbFinGen)
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
  @req flag "H  must be a subgroup of G"
  return (H, emb)
end

function normalizer(G::GrpAbFinGen, H::GrpAbFinGen)
  @req is_subgroup(H, G)[1] "H  must be a subgroup of G"
  return (G, identity_map(G))
end

function normalizer(G::GrpAbFinGen, x::GrpAbFinGenElem)
  @req is_element(x, G) "x does not lie in G"
  return (G, identity_map(G))
end

function normal_closure(G::GrpAbFinGen, H::GrpAbFinGen)
  flag, emb = is_subgroup(H, G)
  @req flag "H  must be a subgroup of G"
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
   @req is_finite(G) "G is not finite"
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
   @req is_finite(G) "G is not finite"
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

trivial_subgroup(G::GrpAbFinGen) = sub(G, GrpAbFinGenElem[])

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
   @req is_finite(G) "G is not finite"
   @req is_prime(p) "p is not a prime"
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
   @req is_finite(G) "G is not finite"
   result = GrpAbFinGen[]
   for (p, e) in collect(factor(order(G)))
     push!(result, sylow_subgroup(G, p)[1])
   end
   return result
end

# no longer documented, better use `hall_subgroup_reps`
hall_subgroup(G::GrpAbFinGen, P::AbstractVector{<:IntegerUnion}) = hall_subgroup_reps(G, P)[1]

function hall_subgroup_reps(G::GrpAbFinGen, P::AbstractVector{<:IntegerUnion})
   @req is_finite(G) "G is not finite"
   P = unique(P)
   @req all(is_prime, P) "The integers must be prime"
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
   @req is_solvable(G) "G must be solvable"
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
  flag, _, p = is_prime_power_with_data(order(G))
  flag && return true, T(p)
  return false, nothing
end

is_pgroup_with_prime(G::GrpAbFinGen) = is_pgroup_with_prime(ZZRingElem, G)

# Let `v` be a vector of integers.
# This function returns the unique sorted vector `w` of zeros and prime powers
# such that `v` and `w` describe the same abelian group in the sense that
# the direct product of the groups `ZZ /(v[i]*ZZ)` is isomorphic to
# the direct product of the groups `ZZ /(w[i]*ZZ)`.
function abelian_invariants_of_vector(::Type{T}, v::Vector) where T <: IntegerUnion
  invs = T[]
  for elm in v
    if elm == 0
      push!(invs, 0)
    elseif 1 < elm
      append!(invs, [x[1]^x[2] for x in factor(elm)])
    elseif elm < -1
      append!(invs, [x[1]^x[2] for x in factor(-elm)])
    end
  end
  return sort!(invs)
end

# Let `v` be a vector of integers.
# This function returns the unique vector `w` of nonnegative integers
# such that `w[1] != 1`, `w[i]` divides `w[i+1]` for `1 < i < length(w)-1`
# and such that `v` and `w` describe the same abelian group in the sense that
# the direct product of the groups `ZZ /(v[i]*ZZ)` is isomorphic to
# the direct product of the groups `ZZ /(w[i]*ZZ)`.
function elementary_divisors_of_vector(::Type{T}, v::Vector) where T <: IntegerUnion
  invs = T[]
  d = Dict{T, Vector{T}}()
  for elm in v
    if elm == 0
      push!(invs, 0)
    else
      if elm < -1
        elm = -elm
      end
      for (p, e) in factor(elm)
        if haskey(d, p)
          push!(d[p], p^e)
        else
          d[p] = T[p^e]
        end
      end
    end
  end
  l = 0
  ps = keys(d)
  for p in ps
    l = max(l, length(d[p]))
  end
  o = T(1)
  for p in ps
    dp = d[p]
    for i in 1:(l-length(dp))
      push!(dp, o)
    end
    sort!(dp)
  end
  for i in l:-1:1
    e = o
    for p in ps
      e = e * d[p][i]
    end
    push!(invs, e)
  end
  return reverse(invs)
end

abelian_invariants(::Type{T}, G::GrpAbFinGen) where T <: IntegerUnion =
  abelian_invariants_of_vector(T, elementary_divisors(G))

abelian_invariants(G::GrpAbFinGen) = abelian_invariants(ZZRingElem, G)

function abelian_invariants_schur_multiplier(::Type{T}, G::GrpAbFinGen) where T <: IntegerUnion
  # By a theorem of I. Schur,
  # the multiplier of an abelian group with elementary divisors
  # n_1 | n_2 | ... | n_k, with k > 1,
  # has the elementary divisors n_i with multiplicity k-i, for 1 <= i < k.
  invs = elementary_divisors(G)
  res = T[]
  k = length(invs)
  for i in 1:(k-1)
    append!(res, repeat(T[invs[i]], k-i))
  end
  return abelian_invariants_of_vector(T, res)
end

abelian_invariants_schur_multiplier(G::GrpAbFinGen) = abelian_invariants_schur_multiplier(ZZRingElem, G)

nilpotency_class(G::GrpAbFinGen) = (order(G) == 1 ? 0 : 1)

# helper for prime_of_pgroup: this helper is efficient thanks to
# @gapattribute, but we also want prime_of_pgroup to accept an optional
# type argument; so we cannot use @gapattribute directly. Instead we set
# up this auxiliary _prime_of_pgroup which then is called by the real
# prime_of_pgroup.
# TODO: enhance @gapattribute so this is not necessary
function _prime_of_pgroup(G::GrpAbFinGen)
  flag, _, p = is_prime_power_with_data(order(G))
  @req flag "only supported for non-trivial p-groups"
  return p
end

function prime_of_pgroup(::Type{T}, G::GrpAbFinGen) where T <: IntegerUnion
  return T(_prime_of_pgroup(G))
end

prime_of_pgroup(G::GrpAbFinGen) = prime_of_pgroup(ZZRingElem, G)

#TODO There is no underlying GAP attribute, we do not store the prime.
#   has_prime_of_pgroup(G::GrpAbFinGen) = false
#   function set_prime_of_pgroup(G::GrpAbFinGen, p::IntegerUnion)
#     # do nothing
#   end
# The same holds for `fitting_subgroup`, `frattini_subgroup`, `socle`, ...
