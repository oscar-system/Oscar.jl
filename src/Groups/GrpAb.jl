# additional functions for abelian groups (type `FinGenAbGroup`)

# Restrict a morphism f : G -> H to the surjective morphism g : G -> im(G)
function restrict_codomain(f::FinGenAbGroupHom)
  G = domain(f)
  H, Htocd = image(f, false)
  imgs = elem_type(H)[]
  for g in gens(G)
    fl, h = has_preimage_with_preimage(Htocd, f(g))
    @assert fl
    push!(imgs, h)
  end
  return hom(G, H, imgs)
end


# Compute whether `x` is contained in `G`, modulo natural embeddings.
# (This is analogous to `issubset`.)
is_element(x::FinGenAbGroupElem, G::FinGenAbGroup) = issubset(sub([x])[1], G)

_coeff(x::FinGenAbGroupElem) = x.coeff

function is_finite_order(a::FinGenAbGroupElem)
  G, m = snf(a.parent)
  b = preimage(m, a)
  return any(i -> iszero(G.snf[i]) && !iszero(b[i]), 1:ngens(G))
end

function order(::Type{T}, a::FinGenAbGroupElem) where T <: IntegerUnion
   return T(order(a))
end

# order(x::FinGenAbGroupElem) = order(ZZRingElem, x) # see Hecke.jl/src/GrpAb/Elem.jl

function order(::Type{T}, G::FinGenAbGroup) where T <: IntegerUnion
   return T(order(G))
end

#TODO: do we need `one`? (Hecke defines `is_one(x::FinGenAbGroup)`.)
# one(x::FinGenAbGroup) == zero(x)
# one(x::FinGenAbGroupElem) == zero(x)
# Base.inv(x::FinGenAbGroupElem) = -x
# x^n for FinGenAbGroupElem -> n*x

Base.:^(g::FinGenAbGroupElem, x::FinGenAbGroupElem) = g

function comm(x::FinGenAbGroupElem, y::FinGenAbGroupElem)
  @req (x.parent == y.parent) "elements must belong to the same group"
  return zero(parent(x))
end

has_gens(G::FinGenAbGroup) = true

function small_generating_set(G::FinGenAbGroup)
   is_snf(G) && return gens(G)
   S, mp = snf(G)
   return [mp(x) for x in gens(S)]
end

function index(::Type{I}, G::T, H::T) where I <: IntegerUnion where T <: FinGenAbGroup
   @req is_subgroup(H, G)[1] "H must be a subgroup of G"
   f = count(is_zero, snf(G)[1].snf) - count(is_zero, snf(H)[1].snf)
   @req f == 0 "index is supported only for subgroups of finite index"
   return I(divexact(order(torsion_subgroup(G)[1]), order(torsion_subgroup(H)[1]), check = false))
end


################################################################################
#
#   Conjugacy Classes in FinGenAbGroup groups: just wrap elements
#
################################################################################

struct FinGenAbGroupConjClass{T<:FinGenAbGroup, S<:Union{FinGenAbGroupElem,FinGenAbGroup}} <: GroupConjClass{T, S}
   X::T
   repr::S
end

function Base.hash(C::FinGenAbGroupConjClass, h::UInt)
  return Base.hash(Representative(C), h)
end

function Base.show(io::IO, C::FinGenAbGroupConjClass)
  print(io, string(representative(C)), " ^ ", string(C.X))
end

==(a::FinGenAbGroupConjClass, b::FinGenAbGroupConjClass) = representative(a) == representative(b)

Base.length(::Type{T}, C::FinGenAbGroupConjClass) where T <: IntegerUnion = T(1)

Base.length(C::FinGenAbGroupConjClass) = ZZRingElem(1)


################################################################################
#
#   Conjugacy classes of elements
#
################################################################################

conjugacy_class(G::FinGenAbGroup, g::FinGenAbGroupElem) = FinGenAbGroupConjClass(G, g)

Base.eltype(::Type{FinGenAbGroupConjClass{T,S}}) where {T,S} = S

Base.rand(C::FinGenAbGroupConjClass) = representative(C)

Base.rand(rng::Random.AbstractRNG, C::FinGenAbGroupConjClass) = representative(C)

number_of_conjugacy_classes(G::FinGenAbGroup) = order(ZZRingElem, G)

number_of_conjugacy_classes(::Type{T}, G::FinGenAbGroup) where T <: IntegerUnion = order(T, G)

conjugacy_classes(G::FinGenAbGroup) = [FinGenAbGroupConjClass(G, x) for x in G]

is_conjugate(G::FinGenAbGroup, x::FinGenAbGroupElem, y::FinGenAbGroupElem) = (x == y)

function is_conjugate_with_data(G::FinGenAbGroup, x::FinGenAbGroupElem, y::FinGenAbGroupElem)
   x == y ? (true, zero(G)) : (false, nothing)
end

################################################################################
#
#   Conjugacy classes of subgroups
#
################################################################################

conjugacy_class(G::T, H::T) where T <: FinGenAbGroup = FinGenAbGroupConjClass(G, H)

function subgroup_classes(G::FinGenAbGroup; order::T = ZZRingElem(-1)) where T <: IntegerUnion
   @req is_finite(G) "G is not finite"
   if order > 0 && mod(Hecke.order(G), order) != 0
     # `subgroups` would throw an error
     return FinGenAbGroupConjClass{FinGenAbGroup, FinGenAbGroup}[]
   end
   return [conjugacy_class(G, H) for (H, mp) in Hecke.subgroups(G, order = order)]
end

function low_index_subgroup_classes(G::FinGenAbGroup, n::Int)
   @req (n > 0) "index must be positive"
   res = [conjugacy_class(G, G)]
   ord = order(G)
   for i in 2:n
     if mod(ord, i) == 0
       append!(res, [conjugacy_class(G, H) for (H, mp) in Hecke.subgroups(G, index = i)])
     end
   end
   return res
end

function maximal_subgroup_classes(G::FinGenAbGroup)
   @req is_finite(G) "G is not finite"
   primes = [p for (p, e) in factor(order(G))]
   res = typeof(G)[]
   for p in primes
     append!(res, [H for (H, mp) in Hecke.subgroups(G, index = p)])
   end
   return [conjugacy_class(G, H) for H in res]
end

conjugate_group(G::FinGenAbGroup, x::FinGenAbGroupElem) = G

Base.:^(G::FinGenAbGroup, x::FinGenAbGroupElem) = G

function is_conjugate(G::FinGenAbGroup, H::FinGenAbGroup, K::FinGenAbGroup)
   return is_subgroup(H, G)[1] && is_subgroup(K, G)[1] &&
          is_subgroup(H, K)[1] && is_subgroup(K, H)[1]
end

function is_conjugate_with_data(G::FinGenAbGroup, H::FinGenAbGroup, K::FinGenAbGroup)
   if is_subgroup(H, K)[1] && is_subgroup(K, H)[1]
     return true, zero(G)
   else
     return false, nothing
   end
end

is_conjugate_subgroup(G::T, U::T, V::T) where T <: FinGenAbGroup = is_subgroup(V, U)[1]
is_conjugate_subgroup_with_data(G::T, U::T, V::T) where T <: FinGenAbGroup = is_subgroup(V, U)[1], zero(G)

Base.IteratorSize(::Type{<:FinGenAbGroupConjClass}) = Base.HasLength()

Base.iterate(C::FinGenAbGroupConjClass) = iterate(C, 0)

function Base.iterate(C::FinGenAbGroupConjClass, state::Int)
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

function core(G::FinGenAbGroup, H::FinGenAbGroup)
  flag, emb = is_subgroup(H, G)
  @req flag "H  must be a subgroup of G"
  return (H, emb)
end

function normalizer(G::FinGenAbGroup, H::FinGenAbGroup)
  @req is_subgroup(H, G)[1] "H  must be a subgroup of G"
  return (G, identity_map(G))
end

function normalizer(G::FinGenAbGroup, x::FinGenAbGroupElem)
  @req is_element(x, G) "x does not lie in G"
  return (G, identity_map(G))
end

function normal_closure(G::FinGenAbGroup, H::FinGenAbGroup)
  flag, emb = is_subgroup(H, G)
  @req flag "H  must be a subgroup of G"
  return (H, emb)
end

pcore(G::FinGenAbGroup, p::IntegerUnion) = sylow_subgroup(G, p)


################################################################################
#
# Specific Subgroups
#
################################################################################

fitting_subgroup(G::FinGenAbGroup) = (G, identity_map(G))

function frattini_subgroup(G::FinGenAbGroup)
   @req is_finite(G) "G is not finite"
   subgens = FinGenAbGroupElem[]
   for x in gens(G)
     for (p, e) in factor(order(x))
       x = p*x
     end
     if !is_zero(x)
       push!(subgens, x)
     end
   end
   return sub(G, subgens)
end

function socle(G::FinGenAbGroup)
   @req is_finite(G) "G is not finite"
   subgens = FinGenAbGroupElem[]
   for x in gens(G)
     n = 1
     ord = order(x)
     for (p, e) in factor(ord)
       n = p*n
     end
     x = divexact(ord, n)*x
     if !is_zero(x)
       push!(subgens, x)
     end
   end
   return sub(G, subgens)
end

trivial_subgroup(G::FinGenAbGroup) = sub(G, FinGenAbGroupElem[])

solvable_radical(G::FinGenAbGroup) = (G, identity_map(G))


################################################################################
#
# Sylow & Hall Subgroups
#
################################################################################

#TODO: how to compute complements?
# complement_classes(G::T, N::T) where T <: FinGenAbGroup
# complement_system(G::FinGenAbGroup)

function sylow_system(G::FinGenAbGroup)
   @req is_finite(G) "G is not finite"
   result = FinGenAbGroup[]
   for (p, e) in factor(order(G))
     push!(result, sylow_subgroup(G, p)[1])
   end
   return result
end

function hall_subgroup_classes(G::FinGenAbGroup, P::AbstractVector{<:IntegerUnion})
   @req is_finite(G) "G is not finite"
   P = unique(P)
   @req all(is_prime, P) "The integers must be prime"
   subgens = FinGenAbGroupElem[]
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
   return [conjugacy_class(G, sub(G, subgens)[1])]
end

function hall_system(G::FinGenAbGroup)
   @req is_solvable(G) "G must be solvable"
   primes = [p for (p, e) in factor(order(G))]
   result = FinGenAbGroup[]
   for P in subsets(Set(primes))
     push!(result, representative(hall_subgroup_classes(G, collect(P))[1]))
   end
   return result
end


################################################################################
#
#   G-Sets of `FinGenAbGroup`s
#
################################################################################

function action_homomorphism(Omega::GSetByElements{FinGenAbGroup, S}) where S
  A = acting_group(Omega)

  # Compute a permutation group `G` isomorphic with `A`.
  phi = isomorphism(PermGroup, A)
  G = codomain(phi)

  # Let `G` act on `Omega` as `A` does.
  OmegaG = induce(Omega, inv(phi))

  # Compute the permutation action on `1:length(Omega)`
  # corresponding to the action of `A` on `Omega`.
  return compose(phi, action_homomorphism(OmegaG))
end


################################################################################
#
#   Properties
#
################################################################################

is_almost_simple(G::FinGenAbGroup) = false

is_finitely_generated(G::FinGenAbGroup) = true

is_perfect(G::FinGenAbGroup) = is_trivial(G)

is_pgroup(G::FinGenAbGroup) = is_pgroup_with_prime(G)[1]

is_quasisimple(G::FinGenAbGroup) = false

is_simple(G::FinGenAbGroup) = is_finite(G) && is_prime(order(G))

is_solvable(G::FinGenAbGroup) = true

is_sporadic_simple(G::FinGenAbGroup) = false

function is_pgroup_with_prime(::Type{T}, G::FinGenAbGroup) where T <: IntegerUnion
  is_trivial(G) && return true, nothing
  is_finite(G) || return false, nothing
  flag, _, p = is_prime_power_with_data(order(G))
  flag && return true, T(p)
  return false, nothing
end

is_pgroup_with_prime(G::FinGenAbGroup) = is_pgroup_with_prime(ZZRingElem, G)

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

abelian_invariants(::Type{T}, G::FinGenAbGroup) where T <: IntegerUnion =
  abelian_invariants_of_vector(T, elementary_divisors(G))

abelian_invariants(G::FinGenAbGroup) = abelian_invariants(ZZRingElem, G)

function abelian_invariants_schur_multiplier(::Type{T}, G::FinGenAbGroup) where T <: IntegerUnion
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

abelian_invariants_schur_multiplier(G::FinGenAbGroup) = abelian_invariants_schur_multiplier(ZZRingElem, G)

nilpotency_class(G::FinGenAbGroup) = (order(G) == 1 ? 0 : 1)

# helper for prime_of_pgroup: this helper is efficient thanks to
# @gapattribute, but we also want prime_of_pgroup to accept an optional
# type argument; so we cannot use @gapattribute directly. Instead we set
# up this auxiliary _prime_of_pgroup which then is called by the real
# prime_of_pgroup.
# TODO: enhance @gapattribute so this is not necessary
function _prime_of_pgroup(G::FinGenAbGroup)
  flag, _, p = is_prime_power_with_data(order(G))
  @req flag "only supported for non-trivial p-groups"
  return p
end

function prime_of_pgroup(::Type{T}, G::FinGenAbGroup) where T <: IntegerUnion
  return T(_prime_of_pgroup(G))
end

prime_of_pgroup(G::FinGenAbGroup) = prime_of_pgroup(ZZRingElem, G)

#TODO There is no underlying GAP attribute, we do not store the prime.
#   has_prime_of_pgroup(G::FinGenAbGroup) = false
#   function set_prime_of_pgroup(G::FinGenAbGroup, p::IntegerUnion)
#     # do nothing
#   end
# The same holds for `fitting_subgroup`, `frattini_subgroup`, `socle`, ...
