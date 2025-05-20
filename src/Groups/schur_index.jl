# AbstractAlgebra.set_verbosity_level(:SchurIndices, 1)

"""
    yamada_example(p::Int, d::Int)

Return `G, chi` where `G` is a pc group of order `p*(p-1)*d`
and `chi` is a faithful irreducible character of degree `p-1` of `G`
that has Schur index `d`.

`p` must be a prime integer, and `d` must divide `p-1`.

This example is described in [Yam78](@cite).

# Examples
```jldoctest
julia> G, chi = Oscar.yamada_example(7, 3);

julia> order(G)
126

julia> schur_index(chi)
3

julia> local_schur_indices(chi)
1-element Vector{Pair{Int64, Int64}}:
 7 => 3
```
"""
function yamada_example(p::Int, d::Int)
  @req is_odd(p) && is_prime(p) "p must be an odd prime integer"
  @req mod(p-1, d) == 0 "d must divide p-1"

  r = Hecke.gen_mod_pk(ZZ(p))
  F = free_group([:x, :y])
  x, y = gens(F)
  rels = [x^p, y^(d*(p-1)), x*y/(y*x^r)] # the relators in the paper are wrong
  G = quo(F, rels)[1]
  P = codomain(isomorphism(PcGroup, G))

  # take any faithful character of degree p-1
  tbl = character_table(P)
  i = findfirst(x -> degree(x) == p-1 &&
                     length(class_positions_of_kernel(x)) == 1, tbl)

  return P, tbl[i]
end 

# See section 5 in [Ung19]
# return ("Q8", 1) or ("Q8", r) or ("QD", r) or ("", 0),
# where the last means that the given group is *not* one of the
# types in the Riese/Schmid paper
#
function _Riese_Schmid_type(G::GAPGroup)
  no_type = ("", 0)  # the group is not among the interesting ones
  is_finite(G) || return no_type
  n = order(G)

  # recognize the quaternion group of order 8
  id_q8 = (8, 4)
  n == 8 && small_group_identification(G) == id_q8 && return ("Q8", 1)

  # check for the structure U:P where U is a normal subgroup of odd prime
  # order r and P is a nontrivial 2-group
  is_even(n) || return no_type
  _, r = remove(n, 2)
  is_prime(r) || return no_type
  U = sylow_subgroup(G, r)[1]
  is_normalized_by(U, G) || return no_type

  # check for type ("Q8", r), i.e.,
  # X = C_P(U) is isomorphic with Q8,
  # P/X is nontrivial and cyclic,
  # P is a central product of X and C_P(X),
  # and |P/X| is the 2-part of the mult. order of 2 modulo r.
  P = sylow_subgroup(G, 2)[1]
  X = centralizer(P, U)[1]
  x = order(X)
  mord = modord(ZZRingElem(2), r)
  twopart, mord = ppio(mord, 2)
  p = order(P)
  if x == 8 && small_group_identification(X) == id_q8 &&
     p != x && x * twopart == p
    Q = quo(P, X)[1]
    if is_cyclic(Q)
      cpx = centralizer(P, X)[1]
      inter = intersect(X, cpx)
      C = sub(P, vcat(gens(X), gens(cpx)))[1]
      is_subset(inter[1], center(P)[1]) && order(C) == p && return ("Q8", r)
    end
  end

  # check for type ("QD", r), i.e.,
  # Z = P' is cyclic of order at least 4,
  # Y = C_P(z) satisfies Y/Z cyclic, where z is an element of order 4 in Z,
  # X = C_P(U) is not abelian,
  # Z has index 2 in X,
  # and the order of 2 modulo r has the form 2^s f where f is a (not nec. odd)
  # integer and |P/X| = 2^s.
  # (We exclude groups where f is odd and X is isomorphic with Q8;
  # these have been detected above with result ("Q8", r),
  # hence we need not check this property.)
  Z = derived_subgroup(P)[1]
  (is_cyclic(Z) && order(Z) >= 4) || return no_type
  z = cyclic_generator(Z)^div(order(Z), 4)
  Y = centralizer(P, z)[1]
  is_cyclic(quo(Y, Z)[1]) || return no_type
  is_abelian(X) && return no_type
  (order(X) == 2*order(Z)) || return no_type
  (mod(twopart, div(order(P), order(X))) == 0) || return no_type
  return ("QD", r)
end


#TODO: to be replaced in Hecke or so
"""
    prime_residues(n::T) where T <: IntegerUnion

Return the `Vector{T}` of all `i` in the range `0:(abs(n)-1)`
that are coprime to `n`.

# Examples
```jldoctest
julia> println(Oscar.prime_residues(20))
[1, 3, 7, 9, 11, 13, 17, 19]

julia> println(Oscar.prime_residues(0))
Int64[]

julia> println(Oscar.prime_residues(1))
[0]
```
"""
function prime_residues(n::T) where T <: IntegerUnion
  return Vector{T}(GAP.Globals.PrimeResidues(GapObj(n)))
end


#############################################################################
#
# deal with extensions of p-adic fields by elements of abelian number fields
#
# For `vals` with conductor N = p^k m, with m coprime to p,
# we describe Q_p(vals) by (N, stab, [Q_p(N):Q_p]),
# where stab is the set of prime residues u modulo N
# such that *u is in Gal(Q_p(N)/Q_p) and chi is fixed by *u
#
# We use that Gal(Q_p(N)/Q_p) is the direct product of the group generated
# by \sigma_p: \zeta_m \mapsto \zeta_m^p, \zeta_{p^k} \mapsto \zeta_{p^k}
# and the group of maps
# \sigma_u: \zeta_m \mapsto \zeta_m, \zeta_{p^k} \mapsto \zeta_{p^k}^u,
# where u is a prime residue modulo p^k.
#
# This means that [Q_p(vals):Q_p] == [Q_p(N):Q_p] / length(stab)
#
const _Q_p_chi_cache = Dict{Tuple{Int, Vector{Int}, Int},
                            Tuple{Int, Vector{Int}, Int}}()

function _Q_p_chi(vals::GapObj, p::ZZRingElem)
  if length(vals) == 0
    F = GAP.Globals.Rationals
  else
    F = GAP.Globals.Field(GAP.Globals.Rationals, vals)
  end
  N = GAP.Globals.Conductor(F)
  if N == 1
    stab = [1]
  else
    stab = Vector{Int}(GAP.Globals.GaloisStabilizer(F))
  end

  return get!(_Q_p_chi_cache, (N, stab, p)) do
    N == 1 && return (N, stab, 1)

    ppart, m = ppio(ZZRingElem(N), p)
    _, a, b = gcdx(ppart, m)
    r = m == 1 ? 1 : modord(p, m)
#T `modord(N, 1)` could return 1 (smallest pos. int. i s.t. 1 divides (N^i-1))
    stab = Int[]
    sigma_p = mod(1 + a*(p-1)*ppart, N)
    res = ppart == 1 ? [1] : prime_residues(ppart)
    for u in res
      for i in 0:(r-1)
        sigma = mod(sigma_p^i * (1 + (u-1)*m*b), N)
        if GAP.Globals.GaloisCyc(vals, GapObj(sigma)) == vals
          push!(stab, sigma)
        end
      end
    end

    return (N, stab, r*euler_phi(ppart))
  end
end


# `elm` is an element of a cyclotomic field in GAP.
# `info` is a description `(N, stab, index)` of the field
# $F \subseteq \Q_p(N)$ such that all elements in $F$ are fixed under the
# Galois automorphisms $*k$, for $k \in stab$.
# (We have $[\Q_p(N):F] = index$, but this is nnot used here.)
function _membership_test(elm::GAP.Obj, info::Tuple{Int, Vector{Int}, Int})
  return mod(info[1], GAPWrap.Conductor(elm)) == 0 &&
         all(k -> GAP.Globals.GaloisCyc(elm, k) == elm, info[2])
end


#############################################################################

"""
    local_schur_indices(chi::GAPGroupClassFunction; cyclic_defect::Vector{Int} = Int[])

Return the array of pairs `p => m_p(chi)` such that `m_p(chi) > 1` holds.
The Schur index over the real field is stored at `p = 0` if applicable.

`cyclic_defect` contains primes `p` such that `chi` is known to belong to
a `p`-block of cyclic defect.

The implementation follows [Ung19](@cite).

# Examples
```jldoctest
julia> G, chi = Oscar.yamada_example(5, 2);

julia> order(G)
40

julia> local_schur_indices(chi)
2-element Vector{Pair{Int64, Int64}}:
 0 => 2
 5 => 2

julia> schur_index(chi)
2
```
"""
function local_schur_indices(chi::GAPGroupClassFunction; cyclic_defect::Vector{Int} = Int[])
  @assert characteristic(chi) == 0
  tbl = parent(chi)
  deg = degree(ZZRingElem, chi)
  chipos = findfirst(isequal(chi), tbl)
  @assert chipos !== nothing "chi must be irreducible"
  @vprintln :SchurIndices 1 "called with character of degree $deg"
  result = Pair{Int, Int}[]

  # Step 1:
  # Let `m` denote the Schur index.
  # Determine m_\R(chi) and use Theorem 2.2 to find a bound u \in \N
  # such that m(chi) divides u.
  # If u = 1 terminate.
  # `m` divides the `chi(1)`.
  ind = indicator(chi)
  if ind == 0
    m_R = 1
  elseif ind == -1
    m_R = 2
    push!(result, 0 => 2)
  else
    m_R = 1
  end

  u = deg
  # `m` divides `|G|/chi(1)`.
  Gorder = order(tbl)
  u = gcd(u, div(Gorder, deg))
  if u == 1
    @vprintln :SchurIndices 1 "m = 1 (divides chi(1) and |G|/chi(1))"
    return result
  end

  # The character field contains a primitive `m`-th root of unity.
  if ind == 0
    # Compute the conductor `N` of the smallest cyclotomic field
    # that contains the character field of `chi`.
    gapfield = GAPWrap.Field(GapObj(chi))
    N = GAPWrap.Conductor(gapfield)
    # Compute the conductor of the largest cyclotomic field
    # that is contained in the character field of `chi`.
    for n in reverse(sort(divisors(N)))
#TODO: better compute with p-parts of N and form the lcm
      if GAPWrap.E(n) in gapfield
        u = gcd(u, Base.isodd(n) ? 2*n : n )
        cond = n
        break
      end
    end
  else
    # The character field is real, hence contains exactly 2nd roots of unity
    cond = 2
    u = gcd(u, 2)
  end

  # `m` divides the exponent of G.
  orders = orders_class_representatives(tbl)
  u = gcd(u, lcm(orders))
  if u == 1
    @vprintln :SchurIndices 1 "m = 1 (divides exp(G))"
    return result
  end

  # `m` divides \phi(2\prod_i p_i) where the product is taken over the
  # distinct prime divisors p_i of |G|.
  primes = sort(prime_divisors(Gorder))
  u = gcd(u, euler_phi(2*prod(primes)))
  if u == 1
    @vprintln :SchurIndices 1 "m = 1 (divides phi(2 prod_i p_i))"
    return result
  end

  # If `G` is a `q`-group and either `q` is odd or the char. field contains
  # a primitive 4-th root of unity then `m = 1`.
  if length(primes) == 1 &&
     ((primes[1] != 2) || (mod(cond, 4) == 0))
    u = 1
    @vprintln :SchurIndices 1 "m = 1 (q-group))"
    return result
  end

  # `m` divides the multiplicity of `chi` in any character that is
  # afforded by a rational representation.
  # (This is Thm. 2.1 in the 1983 paper by Feit.)
  # Here we consider permutation characters.
  for psi in induced_cyclic(tbl)
    u = gcd(u, scalar_product(ZZRingElem, chi, psi))
    if u == 1
      @vprintln :SchurIndices 1 "m = 1 (mult. of chi in cyclic perm. char.)"
      return result
    end
  end
  for name in names_of_fusion_sources(tbl)
    s = character_table(name)
    if s !== nothing
      known, fus = known_class_fusion(s, tbl)
      if known && length(class_positions_of_kernel(fus)) == 1
        psi = trivial_character(s)^(tbl)
        u = gcd(u, scalar_product(ZZRingElem, chi, psi))
        if u == 1
          @vprintln :SchurIndices 1 "m = 1 (mult. of chi in perm. char. from $name)"
          return result
        end
      end
    end
  end

  # If `q` is a prime such that the Sylow `q`-subgroups of G are
  # elementary abelian then `q` does not divide `m`.
  # (Here we need the group for the first time.)
  @req isdefined(tbl, :group) "cannot determine the Schur index with the currently used criteria"
  G = group(tbl)
  for q in primes
    if mod(Gorder, q^2) != 0 ||
       (mod(Gorder, q^3) != 0 && !(q in orders)) ||
       is_elementary_abelian(sylow_subgroup(G, q)[1])
      _, u = ppio(u, q)
    end
  end
  if u == 1
    @vprintln :SchurIndices 1 "m = 1 (el. ab. Sylow subgroups))"
    return result
  end

  # A prime power q > 2 can divide `m` only if there is a prime `p`
  # such that q divides p-1 and G contains an element of order p q.
  for (l, e) in factor(u)
    q = l^e
    if q > 2
      test = true
      while test
        test = false
        if all(p -> mod(p-1, q) != 0 || !(p*q in orders), primes)
          q = div(q, l)
          u = div(u, l)
          if q > 2
            test == true
          end
        end
      end
    end
  end
  if u == 1
    @vprintln :SchurIndices 1 "m = 1 (prime powers q))"
    return result
  end
  @vprintln :SchurIndices 1 "after step 1: bound is $u"

  # Step 2.
  # For each rational prime p dividing the order of G
  # use Theorems 4.4, 4.5, 4.6 and 4.7 to find u_p dividing u
  # such that m_p(chi) divides u_p .
  # Set u to be the LCM of m_\R(chi) and the u_p’s.
  # If u = 1 terminate.
  # (We store for each p the pair (u_p, exact) where u_p is a known multiple
  # of m_p(chi), and exact is `true` or `false`,
  # where `true` means that u_p = m_p(chi).)
  u_p = Dict([p => p == 2 ? (2, false) : (gcd(u, p-1), false) for p in primes])
  # Theorem 4.4
  if primes[1] == 2
    if mod(Gorder, 4) != 0 ||
       (mod(Gorder, 8) != 0 && !(4 in orders)) ||
       is_elementary_abelian(sylow_subgroup(G, 2)[1])
      u_p[2] = (1, true)
    else
      u_p[2] = (gcd(u_p[2][1], 2), false)
    end
  end
  @vprintln :SchurIndices 1 "initialize u_p: $u_p"

  # Theorem 4.5
  # 1. ...
#TODO: Which test characters should we use here?
#      (the natural character? other characters with known triv. m?)
  # 2. If `chi` has `p`-defect zero then `m_p(chi) = 1`
  # 3. If x is a `p`-regular element of `G` such that \psi(x) \in \Q_p(chi)
  #    for all \psi in Irr(B).
  #    Then `m_p(chi)` divides `chi(x)` (in the ring of alg. integers).
  GAPfield = GAPWrap.Field(GapObj(chi))  # Q(chi)
  GAPfield_p = Dict{Int, Tuple{Int, Vector, Int}}() # Q_p(chi), for the p that occur
  ords = orders_class_representatives(tbl)
  irr = map(GapObj, collect(tbl))
  for p in primes
    @vprintln :SchurIndices 1 "step 2, p = $p"
    pbl = block_distribution(tbl, p)
    blockpos = pbl[:block][chipos]
    if mod(div(Gorder, deg), p) != 0
      u_p[p] = (1, true)
    else
      block = findall(isequal(blockpos), pbl[:block])
      # description of \Q_p(chi)
      GAPfield_p[p] = _Q_p_chi(GapObj(chi), p)
      # Run over the p-regular elements (without the identity element).
      pregs = findall(x -> mod(x, p) != 0 && x != 1, ords)
      # we can ignore integer values in `chi`
      vals = filter(x -> isa(x, GapObj), GAP.Obj[x for x in GapObj(chi)])
      for preg in pregs
        # Run over the other characters in the block.
        if all(i -> _membership_test(irr[i][preg], GAPfield_p[p]), block)
          # m_p(chi) divides `chi[preg]` in the ring of algebraic integers;
          # the field element is represented w.r.t. an integral basis,
          # thus it is sufficient to look at its coefficients
          g = GAP.Globals.Gcd(GAP.Globals.COEFFS_CYC(GapObj(chi)[preg]))
          u_p[p] = (gcd(u_p[p][1], g), false)
        else
          # we can ignore integer values
          append!(vals, filter(x -> isa(x, GapObj), [irr[i][preg] for i in block]))
        end
      end

      # Theorem 4.6
      # If the p-block of `chi` has cyclic defect then
      # m_p(chi) = [K(chi):\Q_p(chi)].
      if p in cyclic_defect || is_block_with_cyclic_defect_group(tbl, p, blockpos)
        # K(\chi) describes the field over Q_p that is generated by
        # all values of `chi` and all values of
        # the restrictions of the ordinary irreducibles in the block of `chi`
        # to `p`-regular classes.
        K = _Q_p_chi(GapObj(vals), p)
        d_K = div(K[3], length(K[2])) # [K(\chi):\Q_p]
        d_chi = div(GAPfield_p[p][3], length(GAPfield_p[p][2])) # [\Q_p(chi):\Q_p]
        @assert mod(d_K, d_chi) == 0 "p = $p, d_K = $d_K, d_chi = $d_chi"
        d = div(d_K, d_chi)           # [K(\chi):\Q_p(\chi)]
        @assert mod(u_p[p][1], d) == 0
        u_p[p] = (d, true)
      end
#TODO: use more known results on cyclic defect cases (not stated in the paper)
#TODO: how to apply Theorem 4.7?
    end
    @vprintln :SchurIndices 1 "step 2, bound for p = $p is $(u_p[p])"
  end

  for pair in u_p
    if pair.second[1] == 1
      u_p[pair.first] = (1, true)
    end
  end
  if all(p -> u_p[p][2], primes)
    for p in primes
      pair = u_p[p]
      pair[1] > 1 && push!(result, p => pair[1])
    end
    @vprintln :SchurIndices 1 "after step 2, return $result"
    return result
  end

  # Step 3.
  # For each prime p with u_p > 1, and each prime q with q
  # dividing u_p, do the following steps.
  #   3a) Locate a quasi-elementary subgroup H of G and \eta \in Irr(H) such
  #       that q does not divide [\eta^G, \chi] ∣\Q_p(\chi, \eta):Q_p(\chi)∣;
  #       if such H or \eta does not exist then the q-part of u_p is 1.
  #   3b) Determine m_p(\eta).
  #   3c) Divide u_p by a suitable power of q so that the q-part of u_p
  #       becomes equal to
  #       m_{Q_p(\chi)}(\eta) =
  #        \frac{m_p(\eta}{\gcd( m_p(\eta), |Q_p(\chi, \eta):Q_p(\eta)∣)}.
  ccl = conjugacy_classes(G) # same ordering as columns of `tbl`
  for p in reverse(filter(p -> u_p[p][1] > 1, primes))
    @vprintln :SchurIndices 1 "step 3, p = $p"
    u_p[p][2] && continue
    for (q, e) in factor(u_p[p][1])
      # look at q-quasielementary subgroups
      @vprintln :SchurIndices 1 "step 3, q = $q"
      qregs = findall(x -> mod(x, q) != 0, ords)
      found = false
      for i in reverse(qregs)
        C = sub(G, [representative(ccl[i])])[1]
        N = normalizer(G, C)[1]
        Q = sylow_subgroup(N, q)[1]
        order(Q) > 1 || continue
        # H is a maximally q-quasielementary subgroup
        # with cyclic normal subgroup `C`.
        H = sub(G, vcat(gens(C), gens(Q)))[1]
        @vprintln :SchurIndices 1 "step 3, choose H of order $(order(H))"

        for eta in character_table(H)
          # we look ony at \eta with [\eta^G, \chi] not divisible by q
          mod(scalar_product(chi, eta^G), q) == 0 && continue

          # the field $\Q_p(\chi, \eta)$; again, ignore integer values
          vals = Vector{GapObj}(filter(x -> isa(x, GapObj), [x for x in GapObj(chi)]))
          append!(vals, filter(x -> isa(x, GapObj), [x for x in GapObj(eta)]))
          GAPfield_chieta = _Q_p_chi(GapObj(vals), p)
          @assert mod(GAPfield_chieta[3], GAPfield_p[p][3]) == 0
          mod(div(GAPfield_chieta[3], GAPfield_p[p][3]), q) == 0 && continue
          found = true
          @vprintln :SchurIndices 1 "step 3, choose eta of degree $(degree(eta))"

          # compute `m_p(eta)`, a power of `q`
          if p != q
            # Thm. 4.6 gives the answer.
            pairs = local_schur_indices(eta; cyclic_defect = [Int(p)])
            pairs = filter(x -> x[1] == p, pairs)
            if length(pairs) == 0
              @vprintln :SchurIndices 1 "step 3, m_$p(eta) is trivial (p != q)"
              m_p_eta = 1
            else
              m_p_eta = pairs[1][2][1]
              @vprintln :SchurIndices 1 "step 3, m_$p(eta) = $m_p_eta (p != q)"
            end
          elseif is_odd(p)
            # Use Thm. 4.4
            @vprintln :SchurIndices 1 "step 3, m_$p(eta) = 1 (p = q odd)"
            m_p_eta = 1
          else
            # p = q = 2
            type, r = _Riese_Schmid_type(H)
            m_p_eta = (type == "") ? 1 : 2
            @vprintln :SchurIndices 1 "step 3, m_$p(eta) = $m_p_eta (p = q = 2)"
          end

          # Now we derive the q-part of m_p(chi)
          GAPfield_eta = _Q_p_chi(GapObj([x for x in GapObj(eta)]), p)
          denom = gcd(m_p_eta, div(GAPfield_chieta[3], GAPfield_eta[3]))
          res = div(m_p_eta, denom)
          qpart, qprime = ppio(u_p[p][1], q)
          u_p[p] = (qprime * res, false)
          @vprintln :SchurIndices 1 "step 3, $q-part of m_$p(chi) is $res"
          found && break
        end
        found && break
      end
      if !found
        # We have no maximally q-quasielementary subgroup with
        # a character eta as in the paper,
        # thus we know that the q-part of m_p(chi) is 1.
        @vprintln :SchurIndices 1 "step 3, $q-part of m_$p(chi) is trivial"
        qpart, qprime = ppio(u_p[p][1], q)
        u_p[p] = (qprime, false)
      end
    end
    @vprintln :SchurIndices 1 "step 3, m_$p(chi) is $(u_p[p])"
  end

  # Step 4.
  # At this point each u_p = m_p(\chi).
  # Return m_\R(\chi) and the sequence of pairs `p => u_p` where u_p \not= 1.
  for p in primes
    pair = u_p[p]
    pair[1] > 1 && push!(result, p => pair[1])
  end
  @vprintln :SchurIndices 1 "step 4, return $result"
  return result
end
