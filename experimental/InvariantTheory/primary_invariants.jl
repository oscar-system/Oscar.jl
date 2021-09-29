mutable struct PrimaryInvarsCache{FldT, GrpT, PolyElemT, PolyRingT, ActionT, SingularActionT}
  R::InvRing{FldT, GrpT, PolyElemT, PolyRingT, ActionT, SingularActionT}
  bases::Dict{fmpz, Vector{PolyElemT}}

  function PrimaryInvarsCache(R::InvRing{FldT, GrpT, PolyElemT, PolyRingT, ActionT, SingularActionT}) where {FldT, GrpT, PolyElemT, PolyRingT, ActionT, SingularActionT}
    C = new{FldT, GrpT, PolyElemT, PolyRingT, ActionT, SingularActionT}()
    C.R = R
    C.bases = Dict{fmpz, Vector{PolyElemT}}()
    return C
  end
end

mutable struct VectorSpaceIterator{IteratorT}
  basis_iterator::IteratorT
  rand_bound::Int

  function VectorSpaceIterator(basis_iterator::IteratorT, bound::Int = 10^5) where {IteratorT}
    VSI = new{IteratorT}()
    VSI.basis_iterator = basis_iterator
    VSI.rand_bound = bound
    return VSI
  end
end

function Base.iterate(VSI::VectorSpaceIterator)
  if isempty(VSI.basis_iterator)
    return nothing
  end

  b, s = iterate(VSI.basis_iterator)
  basis_collected = typeof(b)[ b ]
  return b, (s, Int[ length(VSI.basis_iterator) ], basis_collected, 1)
end

# We "iterate" the vector space in 3 phases:
# 1) iterate the basis using VSI.basis_iterator, so return one basis element at
#    a time
# 2) iterate all possible sums of basis elements
# 3) return random linear combinations of the basis elements
#
# state is supposed to be a tuple of length 4 containing:
# - the state of VSI.basis_iterator
# - a Vector{Int} being the state for "phase 2"
# - a Vector containing the (so far) collected basis
# - an Int: the number of the phase
function Base.iterate(VSI::VectorSpaceIterator, state)
  phase = state[4]

  if phase == 1
    bs = iterate(VSI.basis_iterator, state[1])
    if bs == nothing
      phase = 2
    else
      b = bs[1]
      s = bs[2]
      push!(state[3], b)
      return b, (s, state[2], state[3], state[4])
    end
  end

  if phase == 2
    expand = true
    s = state[2]
    if s[end] < length(state[3])
      s[end] += 1
      expand = false
    else
      for i = length(s) - 1:-1:1
        if s[i] + 1 < s[i + 1]
          s[i] += 1
          for j = i + 1:length(s)
            s[j] = s[j - 1] + 1
          end
          expand = false
          break
        end
      end
    end

    if expand
      if length(s) < length(state[3])
        s = collect(1:length(s) + 1)
      else
        phase = 3
      end
    end
    if phase != 3
      return sum([ state[3][i] for i in s ]), (state[1], s, state[3], 2)
    end
  end

  coeffs = rand(-VSI.rand_bound:VSI.rand_bound, length(state[3]))
  return sum([ coeffs[i]*state[3][i] for i = 1:length(state[3]) ]), (state[1], state[2], state[3], 3)
end

function basis(C::PrimaryInvarsCache, d::fmpz)
  return get!(C.bases, d) do
    basis(C.R, d)
  end
end

# If d_1, ..., d_n are degrees of primary invariants, then the Hilbert series
# must be f(t)/\prod_i (1 - t^{d_i}) where f is a polynomial with non-negative
# integer coefficients.
# See Derksen, Kemper "Computational Invariant Theory", pp. 94, 95.
function test_primary_degrees_via_hilbert_series(R::InvRing, degrees::Vector{fmpz})
  @assert !ismodular(R)

  mol = molien_series(R)
  f = numerator(mol)
  g = denominator(mol)
  S = parent(f)
  t = gen(S)
  for d in degrees
    f *= (1 - t^d)
  end
  if !iszero(mod(f, g))
    return false
  end
  h = div(f, g)
  for c in coefficients(h)
    if !isinteger(c)
      return false
    end
    if !(c isa fmpq)
      c = coeff(c, 0)
    end
    if c < 0
      return false
    end
  end
  return true
end

# Returns possible degrees of primary invariants d_1, ..., d_n with
# d_1 \cdots d_n == k*|G|, where G = group(R).
# TODO: Possibly include checks involving earlier runs of the search, see Kemper, p. 181
function candidates_primary_degrees(R::InvRing, k::Int)
  # If we can't factor the group order in reasonable time, we might as well
  # give up.
  factors = Vector{fmpz}()
  for n in [ k, order(group(R)) ]
    fac = factor(n) # TODO: cache this!
    for (p, e) in fac
      for i = 1:e
        push!(factors, p)
      end
    end
  end

  n = degree(group(R))
  parts = Vector{Any}()
  for i = 1:n
    # Actually, we want to multiset partitions here (as the factors may show
    # up more than once). But GAP doesn't have this as far as I'm aware.
    # It's in JuLie.jl though.
    append!(parts, GAP.gap_to_julia(GAP.Globals.PartitionsSet(GAP.Globals.Set(GAP.julia_to_gap(1:length(factors))), i)))
  end

  #expG = exponent(group(R))

  sorted_degrees = Vector{Tuple{fmpz, Vector{fmpz}}}()
  for part in parts
    ds = ones(fmpz, n)
    for i = 1:length(part)
      for j in part[i]
       ds[i + n - length(part)] *= factors[j]
      end
    end
    #if !iszero(mod(lcm(ds), expG))
    #  continue
    #end
    sort!(ds) # To avoid duplicates
    s = sum(ds)
    sds = (s, ds)
    l = searchsortedfirst(sorted_degrees, sds)
    # If we used multiset partitions there shouldn't be duplicates (I think)
    if l <= length(sorted_degrees) && sorted_degrees[l] == sds
      continue
    end
    insert!(sorted_degrees, l, sds)
  end

  degrees = Vector{Vector{fmpz}}()
  for (s, ds) in sorted_degrees
    if test_primary_degrees_via_hilbert_series(R, ds)
      push!(degrees, ds)
    end
  end
  return degrees
end

function check_primary_degrees(RG::InvRing{FldT, GrpT, PolyElemT}, degrees::Vector{fmpz}, invars::Vector{PolyElemT}, k::Int, C::PrimaryInvarsCache) where {FldT, GrpT, PolyElemT}
  R = polynomial_ring(RG)
  n = length(degrees)

  deg_dict = Dict{fmpz, Int}()
  for e in degrees[length(invars) + 1:length(invars) + k]
    deg_dict[e] = get(deg_dict, e, 0) + 1
  end
  for degs in Hecke.subsets(Set(keys(deg_dict)))
    if isempty(invars)
      I = ideal(R, R())
    else
      I = ideal(R, invars)
    end
    for e in degs
      I = I + ideal(R, basis(C, e))
    end
    numbersInds = sum(Int[ deg_dict[e] for e in degs ])
    if dim(I) > n - length(invars) - numbersInds
      return false
    end
  end
  return true
end

function primary_invariants_via_optimal_hsop(RG::InvRing)
  C = PrimaryInvarsCache(RG)
  k = 1
  while true
    degrees = candidates_primary_degrees(RG, k)
    for ds in degrees
    #@info "Checking k = $k"
      invars = elem_type(polynomial_ring(RG))[]
      b, _ = primary_invariants_via_optimal_hsop!(RG, ds, invars, C)
      b && return invars
    end
    k += 1
  end
end

# Kemper "An Algorithm to Calculate Optimal Homogeneous Systems of Parameters", 1999
# Returns a bool b and an integer k.
# b == true iff primary invariants of the given degrees exist. In this case
# invars will contain those invariants.
# k is only needed for recursive calls of the function.
function primary_invariants_via_optimal_hsop!(RG::InvRing{FldT, GrpT, PolyElemT}, degrees::Vector{fmpz}, invars::Vector{PolyElemT}, C::PrimaryInvarsCache) where {FldT, GrpT, PolyElemT}
  k = 0
  n = length(degrees) - length(invars)
  R = polynomial_ring(RG)

  B = Vector{PolyElemT}()
  if length(invars) != length(degrees)
    d = degrees[length(invars) + 1]
    #@time B = basis(C, d)
    B = iterate_basis(RG, Int(d))
  end

  for f in VectorSpaceIterator(B)
    if iszero(f) || f in invars
      continue
    end
    push!(invars, f)
    if k > 0
      #@info "Test needed 1"
      if !check_primary_degrees(RG, degrees, invars, k - 1, C)
        pop!(invars)
        #@info "Test failed"
        continue
      end
      #@info "Test succeeded"
    end
    b, kk = primary_invariants_via_optimal_hsop!(RG, degrees, invars, C)
    if b
      return true, 0
    end
    pop!(invars)
    if kk >= k
      k = kk + 1
      #@info "Test needed 2"
      if !check_primary_degrees(RG, degrees, invars, k, C)
        #@info "Test failed"
        break
        #return false, 0
      end
      #@info "Test succeeded"
    end
  end

  if k <= 1
    #@info "Compute GB"
    if !isempty(invars) && dim(ideal(polynomial_ring(RG), invars)) > n
      k = 0
    else
      k = 1
    end
  end
  if iszero(n) && isone(k)
    return true, k
  end
  return false, k
end

# Shouldn't be here and should be like this
order(G::MatrixGroup{T}) where {T <: Union{nf_elem, fmpq}} = order(isomorphic_group_over_finite_field(G)[1])
exponent(G::MatrixGroup{T}) where {T <: Union{nf_elem, fmpq}} = exponent(isomorphic_group_over_finite_field(G)[1])

function symmetric_matrix_group(n::Int)
  @assert n >= 2
  M = identity_matrix(FlintQQ, n)
  matrices = Vector{typeof(M)}(undef, n - 1)
  for i = 1:n - 1
    M[i, i] = zero(FlintQQ)
    M[i + 1, i + 1] = zero(FlintQQ)
    M[i + 1, i] = one(FlintQQ)
    M[i, i + 1] = one(FlintQQ)
    matrices[i] = M
    M = identity_matrix(FlintQQ, n)
  end
  return matrix_group(matrices)
end
