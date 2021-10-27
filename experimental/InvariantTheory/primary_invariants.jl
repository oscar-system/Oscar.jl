# If d_1, ..., d_n are degrees of primary invariants, then the Hilbert series
# must be f(t)/\prod_i (1 - t^{d_i}) where f is a polynomial with non-negative
# integer coefficients (non-negativity only holds in the non-modular case!).
# See DK15, pp. 94, 95.
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

  factors = Vector{fmpz}()
  for n in [ k, order(group(R)) ]
    # If we can't factor the group order in reasonable time, we might as well
    # give up.
    fac = factor(n) # TODO: cache this!
    for (p, e) in fac
      for i = 1:e
        push!(factors, p)
      end
    end
  end

  # Find all possible tuples (d_1, ..., d_n), such that d_1 \cdots d_n = k*order(group(R))
  n = degree(group(R))
  parts = Vector{Any}()
  for i = 1:n
    # Actually, we want multiset partitions here (as the factors may show
    # up more than once). But GAP doesn't have this as far as I'm aware.
    append!(parts, GAP.gap_to_julia(GAP.Globals.PartitionsSet(GAP.Globals.Set(GAP.julia_to_gap(1:length(factors))), i)))
  end

  sorted_degrees = Vector{Tuple{fmpz, Vector{fmpz}}}()
  for part in parts
    ds = ones(fmpz, n)
    for i = 1:length(part)
      for j in part[i]
       ds[i + n - length(part)] *= factors[j]
      end
    end
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

# Checks whether there exist homogeneous f_1, ..., f_k in RG of degrees
# degrees[length(invars) + 1:length(invars) + k] such that
# RG/< invars, f_1, ..., f_k > has Krull dimension n - k, where n == length(invars).
# If the base field is finite, the answer "true" might be wrong (for theoretical reasons).
# See Kem99, Theorem 2.
function check_primary_degrees(RG::InvRing{FldT, GrpT, PolyElemT}, degrees::Vector{fmpz}, invars::Vector{PolyElemT}, k::Int, iters::Dict{Int, <: VectorSpaceIterator}) where {FldT, GrpT, PolyElemT}
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
      iter = get!(iters, Int(e)) do
        vector_space_iterator(coefficient_ring(RG), iterate_basis(RG, Int(e)))
      end
      I = I + ideal(R, collect_basis(iter))
    end
    numbersInds = sum(Int[ deg_dict[e] for e in degs ])
    if dim(I) > n - length(invars) - numbersInds
      return false
    end
  end
  return true
end

function primary_invariants_via_optimal_hsop(RG::InvRing)
  iters = Dict{Int, VectorSpaceIterator}()
  k = 1
  while true
    degrees = candidates_primary_degrees(RG, k)
    for ds in degrees
      invars = elem_type(polynomial_ring(RG))[]
      b, _ = primary_invariants_via_optimal_hsop!(RG, ds, invars, iters)
      b && return invars
    end
    k += 1
  end
end

# Kemper "An Algorithm to Calculate Optimal Homogeneous Systems of Parameters", 1999, [Kem99]
# Returns a bool b and an integer k.
# b == true iff primary invariants of the given degrees exist. In this case
# invars will contain those invariants.
# k is only needed for recursive calls of the function.
function primary_invariants_via_optimal_hsop!(RG::InvRing{FldT, GrpT, PolyElemT}, degrees::Vector{fmpz}, invars::Vector{PolyElemT}, iters::Dict{Int, <: VectorSpaceIterator}) where {FldT, GrpT, PolyElemT}
  k = 0
  n = length(degrees) - length(invars)
  R = polynomial_ring(RG)

  if length(invars) != length(degrees)
    d = degrees[length(invars) + 1]

    iter = get!(iters, Int(d)) do
      vector_space_iterator(coefficient_ring(RG), iterate_basis(RG, Int(d)))
    end

    for f in iter
      if f in invars
        continue
      end
      push!(invars, f)
      if k > 0
        if !check_primary_degrees(RG, degrees, invars, k - 1, iters)
          pop!(invars)
          continue
        end
      end
      b, kk = primary_invariants_via_optimal_hsop!(RG, degrees, invars, iters)
      if b
        return true, 0
      end
      pop!(invars)
      if kk >= k
        k = kk + 1
        if !check_primary_degrees(RG, degrees, invars, k, iters)
          break
        end
      end
    end
  end

  if k <= 1
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

# Shouldn't be here and shouldn't be like this
order(G::MatrixGroup{T}) where {T <: Union{nf_elem, fmpq}} = order(isomorphic_group_over_finite_field(G)[1])
