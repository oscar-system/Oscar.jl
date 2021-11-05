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

  factors = MSet{fmpz}()
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

  n = degree(group(R))

  # Find all possible tuples (d_1, ..., d_n), such that d_1 \cdots d_n = k*order(group(R))
  degrees = Vector{Vector{fmpz}}()
  for part in iterate_partitions(factors)
    # We actually only need the multiset partitions with at most n parts.
    if length(part) > n
      continue
    end
    ds = ones(fmpz, n)
    for i = 1:length(part)
      for (p, e) in part[i].dict
       ds[i + n - length(part)] *= p^e
      end
    end

    sort!(ds)

    if test_primary_degrees_via_hilbert_series(R, ds)
      push!(degrees, ds)
    end
  end

  sort!(degrees, lt = (x, y) -> sum(x) < sum(y) || x < y)
  return degrees
end

# Checks whether there exist homogeneous f_1, ..., f_k in RG of degrees
# degrees[length(invars) + 1:length(invars) + k] such that
# RG/< invars, f_1, ..., f_k > has Krull dimension n - k, where n == length(invars).
# If the base field is finite, the answer "true" might be wrong (for theoretical reasons).
# See Kem99, Theorem 2.
function check_primary_degrees(RG::InvRing{FldT, GrpT, PolyElemT}, degrees::Vector{fmpz}, invars::Vector{PolyElemT}, k::Int, iters::Dict{Int, <: VectorSpaceIterator}, ideals::Dict{Set{PolyElemT}, Int}) where {FldT, GrpT, PolyElemT}
  R = polynomial_ring(RG)
  n = length(degrees)

  deg_dict = Dict{fmpz, Int}()
  for e in degrees[length(invars) + 1:length(invars) + k]
    deg_dict[e] = get(deg_dict, e, 0) + 1
  end
  for degs in Hecke.subsets(Set(keys(deg_dict)))
    gens_ideal = Set(invars)
    for e in degs
      iter = get!(iters, Int(e)) do
        vector_space_iterator(coefficient_ring(RG), iterate_basis(RG, Int(e)))
      end
      union!(gens_ideal, collect_basis(iter))
    end
    numbersInds = !isempty(degs) ? sum(deg_dict[e] for e in degs) : 0
    if isempty(gens_ideal)
      push!(gens_ideal, R())
    end
    dimI = get!(ideals, gens_ideal) do
      dim(ideal(R, collect(gens_ideal)))
    end
    if dimI > n - length(invars) - numbersInds
      return false
    end
  end
  return true
end

function primary_invariants_via_optimal_hsop(RG::InvRing)
  iters = Dict{Int, VectorSpaceIterator}()
  ideals = Dict{Set{elem_type(polynomial_ring(RG))}, Int}()
  k = 1
  while true
    degrees = candidates_primary_degrees(RG, k)
    for ds in degrees
      invars = elem_type(polynomial_ring(RG))[]
      b, _ = primary_invariants_via_optimal_hsop!(RG, ds, invars, iters, ideals)
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
function primary_invariants_via_optimal_hsop!(RG::InvRing{FldT, GrpT, PolyElemT}, degrees::Vector{fmpz}, invars::Vector{PolyElemT}, iters::Dict{Int, <: VectorSpaceIterator}, ideals::Dict{Set{PolyElemT}, Int}, k::Int = 0) where {FldT, GrpT, PolyElemT}
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
        if !check_primary_degrees(RG, degrees, invars, k - 1, iters, ideals)
          pop!(invars)
          continue
        end
      end
      b, kk = primary_invariants_via_optimal_hsop!(RG, degrees, invars, iters, ideals, k - 1)
      if b
        return true, 0
      end
      pop!(invars)
      if kk >= k
        k = kk + 1
        if !check_primary_degrees(RG, degrees, invars, k, iters, ideals)
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
