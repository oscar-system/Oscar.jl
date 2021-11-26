# If d_1, ..., d_n are degrees of primary invariants, then the Hilbert series
# must be f(t)/\prod_i (1 - t^{d_i}) where f is a polynomial with non-negative
# integer coefficients (non-negativity only holds in the non-modular case!).
# See DK15, pp. 94, 95.
function test_primary_degrees_via_hilbert_series(R::InvRing, degrees::Vector{Int})
  @assert !ismodular(R)

  mol = molien_series(R)
  f = numerator(mol)
  g = denominator(mol)
  for d in degrees
    # multiply f by 1 - t^d
    f -= shift_left(f, d)
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
function candidates_primary_degrees(R::InvRing, k::Int, bad_prefixes::Vector{Vector{Int}} = Vector{Vector{Int}}())

  factors = MSet{Int}()
  for n in [ k, order(group(R)) ]
    # If we can't factor the group order in reasonable time, we might as well
    # give up.
    fac = factor(n)
    for (p, e) in fac
      for i = 1:e
        push!(factors, Int(p))
      end
    end
  end

  n = degree(group(R))

  # Find all possible tuples (d_1, ..., d_n), such that d_1 \cdots d_n = k*order(group(R))
  degrees = Vector{Vector{Int}}()
  for part in iterate_partitions(factors)
    # We actually only need the multiset partitions with at most n parts.
    if length(part) > n
      continue
    end
    ds = ones(Int, n)
    for i = 1:length(part)
      for (p, e) in part[i].dict
       ds[i + n - length(part)] *= p^e
      end
    end

    sort!(ds)

    skip = false
    for prefix in bad_prefixes
      if ds[1:length(prefix)] == prefix
        skip = true
        break
      end
    end
    skip ? continue : nothing

    for d in ds
      # Check whether there exist invariants of this degree.
      if dimension_via_molien_series(fmpz, R, d) == 0
        skip = true
        break
      end
    end
    skip ? continue : nothing

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
function check_primary_degrees(RG::InvRing{FldT, GrpT, PolyElemT}, degrees::Vector{Int}, invars::Vector{PolyElemT}, k::Int, iters::Dict{Int, <: VectorSpaceIterator}, ideals::Dict{Set{PolyElemT}, Int}) where {FldT, GrpT, PolyElemT}
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
    dimI = get!(ideals, gens_ideal) do
      dim(ideal(R, collect(gens_ideal)))
    end
    if dimI > n - length(invars) - numbersInds
      return false
    end
  end
  return true
end

function primary_invariants_via_optimal_hsop(RG::InvRing, ensure_minimality::Int = 0)
  iters = Dict{Int, VectorSpaceIterator}()
  ideals = Dict{Set{elem_type(polynomial_ring(RG))}, Int}()
  bad_prefixes = Vector{Vector{Int}}()
  l = 1
  while true
    degrees = candidates_primary_degrees(RG, l, bad_prefixes)
    for ds in degrees
      invars = elem_type(polynomial_ring(RG))[]
      b, k = primary_invariants_via_optimal_hsop!(RG, ds, invars, iters, ideals, ensure_minimality, 0)
      b && return invars
      push!(bad_prefixes, ds[1:k])
    end
    l += 1
  end
end

# Kemper "An Algorithm to Calculate Optimal Homogeneous Systems of Parameters", 1999, [Kem99]
# Returns a bool b and an integer k.
# b == true iff primary invariants of the given degrees exist. In this case
# invars will contain those invariants.
# k is only needed for recursive calls of the function.
function primary_invariants_via_optimal_hsop!(RG::InvRing{FldT, GrpT, PolyElemT}, degrees::Vector{Int}, invars::Vector{PolyElemT}, iters::Dict{Int, <: VectorSpaceIterator}, ideals::Dict{Set{PolyElemT}, Int}, ensure_minimality::Int = 0, k::Int = 0) where {FldT, GrpT, PolyElemT}

  n = length(degrees) - length(invars)
  R = polynomial_ring(RG)

  if length(invars) != length(degrees)
    d = degrees[length(invars) + 1]

    iter = get!(iters, Int(d)) do
      vector_space_iterator(coefficient_ring(RG), iterate_basis(RG, Int(d)))
    end

    attempts = 1
    for f in iter
      if f in invars
        continue
      end
      if ensure_minimality > 0 && attempts > ensure_minimality
        break
      end
      attempts += 1
      push!(invars, f)
      if k > 0
        if !check_primary_degrees(RG, degrees, invars, k - 1, iters, ideals)
          pop!(invars)
          continue
        end
      end
      b, kk = primary_invariants_via_optimal_hsop!(RG, degrees, invars, iters, ideals, ensure_minimality, k - 1)
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
    k = 1
    dimI = get!(ideals, Set(invars)) do
      dim(ideal(polynomial_ring(RG), invars))
    end
    if dimI > n
      k = 0
    end
  end
  if iszero(n) && isone(k)
    return true, k
  end
  return false, k
end
