export primary_invariants

# If d_1, ..., d_n are degrees of primary invariants, then the Hilbert series
# must be f(t)/\prod_i (1 - t^{d_i}) where f is a polynomial with integer
# coefficients respectively non-negative integer coefficients in the non-modular
# case.
# See DK15, pp. 94, 95.
function test_primary_degrees_via_hilbert_series(R::InvRing, degrees::Vector{Int})
  fl, h = _reduce_hilbert_series_by_primary_degrees(R, degrees)
  if !fl
    return false
  end
  for c in coefficients(h)
    if !isinteger(c)
      return false
    end
    if !is_modular(R) && c < 0
      return false
    end
  end
  return true
end

function reduce_hilbert_series_by_primary_degrees(R::InvRing, chi::Union{GAPGroupClassFunction, Nothing} = nothing)
  fl, h = _reduce_hilbert_series_by_primary_degrees(R, [ total_degree(f.f) for f in primary_invariants(R) ], chi)
  @assert fl
  return h
end

function _reduce_hilbert_series_by_primary_degrees(R::InvRing, degrees::Vector{Int}, chi::Union{GAPGroupClassFunction, Nothing} = nothing)
  mol = molien_series(R, chi)
  f = numerator(mol)
  g = denominator(mol)
  for d in degrees
    # multiply f by 1 - t^d
    f -= shift_left(f, d)
  end
  return divides(f, g)
end

# Return possible degrees of primary invariants d_1, ..., d_n with
# d_1 \cdots d_n == k*|G|, where G = group(R).
# (Note that |G| must divide d_1 \cdots d_n, see DK15, Prop. 3.5.5.)
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

    if mod(lcm(ds), exponent(group(R))) != 0
      continue
    end

    for d in ds
      # Check whether there exist invariants of this degree.
      if dimension_via_molien_series(fmpz, R, d) == 0
        skip = true
        break
      end
    end
    skip ? continue : nothing

    if is_molien_series_implemented(R)
      if !test_primary_degrees_via_hilbert_series(R, ds)
        continue
      end
    end

    push!(degrees, ds)
  end

  sort!(degrees, lt = (x, y) -> sum(x) < sum(y) || x < y)
  return degrees
end

# Checks whether there exist homogeneous f_1, ..., f_k in RG of degrees
# degrees[length(invars) + 1:length(invars) + k] such that
# RG/< invars, f_1, ..., f_k > has Krull dimension n - k, where n == length(invars).
# If the base field is finite, the answer "true" might be wrong (for theoretical reasons).
# See Kem99, Theorem 2.
function check_primary_degrees(RG::InvRing{FldT, GrpT, PolyElemT}, degrees::Vector{Int}, invars::Vector{PolyElemT}, k::Int, iters::Dict{Int, <: VectorSpaceIterator}, ideals::Dict{Set{PolyElemT}, Tuple{MPolyIdeal{PolyElemT}, Int}}) where {FldT, GrpT, PolyElemT}
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
    I, dimI = get!(ideals, gens_ideal) do
      I = ideal(R, collect(gens_ideal))
      d = dim(I)
      return (I, d)
    end
    if dimI > n - length(invars) - numbersInds
      return false
    end
  end
  return true
end

function _primary_invariants_via_optimal_hsop(RG::InvRing; ensure_minimality::Int = 0, degree_bound::Int = 1, primary_degrees::Vector{Int} = Int[])
  iters = Dict{Int, VectorSpaceIterator}()
  ideals = Dict{Set{elem_type(polynomial_ring(RG))}, Tuple{MPolyIdeal{elem_type(polynomial_ring(RG))}, Int}}()

  if !isempty(primary_degrees)
    @assert length(primary_degrees) == degree(group(RG))
    invars_cache = PrimaryInvarsCache{elem_type(polynomial_ring(RG))}()
    b, k = primary_invariants_via_optimal_hsop!(RG, primary_degrees, invars_cache, iters, ideals, ensure_minimality, 0)
    if !b
      error("No primary invariants of the given degrees exist")
    end
    return invars_cache
  end

  bad_prefixes = Vector{Vector{Int}}()
  l = degree_bound
  while true
    degrees = candidates_primary_degrees(RG, l, bad_prefixes)
    for ds in degrees
      invars_cache = PrimaryInvarsCache{elem_type(polynomial_ring(RG))}()
      b, k = primary_invariants_via_optimal_hsop!(RG, ds, invars_cache, iters, ideals, ensure_minimality, 0)
      b && return invars_cache
      push!(bad_prefixes, ds[1:k])
    end
    l += 1
  end
end

# Kemper "An Algorithm to Calculate Optimal Homogeneous Systems of Parameters", 1999, [Kem99]
# Return a bool b and an integer k.
# b == true iff primary invariants of the given degrees exist. In this case
# invars_cache will contain those invariants.
# k is only needed for recursive calls of the function.
function primary_invariants_via_optimal_hsop!(RG::InvRing{FldT, GrpT, PolyElemT}, degrees::Vector{Int}, invars_cache::PrimaryInvarsCache{PolyElemT}, iters::Dict{Int, <: VectorSpaceIterator}, ideals::Dict{Set{PolyElemT}, Tuple{MPolyIdeal{PolyElemT}, Int}}, ensure_minimality::Int = 0, k::Int = 0) where {FldT, GrpT, PolyElemT}

  n = length(degrees) - length(invars_cache.invars)
  R = polynomial_ring(RG)

  if length(invars_cache.invars) != length(degrees)
    d = degrees[length(invars_cache.invars) + 1]

    iter = get!(iters, Int(d)) do
      vector_space_iterator(coefficient_ring(RG), iterate_basis(RG, Int(d)))
    end

    attempts = 1
    for f in iter
      if f in invars_cache.invars
        continue
      end
      if ensure_minimality > 0 && attempts > ensure_minimality
        break
      end
      attempts += 1
      push!(invars_cache.invars, f)
      if k > 0
        if !check_primary_degrees(RG, degrees, invars_cache.invars, k - 1, iters, ideals)
          pop!(invars_cache.invars)
          continue
        end
      end
      b, kk = primary_invariants_via_optimal_hsop!(RG, degrees, invars_cache, iters, ideals, ensure_minimality, max(0, k - 1))
      if b
        return true, 0
      end
      pop!(invars_cache.invars)
      if kk >= k
        k = kk + 1
        if !check_primary_degrees(RG, degrees, invars_cache.invars, k, iters, ideals)
          break
        end
      end
    end
  end

  if k <= 1
    k = 1
    I, dimI = get!(ideals, Set(invars_cache.invars)) do
      I = ideal(polynomial_ring(RG), invars_cache.invars)
      d = dim(I)
      return (I, d)
    end
    if dimI > n
      k = 0
    end
  end
  if is_zero(n) && is_one(k)
    invars_cache.ideal = ideals[Set(invars_cache.invars)][1]
    return true, k
  end
  return false, k
end

@doc Markdown.doc"""
    primary_invariants(IR::InvRing;
      ensure_minimality::Int = 0, degree_bound::Int = 1,
      primary_degrees::Vector{Int} = Int[])

Return a system of primary invariants for `IR` as a `Vector` sorted by increasing
degree. The result is cached, so calling this function again with argument `IR`
will be fast and give the same result.

The primary invariants are computed using the algorithm in [Kem99](@cite).

The product of the degrees $d_1,\dots, d_n$ of the returned primary invariants
is guaranteed to be minimal among all possible sets of primary invariants.

Expert users (or users happy to experiment) may enter the following keyword arguments to
speed up the computation. Note that all of these options are ignored if there are already
primary invariants cached.
If admissible degrees $d_1,\dots, d_n$ for a system of primary invariants are known a
priori, these degrees can be specified by `primary_degrees = [d_1, ..., d_n]`.
Note that an error is raised if in fact no primary invariants of the given degrees exist.
An a priori known number $k \geq 1$ with $d_1\cdots d_n \geq k \cdot |G|$, where
$G$ is the underlying group, can be specified by `degree_bound = k`. The default value is
`degree_bound = 1`.
In some situations, the runtime of the algorithm might be improved by assigning
a positive integer to `ensure_minimality`. This leads to an early cancellation of
loops in the algorithm and the described minimality of the degrees is not
guaranteed anymore. A smaller (positive) value of `ensure_minimality` corresponds
to an earlier cancellation. However, the default value `ensure_minimality = 0`
corresponds to no cancellation.

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a");

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0]);

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1]);

julia> G = MatrixGroup(3, K, [M1, M2]);

julia> IR = invariant_ring(G);

julia> primary_invariants(IR)
3-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 x[1]*x[2]*x[3]
 x[1]^3 + x[2]^3 + x[3]^3
 x[1]^3*x[2]^3 + x[1]^3*x[3]^3 + x[2]^3*x[3]^3

julia> IR = invariant_ring(G); # "New" ring to avoid caching

julia> primary_invariants(IR, primary_degrees = [ 3, 6, 6 ])
3-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 x[1]*x[2]*x[3]
 x[1]^3*x[2]^3 + x[1]^3*x[3]^3 + x[2]^3*x[3]^3
 x[1]^6 + x[2]^6 + x[3]^6

```
"""
function primary_invariants(RG::InvRing; ensure_minimality::Int = 0, degree_bound::Int = 1, primary_degrees::Vector{Int} = Int[])
  if !isdefined(RG, :primary)
    RG.primary = _primary_invariants_via_optimal_hsop(RG; ensure_minimality = ensure_minimality, degree_bound = degree_bound, primary_degrees = primary_degrees)
  end
  return copy(RG.primary.invars)
end

# Access the (possibly) cached ideal generated by the primary invariants
function ideal_of_primary_invariants(RG::InvRing)
  _ = primary_invariants(RG)
  if !isdefined(RG.primary, :ideal)
    RG.primary.ideal = ideal(polynomial_ring(RG), RG.primary.invars)
  end
  return RG.primary.ideal
end
