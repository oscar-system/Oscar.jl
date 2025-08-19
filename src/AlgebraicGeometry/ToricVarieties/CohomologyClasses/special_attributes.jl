@doc raw"""
    cohomology_ring(v::NormalToricVarietyType; completeness_check::Bool = true)

!!! note "Simplicial and complete toric varieties"
    This function assumes that the toric variety is both **simplicial** and **complete**.
    Completeness checks may be slow; to skip them, use the optional keyword argument `completeness_check = false`.

# Examples
```jldoctest
julia> p2 = projective_space(NormalToricVariety, 2);

julia> ngens(cohomology_ring(p2))
3
```
"""
@attr Any function cohomology_ring(v::NormalToricVarietyType; completeness_check::Bool = true)
  @req is_simplicial(v) "The cohomology ring is only supported for simplicial and complete toric varieties"
  if completeness_check
    @req is_complete(v) "The cohomology ring is only supported for simplicial and complete toric varieties"
  end
  R, _ = graded_polynomial_ring(coefficient_ring(v), coordinate_names(v); cached=false)
  linear_relations = ideal_of_linear_relations(R, v)
  stanley_reisner = stanley_reisner_ideal(R, v)
  return quo(R, linear_relations + stanley_reisner)[1]
end


@doc raw"""
    volume_form(v::NormalToricVariety)

Construct the volume form of the normal toric toric variety `v`.

# Examples
```jldoctest
julia> polynomial(volume_form(projective_space(NormalToricVariety, 2)))
x3^2

julia> polynomial(volume_form(del_pezzo_surface(NormalToricVariety, 3)))
-e3^2

julia> polynomial(volume_form(hirzebruch_surface(NormalToricVariety, 5)))
1//5*x2^2
```
"""
@attr CohomologyClass function volume_form(v::NormalToricVariety)
    mc = ray_indices(maximal_cones(v))
    exponents = [ZZRingElem(mc[1, i]) for i in 1:length(mc[1, :])]
    indets = gens(cohomology_ring(v))
    poly = prod(indets[k]^exponents[k] for k in 1:length(exponents))
    @req !iszero(poly) && degree(poly)[1] == dim(v) "The volume class does not exist"
    return CohomologyClass(v, poly, true)
end


@attr Any function _intersection_form_via_exponents(v::NormalToricVariety)
    # extract the cohomology classes corresponding to the torus-invariant prime divisors
    generators = [cohomology_class(d) for d in torusinvariant_prime_divisors(v)]
    
    # find combinations of those classes that we have to integrate
    S, _ = graded_polynomial_ring(QQ, "g#" => 1:length(generators); cached=false)
    hc = homogeneous_component(S, [dim(v)])
    monoms = [hc[2](x) for x in gens(hc[1])]
    combinations = reduce(vcat, [[[ZZRingElem(l) for l in k] for k in AbstractAlgebra.exponent_vectors(m)] for m in monoms])
    
    # perform the integrals
    intersection_dict = Dict{ZZMatrix, QQFieldElem}()
    for expos in combinations
        cc = prod(generators[k]^expos[k] for k in 1:length(generators))
        intersection_dict[matrix(ZZ, [expos])] = integrate(cc)
    end
    
    # return the dictionary
    return intersection_dict
end


@doc raw"""
    intersection_form(v::NormalToricVariety)

Compute the intersection numbers among the cohomology classes
associated to the torusinvariant prime divisors of the normal toric toric variety `v`.

# Examples
```jldoctest
julia> F3 = hirzebruch_surface(NormalToricVariety, 3)
Normal toric variety

julia> length(intersection_form(F3))
10
```
"""
function intersection_form(v::NormalToricVariety)
    intersection_dict = _intersection_form_via_exponents(v)
    monoms_of_prime_divisors = gens(cox_ring(v))
    intersection_dict_for_user = Dict{MPolyRingElem, QQFieldElem}()
    for (expos, v) in intersection_dict
        monom = prod(monoms_of_prime_divisors[k]^expos[k] for k in 1:length(monoms_of_prime_divisors))
        intersection_dict_for_user[monom] = v
    end
    return intersection_dict_for_user
end


@doc raw"""
    chern_class(v::NormalToricVariety, k::Int; completeness_check::Bool = true)

Compute the `k`-th Chern class of the tangent bundle of a normal toric variety.  
The algorithm is based on Proposition 13.1.2 in [CLS11](@cite).

!!! note "Smooth and complete toric varieties"
    This function assumes that the toric variety is both **smooth** and **complete**.
    Completeness checks may be slow; to skip them, use the optional keyword argument `completeness_check = false`.

# Examples
```jldoctest
julia> F3 = hirzebruch_surface(NormalToricVariety, 3)
Normal toric variety

julia> chern_class(F3, 0)
Cohomology class on a normal toric variety given by 1

julia> chern_class(F3, 1, completeness_check = false)
Cohomology class on a normal toric variety given by t1 + x1 + t2 + x2

julia> integrate(chern_class(F3, 2), completeness_check = false)
4
```
"""
function chern_class(v::NormalToricVariety, k::Int; completeness_check::Bool = true)
  # Consistency checks
  @req k >= 0 "Chern class index must be non-negative"
  @req k <= dim(v) "Chern class index must not exceed dimension of the toric variety"

  # If thus far, no non-trivial Chern classes have been computed for this toric variety, add an "empty" vector and catch degenerate cases
  if !has_attribute(v, :chern_classes)
    cs = Dict{Int64, CohomologyClass}()
    cs[0] = cohomology_class(v, one(cohomology_ring(v; completeness_check)))
    cs[1] = cohomology_class(v, sum(gens(cohomology_ring(v; completeness_check))))
    set_attribute!(v, :chern_classes, cs)
    if k == 0
      return cs[0]
    elseif k == 1
      return cs[1]
    end
  end
  
  # Check if the Chern class in question is known
  cs = get_attribute(v, :chern_classes)::Dict{Int64, CohomologyClass}
  if haskey(cs, k)
    return cs[k]
  end

  # Check if we can compute the Chern classes for this toric variety
  @req is_smooth(v) "The Chern classes of the tangent bundle are only supported for smooth and complete toric varieties"
  if completeness_check
    @req is_complete(v) "The Chern classes of the tangent bundle are only supported for smooth and complete toric varieties"
  end
  
  # Preparation to check if we can discard a set in the following iteration that computes the Chern class
  mnf = _minimal_nonfaces(v)
  indices = [Set(Vector{Int}(Polymake.row(mnf,i))) for i in 1:Polymake.nrows(mnf)]
  function can_be_ignored(my_set)
    for k in 1:length(indices)
      if is_subset(indices[k], my_set)
        return true
      end
    end
    return false
  end

  # Compute, set and return the desired Chern class
  c_ring = cohomology_ring(v; completeness_check)
  b_ring = base_ring(c_ring)
  number_of_my_gens = ngens(b_ring)
  my_builder = MPolyBuildCtx(b_ring)
  for t in subsets(Set([i for i in 1:ngens(b_ring)]), k)
    can_be_ignored(t) && continue
    exps = zeros(Int64, number_of_my_gens)
    for a in t
      exps[a] = 1
    end
    push_term!(my_builder, one(QQ), exps)
  end
  desired_class = CohomologyClass(v, cohomology_ring(v; completeness_check)(finish(my_builder)), true)
  cs[k] = desired_class
  set_attribute!(v, :chern_classes, cs)
  return cs[k]
end


@doc raw"""
    chern_classes(v::NormalToricVariety; completeness_check::Bool = true)

Compute all Chern classes of the tangent bundle of a normal toric variety.
The algorithm is based on Proposition 13.1.2 in [CLS11](@cite).

!!! note "Smooth and complete toric varieties"
    This function assumes that the toric variety is both **smooth** and **complete**.
    Completeness checks may be slow; to skip them, use the optional keyword argument `completeness_check = false`.

# Examples
```jldoctest
julia> F3 = hirzebruch_surface(NormalToricVariety, 3)
Normal toric variety

julia> cs = chern_classes(F3);

julia> cs[0]
Cohomology class on a normal toric variety given by 1

julia> cs[1]
Cohomology class on a normal toric variety given by t1 + x1 + t2 + x2

julia> integrate(cs[2])
4
```
"""
function chern_classes(v::NormalToricVariety; completeness_check::Bool = true)
  @req is_smooth(v) "The Chern classes of the tangent bundle are only supported for smooth and complete toric varieties"
  if completeness_check
    @req is_complete(v) "The Chern classes of the tangent bundle are only supported for smooth and complete toric varieties"
  end
  for k in 0:dim(v)
    chern_class(v, k; completeness_check)
  end
  return get_attribute(v, :chern_classes)::Dict{Int64, CohomologyClass}
end


@doc raw"""
    basis_of_h4(v::NormalToricVariety; completeness_check::Bool = true)

Compute a monomial basis of the cohomology group $H^4(X, \mathbb{Q})$. The algorithm,
based on Theorem 12.4.1 in [CLS11](@cite), truncates the cohomology ring to degree 2.

!!! note "Simplicial and complete toric varieties"
    This function assumes that the toric variety is both **simplicial** and **complete**.
    Completeness checks may be slow; to skip them, use the optional keyword argument `completeness_check = false`.

# Examples
```jldoctest
julia> Y1 = hirzebruch_surface(NormalToricVariety, 2)
Normal toric variety

julia> Y2 = hirzebruch_surface(NormalToricVariety, 2)
Normal toric variety

julia> Y = Y1 * Y2
Normal toric variety

julia> h4_basis = basis_of_h4(Y)
6-element Vector{CohomologyClass}:
 Cohomology class on a normal toric variety given by yx2^2
 Cohomology class on a normal toric variety given by xx2*yx2
 Cohomology class on a normal toric variety given by xx2*yt2
 Cohomology class on a normal toric variety given by xx2^2
 Cohomology class on a normal toric variety given by xt2*yx2
 Cohomology class on a normal toric variety given by xt2*yt2

julia> betti_number(Y, 4) == length(h4_basis)
true
```
"""
@attr Vector{CohomologyClass} function basis_of_h4(v::NormalToricVariety; completeness_check::Bool = true)
  @req is_simplicial(v) "Computation of basis of H4(X, Q) is currently only supported for simplicial toric varieties"
  if completeness_check
    @req is_complete(v) "Computation of basis of H4(X, Q) is currently only supported for complete toric varieties"
  end
  if dim(v) < 4
    return Vector{CohomologyClass}()
  end
  R = cohomology_ring(v; completeness_check)
  basis_of_h4 = [cohomology_class(v, R(g)) for g in monomial_basis(R, [2])]
  return basis_of_h4
end
