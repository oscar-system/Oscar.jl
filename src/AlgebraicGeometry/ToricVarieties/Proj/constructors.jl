@doc raw"""
    projectivization(E::ToricLineBundle...)

This function computes the projectivization of a direct sum of line bundles or divisors.
Please see [OM78](@cite) for more background information.

# Examples
Let us construct the projective bundles ``X=\mathbb{P}(\mathcal{O}_{\mathbb{P}^1}\oplus\mathcal{O}_{\mathbb{P}^1}(1))`` and ``Y=\mathbb{P}(\mathcal{O}_{\mathbb{P}^1}\oplus\mathcal{O}_{\mathbb{P}^1}(2))``.
```jldoctest
julia> P1 = projective_space(NormalToricVariety, 1);

julia> D0 = toric_divisor(P1, [0,0]);

julia> D1 = toric_divisor(P1, [1,0]);

julia> X = projectivization(D0, D1)
Normal toric variety

julia> L0 = toric_line_bundle(P1, [0]);

julia> L1 = toric_line_bundle(P1, [2]);

julia> Y = projectivization(L0, L1)
Normal toric variety
```
"""
function projectivization(E::ToricLineBundle...)
  return _projectivization_and_total_space(true, [E...])
end

function projectivization(E::ToricDivisor...)
  return _projectivization_and_total_space(true, [E...])
end

function projectivization()
  error("The direct sum is empty")
end

function _projectivization_and_total_space(is_proj::Bool, E::Vector{T}) where T <: Union{ToricDivisor, ToricLineBundle}

  v = toric_variety(E[1])

  is_proj && length(E) == 1 && return v

  @req all(i -> toric_variety(E[i]) === v, eachindex(E)) "The divisors are defined on different toric varieties."

  if is_proj
    PF_fiber = normal_fan(simplex(length(E) - 1))
  else
    PF_fiber = polyhedral_fan(positive_hull(identity_matrix(ZZ, length(E))))
  end
  l = rays(PF_fiber)

  modified_ray_gens = Dict{RayVector{QQFieldElem}, RayVector{QQFieldElem}}()

  for sigma in maximal_cones(v)
    pol_sigma = polarize(sigma)
    for ray in rays(sigma)
      ray in keys(modified_ray_gens) && continue
      modified_ray_gens[ray] = vcat(ray, sum(i -> dot(_m_sigma(sigma, pol_sigma, E[i]), ray) * l[i], eachindex(E)))
    end
    length(keys(modified_ray_gens)) == n_rays(v) && break
  end

  new_maximal_cones = Vector{Vector{Int64}}(undef, n_maximal_cones(v) * n_maximal_cones(PF_fiber))
  index = 1

  for a in 1:n_maximal_cones(v)
    first = [row(ray_indices(maximal_cones(v)), a)...]
    for b in 1:n_maximal_cones(PF_fiber)
      second = [row(ray_indices(maximal_cones(PF_fiber)), b)...] .+ n_rays(v)
      new_maximal_cones[index] = vcat(first, second)
      index += 1
    end
  end

  total_rays_gens = vcat([modified_ray_gens[ray] for ray in rays(v)], [vcat(ray_vector(zeros(Int64, dim(v))), l[i]) for i in eachindex(E)])

  return normal_toric_variety(polyhedral_fan(IncidenceMatrix(new_maximal_cones), total_rays_gens))
end

function _m_sigma(sigma::Cone{QQFieldElem}, pol_sigma::Cone{QQFieldElem}, D::Union{ToricDivisor, ToricLineBundle})

  ans = ray_vector(zeros(QQFieldElem, dim(sigma)))
  coeff = coefficients(isa(D, ToricDivisor) ? D : toric_divisor(D))

  dual_ray = QQFieldElem[]

  for ray in rays(sigma)
    for pol_ray in rays(pol_sigma)
      dot(ray, pol_ray) == 0 && continue
      dual_ray = lcm(denominator.(pol_ray)) * pol_ray
      break
    end
    i = findfirst(==(ray), rays(toric_variety(D)))
    ans -= coeff[i] * dual_ray
  end

  return ans
end

@doc raw"""
    total_space(E::ToricLineBundle...)

This function computes the total space of a direct sum of line bundles or divisors.
Please see [OM78](@cite) for more background information.

# Examples
Let us construct the toric Calabi-Yau varieties given by the total space of ``\mathcal{O}_{\mathbb{P}^1}(2)\oplus\mathcal{O}_{\mathbb{P}^1}(-4)`` and ``\omega_{\mathbb{P}^2}``.
```jldoctest
julia> P1 = projective_space(NormalToricVariety, 1);

julia> L1 = toric_line_bundle(P1, [2]);

julia> L2 = toric_line_bundle(P1, [-4]);

julia> X = total_space(L1, L2)
Normal toric variety

julia> degree(canonical_bundle(X))
0

julia> P2 = projective_space(NormalToricVariety, 2);

julia> D = canonical_divisor(P2);

julia> Y = total_space(D)
Normal toric variety

julia> degree(canonical_bundle(Y))
0
```
"""
function total_space(E::ToricLineBundle...)
  return _projectivization_and_total_space(false, [E...])
end

function total_space(E::ToricDivisor...)
  return _projectivization_and_total_space(false, [E...])
end

function total_space()
  error("The direct sum is empty")
end
