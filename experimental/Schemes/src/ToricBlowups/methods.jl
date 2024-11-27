@doc raw"""
    total_transform(f::AbsSimpleBlowupMorphism, II::IdealSheaf)

Computes the total transform of an ideal sheaf along a blowup.

In particular, this applies in the toric setting. However, note that
currently (October 2023), ideal sheaves are only supported on smooth
toric varieties.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> bl = blow_up(P2, [1, 1])
Toric blowup morphism

julia> S = cox_ring(P2);

julia> x, y, z = gens(S);

julia> I = ideal_sheaf(P2, ideal([x*y]))
Sheaf of ideals
  on normal, smooth toric variety
with restrictions
  1: Ideal (x_1_1*x_2_1)
  2: Ideal (x_2_2)
  3: Ideal (x_1_3)

julia> total_transform(bl, I)
Sheaf of ideals
  on normal toric variety
with restrictions
  1: Ideal (x_1_1*x_2_1^2)
  2: Ideal (x_1_2^2*x_2_2)
  3: Ideal (x_2_3)
  4: Ideal (x_1_4)
```
"""
function total_transform(f::AbsSimpleBlowupMorphism, II::AbsIdealSheaf)
  return pullback(f, II)
end

function total_transform(f::AbsBlowupMorphism, II::AbsIdealSheaf)
  return pullback(f, II)
end

function cox_ring_group_homomorphism(f::ToricBlowupMorphism, g::MPolyDecRingElem{QQFieldElem, QQMPolyRingElem})
  @req parent(g) === cox_ring(codomain(f)) "g must be an element of the Cox ring of the codomain of f"
  d, phi_exponents = _cox_ring_group_homomorphism_data(f)
  S = cox_ring(domain(f))
  new_var = S[index_of_new_ray(f)]

  # Assuming the i-th variable of `cox_ring(X)` is the i-th variable of `cox_ring(Y)`
  S_vars = [S[i] for i in 1:nvars(S)]
  coeff_exps = collect(zip(coefficients(g), exponents(g)))
  make_S_term(c, exps) = c*prod(map(^, S_vars, exps))
  h = S(0)
  for (c, exps) in coeff_exps
    new_exp = div(sum(phi_exponents.*exps), d)
    sum(phi_exponents.*exps) != d && (new_exp += 1)
    h += new_var^new_exp * make_S_term(c, exps)
  end
  return h
end

@attr Tuple{ZZRingElem, Vector{ZZRingElem}} function _cox_ring_group_homomorphism_data(f::ToricBlowupMorphism)
  # f: Y -> X
  X = codomain(f)
  Y = domain(f)
  @req is_normal(X) "Only implemented when the variety is normal"
  @req is_orbifold(X) "Only implemented when the fan is simplicial"
  @req !has_torusfactor(X) "Only implemented when there are no torus factors"
  @req n_rays(Y) == n_rays(X) + 1 "Only implemented when the blowup adds a ray"
  new_ray = rays(domain(f))[index_of_new_ray(f)]
  c = minimal_supercone(X, new_ray)
  U = affine_normal_toric_variety(c)
  d = order(class_group(U))
  A = transpose(matrix(ZZ, rays(c)))
  C = identity_matrix(ZZ, n_rays(c))

  # Making `new_ray` into a singleton SubObjectIterator similar to `rays(c)`
  b_iterator = rays(positive_hull(new_ray))

  # Making `new_ray` into a column vector
  b = transpose(matrix(ZZ, b_iterator))

  # Multiplying by `d` to assure that a solution exists
  b *= d

  vs = solve_mixed(A, b, C)
  n = nvars(cox_ring(X))
  phi_exponents = zeros(ZZRingElem, n)
  for i in 1:n
    ray = Oscar._ray_variable_correspondence(X)[cox_ring(X)[i]]
    j = findfirst(isequal(ray), rays(c))
    if isnothing(j)
      phi_exponents[i] = ZZRingElem(0)
    else
      phi_exponents[i] = ZZRingElem(vs[j])
    end
  end
  return d, phi_exponents
end

function strict_transform(f::ToricBlowupMorphism, I::MPolyIdeal)
  F, ind = _cox_ring_group_homomorphism(f)
  Y = domain(f)
  S = cox_ring(Y)
  new_var = S[index_of_new_ray(f)]
  J = saturation(F(I), ideal(new_var))
  if ind == 1
    return J
  elseif ind > 1
    xs = elem_type(S)[]
    for gen in gens(J)
      new_exponent_vectors = Vector{Int64}[]
      for i in 1:length(terms(gen))
        exp_vec = exponent_vector(gen, i)
        bool, result = divides(exp_vec[index_of_new_ray(f)], ind)
        @assert bool "Given the finite set of generators of `J`, the exponent of `new_var` in every monomial of every member should be divisible by `ind`. If not, first make sure that none of the generators of `J` is divisible by `new_var`."
        exp_vec[index_of_new_ray(f)] = result
        push!(new_exponent_vectors, exp_vec)
      end
      push!(xs, S(collect(coefficients(gen)), new_exponent_vectors))
    end
    return ideal(cox_ring(Y), xs)
  end
end
