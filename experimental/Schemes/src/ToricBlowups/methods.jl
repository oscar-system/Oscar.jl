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

@doc raw"""
    _minimal_supercone(X::NormalToricVariety, r::AbstractVector{<:IntegerUnion})

Given a ray $r$ inside the support of the simplicial fan $Σ$ of a normal toric variety~$X$, not coinciding with any ray in $Σ(1)$, return the unique cone $σ$ in $Σ$ such that $r$ is in the relative interior of $σ$.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> r = [1, 1, 0]
2-element Vector{Int64}:
 1
 1

julia> Oscar._minimal_supercone(X, r)
Polyhedral cone in ambient dimension 3

julia> rays(c)
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0, 0]
 [0, 1, 0]
```
"""
function _minimal_supercone(X::NormalToricVariety, r::AbstractVector{<:Union{IntegerUnion, QQFieldElem}})
  @assert is_orbifold(X) "Only implemented when fan is simplicial"
  contained_in_support_of_fan = false
  mcone = maximal_cones(X)[1]  # initialize `mcone`, fixing its type

  mcones = maximal_cones(X)
  while true
    contained_in_subcone = false
    for c in mcones
      if r in c
        contained_in_support_of_fan = true
        contained_in_subcone = true
        mcone = c
      end
    end
    @assert contained_in_support_of_fan "Ray not contained in the support of the fan"
    !contained_in_subcone && return mcone
    n_rays(mcone) == 1 && return mcone

    # set mcones to be the maximal subcones of mcone
    mcones = [cone(setdiff(collect(rays(mcone)), [ray])) for ray in rays(mcone)]
  end
end

function _cox_ring_homomorphism(f::ToricBlowupMorphism)
  @assert is_normal(codomain(f))
  @assert is_orbifold(codomain(f)) "Only implemented when fan is simplicial"
  Y = domain(f)
  @assert !has_torusfactor(codomain(f)) "Only implemented when there are no torus factors"
  @assert n_rays(Y) == n_rays(codomain(f)) + 1 "Only implemented when the blowup adds a ray"
  @assert grading_group(cox_ring(Y)) === class_group(Y)
  G = class_group(Y)
  n = number_of_generators(G)
  new_ray = rays(domain(f))[index_of_new_ray(f)]
  c = Oscar._minimal_supercone(Y, new_ray)
  U = affine_normal_toric_variety(c)
  d = order(class_group(U))
  # TODO: finish this...


  new_var = cox_ring(Y)[index_of_new_ray(f)]
  M = generator_degrees(cox_ring(Y))
  @assert M[index_of_new_ray(f)][n] < 0 "Assuming the blowup adds a new column vector to `Oscar.generator_degrees(cox_ring(codomain(f)))`, the entry of that vector corresponding to the new ray should be negative. If not, multiply by -1."
  ind = Int64(abs(M[index_of_new_ray(f)][n]))
  for i in 1:n_rays(Y)
    i == index_of_new_ray(f) || @assert M[i][n] >= 0 "Assuming the blowup adds a new column vector to `Oscar.generator_degrees(cox_ring(codomain(f)))` for which the entry corresponding to the new ray is negative, the entries corresponding to all the other rays should be nonnegative. If not, add integer multiples of the other columns until it is."
  end
  if ind == 1
    S = cox_ring(Y)
    inj = hom(S, S, gens(S))
  elseif ind > 1
    # Change grading so that the negative entry would be -1
    xs = ZZRingElem[]
    for i in 1:n-1
      push!(xs, M[index_of_new_ray(f)][i])
    end
    push!(xs, ZZRingElem(-1))
    M_new = deepcopy(M)
    M_new[index_of_new_ray(f)] = G(xs)
    S_ungraded = forget_grading(cox_ring(Y))
    S = grade(S_ungraded, M_new)[1]
    inj = hom(cox_ring(Y), S, gens(S))
  end
  F = hom(cox_ring(codomain(f)), S, [inj(new_var^(M[i][n])*S[i]) for i in 1:n_rays(Y) if i != index_of_new_ray(f)])
  return (F, ind)
end

function strict_transform(f::ToricBlowupMorphism, I::MPolyIdeal)
  F, ind = _cox_ring_homomorphism(f)
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
