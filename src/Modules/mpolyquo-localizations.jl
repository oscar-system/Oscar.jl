clear_denominators(A::MatrixType) where {T<:MPolyQuoLocalizedRingElem, MatrixType<:MatrixElem{T}} = clear_denominators(map_entries(lift, A))

# use the lifts to do computations effectively over the base ring
# rather than the quotient ring.
function clear_denominators(v::FreeModElem{<:MPolyQuoLocalizedRingElem})
  F = parent(v)
  L = base_ring(F)
  R = base_ring(L)
  d = lcm(lifted_denominator.(Vector(v)))
  u = elem_type(R)[]
  for a in Vector(v)
    push!(u, lifted_numerator(a)*div(d, lifted_denominator(a)))
  end
  Fb = base_ring_module(F)
  return sum([a*e for (a, e) in zip(u, gens(Fb))]), d
end

# constructs the matrix for the submodule I*F where 
# F is the free module of rank n over the base ring R 
# and I is the modulus.
function modulus_matrix(L::MPolyQuoLocalizedRing, n::Int)
  I = modulus(L)
  m = length(gens(I))
  R = base_ring(L)
  B = zero(MatrixSpace(R, 0, n))
  for j in 1:n
    A = zero(MatrixSpace(R, m, n))
    for i in 1:m
      A[i, j] = gens(I)[i]
    end
    B = vcat(B, A)
  end
  return B
end

function syz(A::MatrixElem{<:MPolyQuoLocalizedRingElem})
  B, D = clear_denominators(A)
  L = syz(vcat(B, modulus_matrix(base_ring(A), ncols(B))))
  return L*D
end

function ann(b::MatrixType, A::MatrixType) where {T<:MPolyQuoLocalizedRingElem, MatrixType<:MatrixElem{T}}
  R = base_ring(A)
  R === base_ring(b) || error("matrices must be defined over the same ring")
  nrows(b) == 1 || error("only matrices with one row are allowed!")
  B = vcat(A, modulus_matrix(R, ncols(A)))
  m = nrows(B)
  n = ncols(B)
  Aext = vcat(b, B)
  L = syz(Aext)
  return ideal(R, vec(L[:, 1]))
end

function has_solution(A::MatrixType, b::MatrixType) where {T<:MPolyQuoLocalizedRingElem, MatrixType<:MatrixElem{T}}
  S = base_ring(A)
  R = base_ring(S)
  S === base_ring(b) || error("matrices must be defined over the same ring")
  nrows(b) == 1 || error("only matrices with one row are allowed!")
  Aext = vcat(A, modulus_matrix(S, ncols(A)))
  B, D = clear_denominators(Aext)
  c, u = clear_denominators(b)
  (success, y, v) = has_solution(B, c, inverted_set(S))
  success || return (false, zero(MatrixSpace(S, 1, ncols(b))))
  # We have B = D⋅Aext and c = u ⋅ b as matrices. 
  # Now [y z]⋅B = v⋅c ⇔ [y z]⋅D ⋅Aext = v ⋅ u ⋅ b ⇔ v⁻¹ ⋅ u⁻¹ ⋅ [y z] ⋅ D ⋅ Aext = b.
  # Take v⁻¹ ⋅ u⁻¹ ⋅ [y z] ⋅ D to be the solution x of x ⋅ Aext = b. 
  # Then for the first m components x' of x we have x' ⋅ A ≡ b mod I
  x = S(one(R), v*u[1,1])*S(y*D)
  xpart = MatrixSpace(S, 1, nrows(A))
  for i in 1:nrows(A)
    xpart[1, i] = x[1, i]
  end
  return (success, xpart)
end

function pre_saturated_module(M::SubQuo{T}) where {T<:MPolyQuoLocalizedRingElem}
  has_attribute(M, :saturated_module) && return get_attribute(M, :saturated_module)::SubQuo{base_ring_elem_type(T)}
  if !has_attribute(M, :pre_saturated_module)
    S = base_ring(M)
    R = base_ring(S)
    (A, D) = clear_denominators(generator_matrix(M))
    relM = relations_matrix(M)
    (B, E) = clear_denominators(relM)
    mod_mat = modulus_matrix(S, ncols(B))
    B = vcat(B, mod_mat)
    E = vcat(E, zero(MatrixSpace(R, nrows(mod_mat), ncols(E))))
    F = ambient_free_module(M)
    Fb = base_ring_module(F)
    Mb = SubQuo(Fb, A, B)
    set_attribute!(M, :pre_saturation_data_gens, change_base_ring(S, D))
    set_attribute!(M, :pre_saturation_data_rels, change_base_ring(S, E))
    set_attribute!(M, :pre_saturated_module, Mb)
  end
  return get_attribute(M, :pre_saturated_module)::SubQuo{base_ring_elem_type(T)}
end

