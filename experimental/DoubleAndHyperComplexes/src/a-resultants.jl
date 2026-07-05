function a_resultant_complex(support_sets::Vector{T};
    toric_variety::NormalToricVariety=_get_toric_variety(support_sets),
    inner_toric_ctx_object::NewToricCtx=NewToricCtx(toric_variety),
    parent::Ring=_parent_for_resultant(support_sets),
    outer_toric_ctx_object::ToricCtxWithParams=_get_outer_toric_ctx(inner_toric_ctx_object, parent),
    twist::FinGenAbGroupElem=zero(grading_group(cox_ring(toric_variety)))
  ) where {T}

  # partition the generators of `parent` to match the rows of the matrices 
  # in `support_sets`
  coef = Vector{elem_type(parent)}[]
  a = gens(parent)
  offset = 0
  for A in support_sets
    push!(coef, a[offset+1:offset+nrows(A)])
    offset += nrows(A)
  end

  S_ext = graded_ring(outer_toric_ctx_object)
  @assert parent === coefficient_ring(S_ext) "incompatible coefficient ring"
  x = gens(S_ext)

  # Ray generators in fan of toric variety
  U = map(primitive_generator, rays(toric_variety));

  # Coordinates of Cartier divisors in X of the system's Newton polytopes
  div_coords = [map(u -> -minimum( grad*primitive_generator(u) ), rays(toric_variety)) for grad in support_sets];

  # Homogeneous polynomials specified by characters in supports
  f = elem_type(S_ext)[sum(c[j]*prod(x[i]^(dot(A[j,:], U[i]) + dc[i]) for i in 1:n_rays(toric_variety); init=one(S_ext)) 
           for j in 1:size(A, 1); init=zero(S_ext)) for (A, c, dc) in zip(support_sets, coef, div_coords)];
  K = Oscar.HomogKoszulComplex(S_ext, f)
  t = Oscar.ZeroDimensionalComplex(graded_free_module(S_ext, [-twist]))
  Kt = tensor_product(K, t)
  return DirectImageComplex(outer_toric_ctx_object, Kt)
end

### helper functions for the above kwargs
function _parent_for_resultant(support_sets::Vector{T}) where {T}
  # create a polynomial ring for the generic coefficients of the support matrices
  R, a = polynomial_ring(ZZ, vcat([[Symbol("a_$(k)_$(i)") for i in 1:nrows(A)] for (k, A) in enumerate(support_sets)]...))
  return R
end

function _get_toric_variety(support_sets::Vector{T}) where {T}
  # Newton polytopes of support sets
  Q = map(convex_hull, support_sets);
  # Toric variety specified by the normal fan of the Minkowski sum of all Newton polytopes
  return normal_toric_variety(reduce(minkowski_sum, Q));
end

function _get_toric_ctx(support_sets::Vector{T}, X::NormalToricVariety) where {T}
  return 
end

function _get_outer_toric_ctx(inner_ctx::NewToricCtx, R::Ring)
  X = toric_variety(inner_ctx)
  S = cox_ring(X)
  kk = coefficient_ring(R)
  coef_map = MapFromFunc(QQ, R, x->R(kk(x)))
  S_ext, phi = change_base_ring(coef_map, S)
  return ToricCtxWithParams(inner_ctx, phi)
end

