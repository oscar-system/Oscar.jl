struct ModuleSolveTrait <: AbstractAlgebra.Solve.MatrixNormalFormTrait end

AbstractAlgebra.Solve.matrix_normal_form_type(::Union{MPolyRing, MPolyQuoRing, MPolyLocRing, MPolyQuoLocRing}) = ModuleSolveTrait()

function AbstractAlgebra.Solve._can_solve_internal_no_check(::ModuleSolveTrait, A, b, task::Symbol; side::Symbol)
  if side === :right
    fl, _X, _K = AbstractAlgebra.Solve._can_solve_internal_no_check(ModuleSolveTrait(), transpose(A), transpose(b), task; side = :left)
    return fl, transpose(_X), transpose(_K)
  end

  # we solve x * A = b
  R = base_ring(A)

  F = FreeModule(R, nrows(A); cached = false)
  G = FreeModule(R, ncols(A); cached = false)
  h = hom(F, G, A)
  imh, = image(h)
  sol = zero_matrix(R, nrows(b), nrows(A))
  for i in 1:nrows(b)
    c = G(b[i, :])
    if !(c in imh)
      return false
    end
    sol[i, :] = dense_row(coordinates(preimage(h, c)), nrows(A))
  end
  K, KtoF = kernel(h)
  if ngens(K) == 0
    Ke = zero_matrix(R, 0, nrows(A))
  else
    Ke = reduce(vcat, dense_row.(coordinates.(KtoF.(gens(K))), rank(F)))
  end
  return true, sol, Ke
end

struct LaurentSolveTrait <: AbstractAlgebra.Solve.MatrixNormalFormTrait end

AbstractAlgebra.Solve.matrix_normal_form_type(::AbstractAlgebra.Generic.LaurentMPolyWrapRing) = LaurentSolveTrait()

function AbstractAlgebra.Solve._can_solve_internal_no_check(::LaurentSolveTrait, A, b, task::Symbol; side::Symbol)
  # just delegate to the underlying multivariate ring
  R = base_ring(A)
  f = Oscar._polyringquo(R)
  RQ = codomain(f)
  AA = map_entries(f, A)
  bb = map_entries(f, b)
  fl, xx, kk = can_solve_with_solution_and_kernel(AA, bb; side)
  return fl, map_entries(f.inv, xx), map_entries(f.inv, kk)
end
