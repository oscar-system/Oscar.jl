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
      return false, sol, sol
    end
    if task !== :only_check
      sol[i, :] = dense_row(coordinates(preimage(h, c))::sparse_row_type(R), nrows(A))
    end
  end
  if task !== :with_kernel
    return true, sol, sol
  end
  K, KtoF = kernel(h)
  if ngens(K) == 0
    Ke = zero_matrix(R, 0, nrows(A))
  else
    Ke = reduce(vcat, dense_row.((coordinates.(KtoF.(gens(K))))::Vector{sparse_row_type(R)}, rank(F)))
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
  fl::Bool, sol::typeof(AA), K::typeof(AA) = AbstractAlgebra.Solve._can_solve_internal_no_check(AbstractAlgebra.Solve.matrix_normal_form_type(AA), AA, bb, task; side)
  if task === :only_check
    return fl, A, A
  elseif task === :with_solution
    return fl::Bool, map_entries(f.inv, xx)::typeof(A), A
  else
    return fl::Bool, map_entries(f.inv, xx)::typeof(A), map_entries(f.inv, kk)::typeof(b)
  end
end
