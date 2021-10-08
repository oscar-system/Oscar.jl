export solve_non_negative, solve_mixed, solve_ineq


@doc Markdown.doc"""
    solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)

Solves $Ax = b$ under $Cx >= 0$, assumes a finite solution set.
"""
function solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)  # Ax == b && Cx >= 0
    eq = HalfspaceIterator(Matrix{BigInt}(A), vec(Matrix{BigInt}(b)))
    ineq = HalfspaceIterator(Matrix{BigInt}(-C), vec(Matrix{BigInt}(zero_matrix(FlintZZ, nrows(C), 1))))
    P = Polyhedron(ineq, eq)
    inner = interior_lattice_points(P)
    outer = boundary_lattice_points(P)

    result = zero_matrix(FlintZZ, length(inner) + length(outer), ambient_dim(P))
    i = 1
    for v in inner
        result[i, :] = Vector{BigInt}(v.p)
        i += 1
    end
    for v in outer
        result[i, :] = Vector{BigInt}(v.p)
        i += 1
    end
    return result
end


@doc Markdown.doc"""
    solve_ineq(A::fmpz_mat, b::fmpz_mat)

Solves $Ax<=b$, assumes finite set of solutions.
"""
function solve_ineq(A::fmpz_mat, b::fmpz_mat)
    return solve_mixed(zero_matrix(FlintZZ, 0, ncols(A)), zero_matrix(FlintZZ,0,1), hcat(b, -A))
end


@doc Markdown.doc"""
    solve_non_negative(A::fmpz_mat, b::fmpz_mat)

Finds all solutions to $Ax = b$, $x>=0$. Assumes a finite set of solutions.
"""
function solve_non_negative(A::fmpz_mat, b::fmpz_mat)
  return solve_mixed(A, b, identity_matrix(FlintZZ, ncols(A)))
end


@doc Markdown.doc"""
    solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat)

Solves $Ax = b$ under $Cx >= d$, assumes a finite solution set.
"""
function solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat)
  n = ncols(A)
  A = cat(A, identity_matrix(FlintZZ, ncols(d)), dims=(1,2))
  b = vcat(b, identity_matrix(FlintZZ, ncols(d)))
  C = [C -d; zero_matrix(FlintZZ, ncols(d), ncols(C)) identity_matrix(FlintZZ, ncols(d))]
  s = solve_mixed(A, b, C)
  return s[:, 1:n]
end
