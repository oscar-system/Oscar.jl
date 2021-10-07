export solve_non_negative, solve_mixed, solve_ineq

nrows(A::Polymake.MatrixAllocated) = Int(size(A)[1])
ncols(A::Polymake.MatrixAllocated) = Int(size(A)[2])

#TODO: define polytopes and cones properly!!!

function _polytope(; A::fmpz_mat=zero_matrix(FlintZZ, 1, 1), b::fmpz_mat=zero_matrix(FlintZZ, ncols(A), 1), C::fmpz_mat=zero_matrix(FlintZZ, 1, 1))
  if !iszero(A)
    bA = Matrix{BigInt}(hcat(-b, A))
    z = findall(i->!iszero_row(bA, i), 1:nrows(bA))
    zbA = Matrix{BigInt}(bA[z, :])
  else
    zbA = Matrix{BigInt}(undef, 0, 0)
  end
  if !iszero(C)
    z = findall(i->!iszero_row(C, i), 1:nrows(C))
    zI = Matrix{BigInt}(hcat(zero_matrix(FlintZZ, nrows(C), 1), C))[z, :]
  else
    zI = Matrix{BigInt}(undef, 0, 0)
  end
  if length(zbA) == 0
    p =  Polymake.polytope.Polytope(INEQUALITIES = zI)
  else
    if nrows(zI) == 0
      p = Polymake.polytope.Polytope(EQUATIONS = zbA)
    else
      p = Polymake.polytope.Polytope(EQUATIONS = zbA, INEQUALITIES = zI)
    end
  end
  return p
end


function _polytope_new(; A::fmpz_mat=zero_matrix(FlintZZ, 1, 1), b::fmpz_mat=zero_matrix(FlintZZ, ncols(A), 1), C::fmpz_mat=zero_matrix(FlintZZ, 1, 1))
    P = Polyhedron(_polytope(A=A, b=b, C=C))
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
    return _polytope_new(C = hcat(b, -A))
end


@doc Markdown.doc"""
    solve_non_negative(A::fmpz_mat, b::fmpz_mat)

Finds all solutions to $Ax = b$, $x>=0$. Assumes a finite set of solutions.
"""
function solve_non_negative(A::fmpz_mat, b::fmpz_mat)
  return _polytope_new(A = A, b = b, C = identity_matrix(FlintZZ, ncols(A)))
end


@doc Markdown.doc"""
    solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)

Solves $Ax = b$ under $Cx >= 0$, assumes a finite solution set.
"""
function solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)  # Ax == b && Cx >= 0
    return _polytope_new(A = A, b = b, C = C)
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
