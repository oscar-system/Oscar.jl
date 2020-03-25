export solve_non_negative, solve_mixed, solve_ineq

nrows(A::Polymake.MatrixAllocated) = Int(size(A)[1])
ncols(A::Polymake.MatrixAllocated) = Int(size(A)[2])

function _polytope(; A::fmpz_mat=zero_matrix(FlintZZ, 1, 1), b::fmpz_mat=zero_matrix(FlintZZ, ncols(A), 1), C::fmpz_mat=zero_matrix(FlintZZ, 1, 1))
  if !iszero(A)
    bA = Array{BigInt, 2}(hcat(-b, A))
    z = findall(i->!iszero_row(bA, i), 1:nrows(bA))
    zbA = Array{BigInt, 2}(bA[z, :])
  else
    zbA = Array{BigInt, 2}(undef, 0, 0)
  end
  if !iszero(C)
    z = findall(i->!iszero_row(C, i), 1:nrows(C))
    zI = Array{BigInt, 2}(hcat(zero_matrix(FlintZZ, nrows(C), 1), C))[z, :]
  else
    zI = Array{BigInt, 2}(undef, 0, 0)
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



@doc Markdown.doc"""
    solve_ineq(A::fmpz_mat, b::fmpz_mat)

Solves $Ax<=b$, assumes finite set of solutions.
"""
function solve_ineq(A::fmpz_mat, b::fmpz_mat)
  p = _polytope(C = hcat(b, -A))
  p.BOUNDED || error("not a bounded set")
  inner = p.INTERIOR_LATTICE_POINTS
  out = p.BOUNDARY_LATTICE_POINTS

  res = zero_matrix(FlintZZ, nrows(inner) + nrows(out), ncols(A))
  for i=1:nrows(out)
    @assert out[i,1] == 1
    for j=1:ncols(A)
      res[i,j] = out[i, j+1]
    end
  end
  for i=1:nrows(inner)
    @assert inner[i,1] == 1
    for j=1:ncols(A)
      res[i+nrows(out), j] = inner[i, j+1]
    end
  end
  return res
end

@doc Markdown.doc"""
    solve_non_negative(A::fmpz_mat, b::fmpz_mat)

Finds all solutions to $Ax = b$, $x>=0$. Assumes a finite set of solutions.
"""
function solve_non_negative(A::fmpz_mat, b::fmpz_mat)
  p = _polytope(A = A, b = b, C = identity_matrix(FlintZZ, ncols(A)))
  p.BOUNDED || error("not a bounded set")
  inner = p.INTERIOR_LATTICE_POINTS
  out = p.BOUNDARY_LATTICE_POINTS

  res = zero_matrix(FlintZZ, nrows(inner) + nrows(out), ncols(A))
  for i=1:nrows(out)
    @assert out[i,1] == 1
    for j=1:ncols(A)
      res[i,j] = out[i, j+1]
    end
  end
  for i=1:nrows(inner)
    @assert inner[i,1] == 1
    for j=1:ncols(A)
      res[i+nrows(out), j] = inner[i, j+1]
    end
  end
  return res
end

@doc Markdown.doc"""
    solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)

Solves $Ax = b$ under $Cx >= 0$, assumes a finite solution set.
"""
function solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)  # Ax == b && Cx >= 0
  p = _polytope(A = A, b = b, C = C)
  p.BOUNDED || error("not a bounded set")
  inner = p.INTERIOR_LATTICE_POINTS
  out = p.BOUNDARY_LATTICE_POINTS

  res = zero_matrix(FlintZZ, nrows(inner) + nrows(out), ncols(A))
  for i=1:nrows(out)
    if out[i,1] != 1
      println("unbounded polytope!!")
      global last_in = (A, b, C)
    end
    @assert out[i,1] == 1
    for j=1:ncols(A)
      res[i,j] = out[i, j+1]
    end
  end
  for i=1:nrows(inner)
    if inner[i,1] != 1
      println("unbounded polytope!!")
      global last_in = (A, b, C)
    end
    @assert inner[i,1] == 1
    for j=1:ncols(A)
      res[i+nrows(out), j] = inner[i, j+1]
    end
  end
  return res
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
