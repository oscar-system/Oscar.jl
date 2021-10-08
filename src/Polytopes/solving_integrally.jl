export solve_non_negative, solve_mixed, solve_ineq


@doc Markdown.doc"""
    solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat)

Solve $Ax = b$ under $Cx >= d$, assumes a finite solution set.

# Examples
Find all $(x_1, x_2)\in\mathbb{Z}^2$ such that $x_1+x_2=7$, $x_1\ge 2$, and $x_2\ge 3$.
```jldoctest
julia> A = fmpz_mat([1 1]);

julia> b = zero_matrix(FlintZZ, 1,1); b[1,1]=7;

julia> C = fmpz_mat([1 0; 0 1]);

julia> d = zero_matrix(FlintZZ,2,1); d[1,1]=2; d[2,1]=3;

julia> solve_mixed(A,b,C,d)
[3   4]
[2   5]
[4   3]
```
"""
function solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat)
    eq = HalfspaceIterator(Matrix{BigInt}(A), vec(Matrix{BigInt}(b)))
    ineq = HalfspaceIterator(Matrix{BigInt}(-C), vec(Matrix{BigInt}(-d)))
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
    solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)

Solve $Ax = b$ under $Cx >= 0$, assumes a finite solution set.

# Examples
Find all $(x_1, x_2)\in\mathbb{Z}^2_{\ge 0}$ such that $x_1+x_2=3$.
```jldoctest
julia> A = fmpz_mat([1 1]);

julia> b = zero_matrix(FlintZZ, 1,1); b[1,1]=3;

julia> C = fmpz_mat([1 0; 0 1]);

julia> solve_mixed(A,b,C)
[1   2]
[2   1]
[0   3]
[3   0]
```
"""
function solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)  # Ax == b && Cx >= 0
    return solve_mixed(A, b, C, zero_matrix(FlintZZ, nrows(C), 1))
end


@doc Markdown.doc"""
    solve_ineq(A::fmpz_mat, b::fmpz_mat)

Solve $Ax<=b$, assumes finite set of solutions.

# Examples
The following gives the vertices of the square.
```jldoctest
julia> A = fmpz_mat([1 0; 0 1; -1 0; 0 -1]);

julia> b = zero_matrix(FlintZZ, 4,1); b[1,1]=1; b[2,1]=1; b[3,1]=0; b[4,1]=0;

julia> solve_ineq(A,b)
[0   0]
[0   1]
[1   0]
[1   1]
```
"""
function solve_ineq(A::fmpz_mat, b::fmpz_mat)
    return solve_mixed(zero_matrix(FlintZZ, 0, ncols(A)), zero_matrix(FlintZZ,0,1), -A, -b)
end


@doc Markdown.doc"""
    solve_non_negative(A::fmpz_mat, b::fmpz_mat)

Find all solutions to $Ax = b$, $x>=0$. Assumes a finite set of solutions.

# Examples
Find all $(x_1, x_2)\in\mathbb{Z}^2_{\ge 0}$ such that $x_1+x_2=3$.
```jldoctest
julia> A = fmpz_mat([1 1]);

julia> b = zero_matrix(FlintZZ, 1,1); b[1,1]=3;

julia> solve_non_negative(A,b)
[1   2]
[2   1]
[0   3]
[3   0]
```
"""
function solve_non_negative(A::fmpz_mat, b::fmpz_mat)
  return solve_mixed(A, b, identity_matrix(FlintZZ, ncols(A)))
end


