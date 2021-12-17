export solve_non_negative, solve_mixed, solve_ineq


@doc Markdown.doc"""
    solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat)

Solve $Ax = b$ under $Cx >= d$, assumes a finite solution set.

# Examples
Find all $(x_1, x_2)\in\mathbb{Z}^2$ such that $x_1+x_2=7$, $x_1\ge 2$, and $x_2\ge 3$.
The solutions are the rows of the output.
Note that the output can be permuted, hence we sort it.
```jldoctest
julia> A = fmpz_mat([1 1]);

julia> b = zero_matrix(FlintZZ, 1,1); b[1,1]=7;

julia> C = fmpz_mat([1 0; 0 1]);

julia> d = zero_matrix(FlintZZ,2,1); d[1,1]=2; d[2,1]=3;

julia> sortslices(Matrix{BigInt}(solve_mixed(A, b, C, d)), dims=1)
3×2 Matrix{BigInt}:
 2  5
 3  4
 4  3
```
"""
function solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat)
    eq = (Matrix{BigInt}(A), vec(Matrix{BigInt}(b)))
    ineq = (Matrix{BigInt}(-C), vec(Matrix{BigInt}(-d)))
    P = Polyhedron(ineq, eq)
    LP = lattice_points(P)
    return transpose(matrix(ZZ, ambient_dim(P), length(LP), hcat(LP...)))
end


@doc Markdown.doc"""
    solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat)

Solve $Ax = b$ under $Cx >= 0$, assumes a finite solution set.

# Examples
Find all $(x_1, x_2)\in\mathbb{Z}^2_{\ge 0}$ such that $x_1+x_2=3$.
The solutions are the rows of the output.
Note that the output can be permuted, hence we sort it.
```jldoctest
julia> A = fmpz_mat([1 1]);

julia> b = zero_matrix(FlintZZ, 1,1); b[1,1]=3;

julia> C = fmpz_mat([1 0; 0 1]);

julia> sortslices(Matrix{BigInt}(solve_mixed(A, b, C)), dims=1)
4×2 Matrix{BigInt}:
 0  3
 1  2
 2  1
 3  0
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
The solutions are the rows of the output.
Note that the output can be permuted, hence we sort it.
```jldoctest
julia> A = fmpz_mat([1 0; 0 1; -1 0; 0 -1]);

julia> b = zero_matrix(FlintZZ, 4,1); b[1,1]=1; b[2,1]=1; b[3,1]=0; b[4,1]=0;

julia> sortslices(Matrix{BigInt}(solve_ineq(A, b)), dims=1)
4×2 Matrix{BigInt}:
 0  0
 0  1
 1  0
 1  1
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
The solutions are the rows of the output.
Note that the output can be permuted, hence we sort it.
```jldoctest
julia> A = fmpz_mat([1 1]);

julia> b = zero_matrix(FlintZZ, 1,1); b[1,1]=3;

julia> sortslices(Matrix{BigInt}(solve_non_negative(A, b)), dims=1)
4×2 Matrix{BigInt}:
 0  3
 1  2
 2  1
 3  0
```
"""
function solve_non_negative(A::fmpz_mat, b::fmpz_mat)
  return solve_mixed(A, b, identity_matrix(FlintZZ, ncols(A)))
end


