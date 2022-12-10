export solve_non_negative, solve_mixed, solve_ineq


function solve_mixed(as::Type{SubObjectIterator{PointVector{fmpz}}}, A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat)
    ncols(A) == ncols(C) || throw(ArgumentError("solve_mixed(A,b,C,d): A and C must have the same number of columns."))
    nrows(A) == nrows(b) || throw(ArgumentError("solve_mixed(A,b,C,d): A and b must have the same number of rows."))
    nrows(C) == nrows(d) || throw(ArgumentError("solve_mixed(A,b,C,d): C and d must have the same number of rows."))
    ncols(b) == 1 || throw(ArgumentError("solve_mixed(A,b,C,d): b must be a matrix with a single column."))
    ncols(d) == 1 || throw(ArgumentError("solve_mixed(A,b,C,d): d must be a matrix with a single column."))
    P = Polyhedron((-C, vec(-d)), (A, vec(b)))
    return lattice_points(P)
end

function solve_mixed(as::Type{fmpz_mat}, A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat)
    LP = solve_mixed(SubObjectIterator{PointVector{fmpz}}, A, b, C, d)
    return matrix(ZZ, LP)
end


@doc Markdown.doc"""
    solve_mixed(as::Type{T}, A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat) where {T}

Solve $Ax = b$ under $Cx >= d$, assumes a finite solution set.

The output type may be specified in the variable `as`:
- `fmpz_mat` (default) a matrix with integers is returned.
- `SubObjectIterator{PointVector{fmpz}}` an iterator over integer points is returned.

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

julia> typeof(solve_mixed(A, b, C, d))
fmpz_mat

julia> typeof(solve_mixed(fmpz_mat, A, b, C, d))
fmpz_mat

julia> typeof(solve_mixed(SubObjectIterator{PointVector{fmpz}}, A, b, C, d))
SubObjectIterator{PointVector{fmpz}}
```
"""
solve_mixed(as::Type{T}, A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat) where {T} = solve_mixed(T, A, b, C, d)
solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat, d::fmpz_mat) = solve_mixed(fmpz_mat, A, b, C, d)


@doc Markdown.doc"""
    solve_mixed(as::Type{T}, A::fmpz_mat, b::fmpz_mat, C::fmpz_mat) where {T}

Solve $Ax = b$ under $Cx >= 0$, assumes a finite solution set.

The output type may be specified in the variable `as`:
- `fmpz_mat` (default) a matrix with integers is returned.
- `SubObjectIterator{PointVector{fmpz}}` an iterator over integer points is returned.

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

julia> typeof(solve_mixed(A, b, C))
fmpz_mat

julia> typeof(solve_mixed(fmpz_mat, A, b, C))
fmpz_mat

julia> typeof(solve_mixed(SubObjectIterator{PointVector{fmpz}}, A, b, C))
SubObjectIterator{PointVector{fmpz}}
```
"""
solve_mixed(as::Type{T}, A::fmpz_mat, b::fmpz_mat, C::fmpz_mat) where {T} = solve_mixed(T, A, b, C, zero_matrix(FlintZZ, nrows(C), 1))
solve_mixed(A::fmpz_mat, b::fmpz_mat, C::fmpz_mat) = solve_mixed(fmpz_mat, A, b, C, zero_matrix(FlintZZ, nrows(C), 1))



@doc Markdown.doc"""
    solve_ineq(as::Type{T}, A::fmpz_mat, b::fmpz_mat) where {T}

Solve $Ax<=b$, assumes finite set of solutions.

The output type may be specified in the variable `as`:
- `fmpz_mat` (default) a matrix with integers is returned.
- `SubObjectIterator{PointVector{fmpz}}` an iterator over integer points is returned.

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

julia> typeof(solve_ineq(A,b))
fmpz_mat

julia> typeof(solve_ineq(fmpz_mat, A,b))
fmpz_mat

julia> typeof(solve_ineq(SubObjectIterator{PointVector{fmpz}}, A,b))
SubObjectIterator{PointVector{fmpz}}
```
"""
solve_ineq(as::Type{T}, A::fmpz_mat, b::fmpz_mat) where {T} = solve_mixed(T, zero_matrix(FlintZZ, 0, ncols(A)), zero_matrix(FlintZZ,0,1), -A, -b)
solve_ineq(A::fmpz_mat, b::fmpz_mat) = solve_ineq(fmpz_mat, A, b)



@doc Markdown.doc"""
    solve_non_negative(as::Type{T}, A::fmpz_mat, b::fmpz_mat) where {T}

Find all solutions to $Ax = b$, $x>=0$. Assumes a finite set of solutions.

The output type may be specified in the variable `as`:
- `fmpz_mat` (default) a matrix with integers is returned.
- `SubObjectIterator{PointVector{fmpz}}` an iterator over integer points is returned.

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

julia> typeof(solve_non_negative(A,b))
fmpz_mat

julia> typeof(solve_non_negative(fmpz_mat, A,b))
fmpz_mat

julia> typeof(solve_non_negative(SubObjectIterator{PointVector{fmpz}}, A,b))
SubObjectIterator{PointVector{fmpz}}
```
"""
solve_non_negative(as::Type{T}, A::fmpz_mat, b::fmpz_mat) where {T} = solve_mixed(T, A, b, identity_matrix(FlintZZ, ncols(A)))
solve_non_negative(A::fmpz_mat, b::fmpz_mat) = solve_non_negative(fmpz_mat, A, b)


