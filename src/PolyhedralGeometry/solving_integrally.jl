function solve_mixed(
  as::Type{SubObjectIterator{PointVector{ZZRingElem}}},
  A::ZZMatrix,
  b::ZZMatrix,
  C::ZZMatrix,
  d::ZZMatrix;
  permit_unbounded=false,
  check::Bool=true,
)
  @req ncols(A) == ncols(C) "solve_mixed(A,b,C,d): A and C must have the same number of columns."
  @req nrows(A) == nrows(b) "solve_mixed(A,b,C,d): A and b must have the same number of rows."
  @req nrows(C) == nrows(d) "solve_mixed(A,b,C,d): C and d must have the same number of rows."
  @req ncols(b) == 1 "solve_mixed(A,b,C,d): b must be a matrix with a single column."
  @req ncols(d) == 1 "solve_mixed(A,b,C,d): d must be a matrix with a single column."

  permit_unbounded &&
    return (pm_object(polyhedron((-C, _vec(-d)), (A, _vec(b)))).LATTICE_POINTS_GENERATORS)[1][
      :, 2:end
    ]

  P = polyhedron((-C, _vec(-d)), (A, _vec(b)); is_bounded=check ? nothing : true)
  return lattice_points(P; check)
end

function solve_mixed(
  as::Type{ZZMatrix},
  A::ZZMatrix,
  b::ZZMatrix,
  C::ZZMatrix,
  d::ZZMatrix;
  permit_unbounded=false,
  check::Bool=true,
)
  LP = solve_mixed(
    SubObjectIterator{PointVector{ZZRingElem}}, A, b, C, d; permit_unbounded, check
  )
  return matrix(ZZ, LP)
end

@doc raw"""
    solve_mixed(as::Type{T}, A::ZZMatrix, b::ZZMatrix, C::ZZMatrix, d::ZZMatrix) where {T}

Solve $Ax = b$ under $Cx >= d$, assumes a finite solution set.

The output type may be specified in the variable `as`:
- `ZZMatrix` (default) a matrix with integers is returned. The solutions are
  the (transposed) rows of the output.
- `SubObjectIterator{PointVector{ZZRingElem}}` an iterator over integer points
  is returned.

# Examples
Find all $(x_1, x_2)\in\mathbb{Z}^2$ such that $x_1+x_2=7$, $x_1\ge 2$, and $x_2\ge 3$.
Note that the output can be permuted, hence we sort it.
```jldoctest
julia> A = ZZMatrix([1 1]);

julia> b = zero_matrix(ZZ, 1,1); b[1,1]=7;

julia> C = ZZMatrix([1 0; 0 1]);

julia> d = zero_matrix(ZZ,2,1); d[1,1]=2; d[2,1]=3;

julia> sortslices(Matrix{BigInt}(solve_mixed(A, b, C, d)), dims=1)
3×2 Matrix{BigInt}:
 2  5
 3  4
 4  3

julia> typeof(solve_mixed(A, b, C, d))
ZZMatrix

julia> typeof(solve_mixed(ZZMatrix, A, b, C, d))
ZZMatrix

julia> it = solve_mixed(SubObjectIterator{PointVector{ZZRingElem}}, A, b, C);

julia> typeof(it)
SubObjectIterator{PointVector{ZZRingElem}}

julia> for x in it
       print(A*x," ")
       end
[7] [7] [7] [7] [7] [7] [7] [7] 
```
"""
solve_mixed(
  as::Type{T}, A::ZZMatrix, b::ZZMatrix, C::ZZMatrix, d::ZZMatrix; permit_unbounded=false
) where {T} = solve_mixed(T, A, b, C, d; permit_unbounded)
solve_mixed(A::ZZMatrix, b::ZZMatrix, C::ZZMatrix, d::ZZMatrix; permit_unbounded=false) =
  solve_mixed(
    ZZMatrix, A, b, C, d; permit_unbounded
  )

@doc raw"""
    solve_mixed(as::Type{T}, A::ZZMatrix, b::ZZMatrix, C::ZZMatrix) where {T}

Solve $Ax = b$ under $Cx >= 0$, assumes a finite solution set.

The output type may be specified in the variable `as`:
- `ZZMatrix` (default) a matrix with integers is returned. The solutions are
  the (transposed) rows of the output.
- `SubObjectIterator{PointVector{ZZRingElem}}` an iterator over integer points
  is returned.

# Examples
Find all $(x_1, x_2)\in\mathbb{Z}^2_{\ge 0}$ such that $x_1+x_2=3$.
Note that the output can be permuted, hence we sort it.
```jldoctest
julia> A = ZZMatrix([1 1]);

julia> b = zero_matrix(ZZ, 1,1); b[1,1]=3;

julia> C = ZZMatrix([1 0; 0 1]);

julia> sortslices(Matrix{BigInt}(solve_mixed(A, b, C)), dims=1)
4×2 Matrix{BigInt}:
 0  3
 1  2
 2  1
 3  0

julia> typeof(solve_mixed(A, b, C))
ZZMatrix

julia> typeof(solve_mixed(ZZMatrix, A, b, C))
ZZMatrix

julia> it = solve_mixed(SubObjectIterator{PointVector{ZZRingElem}}, A, b, C);

julia> typeof(it)
SubObjectIterator{PointVector{ZZRingElem}}

julia> for x in it
       print(A*x," ")
       end
[3] [3] [3] [3] 
```
"""
solve_mixed(
  as::Type{T}, A::ZZMatrix, b::ZZMatrix, C::ZZMatrix; permit_unbounded=false,
  check::Bool=true,
) where {T} = solve_mixed(T, A, b, C, zero_matrix(ZZ, nrows(C), 1); permit_unbounded, check)

solve_mixed(
  A::ZZMatrix, b::ZZMatrix, C::ZZMatrix; permit_unbounded=false, check::Bool=true
) = solve_mixed(
  ZZMatrix, A, b, C, zero_matrix(ZZ, nrows(C), 1); permit_unbounded, check
)

@doc raw"""
    solve_ineq(as::Type{T}, A::ZZMatrix, b::ZZMatrix) where {T}

Solve $Ax<=b$, assumes finite set of solutions.

The output type may be specified in the variable `as`:
- `ZZMatrix` (default) a matrix with integers is returned.
- `SubObjectIterator{PointVector{ZZRingElem}}` an iterator over integer points is returned.

# Examples
The following gives the vertices of the square.
The solutions are the rows of the output.
Note that the output can be permuted, hence we sort it.
```jldoctest
julia> A = ZZMatrix([1 0; 0 1; -1 0; 0 -1]);

julia> b = zero_matrix(ZZ, 4,1); b[1,1]=1; b[2,1]=1; b[3,1]=0; b[4,1]=0;

julia> sortslices(Matrix{BigInt}(solve_ineq(A, b)), dims=1)
4×2 Matrix{BigInt}:
 0  0
 0  1
 1  0
 1  1

julia> typeof(solve_ineq(A,b))
ZZMatrix

julia> typeof(solve_ineq(ZZMatrix, A,b))
ZZMatrix

julia> typeof(solve_ineq(SubObjectIterator{PointVector{ZZRingElem}}, A,b))
SubObjectIterator{PointVector{ZZRingElem}}
```
"""
solve_ineq(as::Type{T}, A::ZZMatrix, b::ZZMatrix; permit_unbounded=false) where {T} =
  solve_mixed(
    T,
    zero_matrix(ZZ, 0, ncols(A)),
    zero_matrix(ZZ, 0, 1),
    -A,
    -b;
    permit_unbounded,
  )
solve_ineq(A::ZZMatrix, b::ZZMatrix; permit_unbounded=false) = solve_ineq(
  ZZMatrix, A, b; permit_unbounded
)

@doc raw"""
    solve_non_negative(as::Type{T}, A::ZZMatrix, b::ZZMatrix) where {T}

Find all solutions to $Ax = b$, $x>=0$. Assumes a finite set of solutions.

The output type may be specified in the variable `as`:
- `ZZMatrix` (default) a matrix with integers is returned.
- `SubObjectIterator{PointVector{ZZRingElem}}` an iterator over integer points is returned.

# Examples
Find all $(x_1, x_2)\in\mathbb{Z}^2_{\ge 0}$ such that $x_1+x_2=3$.
The solutions are the rows of the output.
Note that the output can be permuted, hence we sort it.
```jldoctest
julia> A = ZZMatrix([1 1]);

julia> b = zero_matrix(ZZ, 1,1); b[1,1]=3;

julia> sortslices(Matrix{BigInt}(solve_non_negative(A, b)), dims=1)
4×2 Matrix{BigInt}:
 0  3
 1  2
 2  1
 3  0

julia> typeof(solve_non_negative(A,b))
ZZMatrix

julia> typeof(solve_non_negative(ZZMatrix, A,b))
ZZMatrix

julia> typeof(solve_non_negative(SubObjectIterator{PointVector{ZZRingElem}}, A,b))
SubObjectIterator{PointVector{ZZRingElem}}
```
"""
solve_non_negative(
  as::Type{T}, A::ZZMatrix, b::ZZMatrix; permit_unbounded=false
) where {T} = solve_mixed(T, A, b, identity_matrix(ZZ, ncols(A)); permit_unbounded)
solve_non_negative(A::ZZMatrix, b::ZZMatrix; permit_unbounded=false) = solve_non_negative(
  ZZMatrix, A, b; permit_unbounded
)
