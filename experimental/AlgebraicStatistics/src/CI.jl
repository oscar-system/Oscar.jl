# -*- Conditional independence statements -*-

export CIStmt, ci_stmt, @CI_str, ci_statements, make_elementary

struct CIStmt
  I::Vector{<:VarName}
  J::Vector{<:VarName}
  K::Vector{<:VarName}
end

@doc raw"""
    ci_stmt(I::Vector{<:VarName}, J::Vector{<:VarName}, K::Vector{<:VarName}; symmetric=true, semigraphoid=true)

A conditional independence statement asserting that `I` is independent
of `J` given `K`. These parameters are lists of names of random variables.
The sets `I` and `J` must be disjoint as this package cannot yet deal
with functional dependencies.

If `symmetric` is `true`, CI statements are assumed to be symmetric in
their `I` and `J` components. The constructor then reorders the arguments
to make the `I` field lexicographically smaller than the `J` to ensure
that comparisons and hashing respect the symmetry.

If `semigraphoid` is set to `true`, the constructor also removes elements
in the intersection of `I` and `K` from `I` (and symetrically removes the
intersection of `J` and `K` from `J`).

As all three fields are sets, each of them may be deduplicated and sorted
to ensure consistent comparison and hashing.

## Examples

``` jldoctest
julia> ci_stmt(["A"], ["B"], ["X"])
[A _||_ B | X]

julia> ci_stmt(["1"], ["2", "3"], ["4", "5"])
[1 _||_ {2, 3} | {4, 5}]
```
"""
function ci_stmt(I::Vector{<:VarName}, J::Vector{<:VarName}, K::Vector{<:VarName}; symmetric=true, semigraphoid=true)
  if length(intersect(I, J)) > 0
    error("Functional dependence statements are not yet implemented")
  end
  if symmetric && I > J
    I, J = J, I
  end
  if semigraphoid
    I = setdiff(I, K)
    J = setdiff(J, K)
  end
  CIStmt(sort(unique(I)), sort(unique(J)), sort(unique(K)))
end

@doc raw"""
    Base.:(==)(lhs::CIStmt, rhs::CIStmt)

Compares `CIStmt`s for identity in all their three fields.
"""
Base.:(==)(lhs::CIStmt, rhs::CIStmt) =
  lhs.I == rhs.I && lhs.J == rhs.J && lhs.K == rhs.K

@doc raw"""
    Base.hash(stmt:;CIStmt, h::UInt)

Computes the hash of a `CIStmt`.
"""
Base.hash(stmt::CIStmt, h::UInt) =
  foldr(hash, stmt.I, stmt.J, stmt.K; init=hash(CIStmt, h))

@doc raw"""
    CI"I...,J...|K..."

A literal syntax for denoting CI statements is provided for cases in which
all variable names consist of a single character. If `I` and `J` only consist
of a single element, then even the comma may be omitted. Once the three sets
are extracted, `ci_stmt` is called.

## Examples

``` jldoctest
julia> CI"AB|X"
[A _||_ B | X]

julia> CI"1,23|5424"
[1 _||_ 3 | {2, 4, 5}]
```
"""
macro CI_str(str)
  # General syntax "12,34|567"
  m = match(r"^(.+),(.+)[|](.*)$", str)
  # Short syntax for elementary statements "12|345"
  if m == nothing
    m = match(r"^(.)(.)[|](.*)$", str)
  end
  # Give up
  if m == nothing
    throw(ArgumentError(str * " is not a CI statement"))
  end
  parse_arg(s) = unique([string(c) for c in s])
  I, J, K = map(s -> parse_arg(s), m)
  return ci_stmt(I, J, K)
end

function Base.show(io::IO, stmt::CIStmt)
  fmt(K) = length(K) == 1 ? string(K[1]) : "{" * join([string(x) for x in K], ", ") * "}"
  if Oscar.is_unicode_allowed()
    print(io, "[$(fmt(stmt.I)) ⫫ $(fmt(stmt.J)) | $(fmt(stmt.K))]")
  else
    print(io, "[$(fmt(stmt.I)) _||_ $(fmt(stmt.J)) | $(fmt(stmt.K))]")
  end
end

@doc raw"""
    ci_statements(random_variables::Vector{<:VarName})

Return a list of all elementary CI statements over a given set of
variable names. A `CIStmt(I, J, K)` is elementary if both `I` and
`J` have only one element.

As a consequence of the semigraphoid properties, these statements
are enough to describe the entire CI structure of a probability
distribution.

## Examples

``` jldoctest
julia> ci_statements(["A", "B", "X", "Y"])
24-element Vector{CIStmt}:
 [A _||_ Y | {}]
 [A _||_ Y | B]
 [A _||_ Y | X]
 [A _||_ Y | {B, X}]
 [B _||_ Y | {}]
 [B _||_ Y | A]
 [B _||_ Y | X]
 [B _||_ Y | {A, X}]
 [X _||_ Y | {}]
 [X _||_ Y | A]
 ⋮
 [A _||_ X | {B, Y}]
 [B _||_ X | {}]
 [B _||_ X | A]
 [B _||_ X | Y]
 [B _||_ X | {A, Y}]
 [A _||_ B | {}]
 [A _||_ B | X]
 [A _||_ B | Y]
 [A _||_ B | {X, Y}]
```
"""
function ci_statements(random_variables::Vector{<:VarName})
  N = collect(1:length(random_variables))
  stmts = Vector{CIStmt}()
  for ij in subsets(N, 2)
    M = setdiff(N, ij)
    for l in 0:length(M)
      for L in subsets(M, l)
        push!(stmts, CIStmt(
          [random_variables[ij[1]]],
          [random_variables[ij[2]]],
           random_variables[L]
        ))
      end
    end
  end
  return stmts
end

@doc raw"""
    make_elementary(stmt::CIStmt; semigaussoid=false)

Convert a `CIStmt` into an equivalent list of `CIStmt`s all of which
are elementary. The default operation assumes the semigraphoid axioms
and converts ``[I \mathrel{⫫} J \mid K]`` into the list consisting of
``[i \mathrel{⫫} j \mid L]`` for all ``i \in I``, ``j \in J`` and ``L``
in the interval ``K \subseteq L \subseteq (I \cup J \cup K) \setminus \{i,j\}``.

If `semigaussoid` is true, the stronger semigaussoid axioms are
assumed and `L` in the above procedure does not range in the interval
above `K` but is always fixed to `K`. Semigaussoids are also known as
*compositional graphoids*.

## Examples

``` jldoctest
julia> make_elementary(CI"12,34|56")
16-element Vector{CIStmt}:
 [1 _||_ 3 | {5, 6}]
 [1 _||_ 3 | {5, 6, 2}]
 [1 _||_ 3 | {5, 6, 4}]
 [1 _||_ 3 | {5, 6, 2, 4}]
 [1 _||_ 4 | {5, 6}]
 [1 _||_ 4 | {5, 6, 2}]
 [1 _||_ 4 | {5, 6, 3}]
 [1 _||_ 4 | {5, 6, 2, 3}]
 [2 _||_ 3 | {5, 6}]
 [2 _||_ 3 | {5, 6, 1}]
 [2 _||_ 3 | {5, 6, 4}]
 [2 _||_ 3 | {5, 6, 1, 4}]
 [2 _||_ 4 | {5, 6}]
 [2 _||_ 4 | {5, 6, 1}]
 [2 _||_ 4 | {5, 6, 3}]
 [2 _||_ 4 | {5, 6, 1, 3}]

julia> make_elementary(CI"12,34|56"; semigaussoid=true)
4-element Vector{CIStmt}:
 [1 _||_ 3 | {5, 6}]
 [1 _||_ 4 | {5, 6}]
 [2 _||_ 3 | {5, 6}]
 [2 _||_ 4 | {5, 6}]
```
"""
function make_elementary(stmt::CIStmt; semigaussoid=false)
  N = union(stmt.I, stmt.J, stmt.K)
  elts = Vector{CIStmt}()
  for i in stmt.I
    for j in stmt.J
      # If we may suppose the semigaussoid axioms, the system becomes simpler.
      M = semigaussoid ? typeof(stmt.K)([]) : setdiff(N, [i, j, stmt.K...])
      for l in 0:length(M)
        for L in subsets(M, l)
          push!(elts, CIStmt([i], [j], union(stmt.K, L)))
        end
      end
    end
  end
  return elts
end
