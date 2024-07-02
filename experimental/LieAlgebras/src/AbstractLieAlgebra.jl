@attributes mutable struct AbstractLieAlgebra{C<:FieldElem} <: LieAlgebra{C}
  R::Field
  dim::Int
  struct_consts::Matrix{SRow{C}}
  s::Vector{Symbol}

  # only set if known
  root_system::RootSystem

  function AbstractLieAlgebra{C}(
    R::Field,
    struct_consts::Matrix{SRow{C}},
    s::Vector{Symbol};
    check::Bool=true,
  ) where {C<:FieldElem}
    (n1, n2) = size(struct_consts)
    @req n1 == n2 "Invalid structure constants dimensions."
    dimL = n1
    @req length(s) == dimL "Invalid number of basis element names."
    if check
      @req all(
        r -> all(e -> parent(last(e)) === R, r), struct_consts
      ) "Invalid structure constants."
      @req all(iszero, struct_consts[i, i] for i in 1:dimL) "Not anti-symmetric."
      @req all(
        iszero, struct_consts[i, j] + struct_consts[j, i] for i in 1:dimL, j in 1:dimL
      ) "Not anti-symmetric."
      @req all(
        iszero,
        sum(
          struct_consts[i, j][k] * struct_consts[k, l] +
          struct_consts[j, l][k] * struct_consts[k, i] +
          struct_consts[l, i][k] * struct_consts[k, j] for k in 1:dimL
        ) for i in 1:dimL, j in 1:dimL, l in 1:dimL
      ) "Jacobi identity does not hold."
    end
    return new{C}(R, dimL, struct_consts, s)
  end
end

struct AbstractLieAlgebraElem{C<:FieldElem} <: LieAlgebraElem{C}
  parent::AbstractLieAlgebra{C}
  mat::MatElem{C}
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{AbstractLieAlgebraElem{C}}) where {C<:FieldElem} = AbstractLieAlgebra{C}

elem_type(::Type{AbstractLieAlgebra{C}}) where {C<:FieldElem} = AbstractLieAlgebraElem{C}

parent(x::AbstractLieAlgebraElem) = x.parent

coefficient_ring(L::AbstractLieAlgebra{C}) where {C<:FieldElem} = L.R::parent_type(C)

dim(L::AbstractLieAlgebra) = L.dim

###############################################################################
#
#   Root system getters
#
###############################################################################

has_root_system(L::LieAlgebra) = isdefined(L, :root_system)

function root_system(L::LieAlgebra)
  @req has_root_system(L) "No root system known."
  return L.root_system
end

@doc raw"""
    chevalley_basis(L::AbstractLieAlgebra{C}) -> NTuple{3,Vector{AbstractLieAlgebraElem{C}}}

Return the Chevalley basis of the Lie algebra `L` in three vectors, stating first the positive root vectors, 
then the negative root vectors, and finally the basis of the Cartan subalgebra. The order of root vectors corresponds
to the order of the roots in the root system.
"""
function chevalley_basis(L::AbstractLieAlgebra)
  @req has_root_system(L) "No root system known."
  # TODO: once there is root system detection, this function needs to be updated to indeed return the Chevalley basis

  npos = n_positive_roots(root_system(L))
  b = basis(L)
  # root vectors
  r_plus = b[1:npos]
  r_minus = b[(npos + 1):(2 * npos)]
  # basis for cartan algebra
  h = b[(2 * npos + 1):dim(L)]
  return (r_plus, r_minus, h)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, mime::MIME"text/plain", L::AbstractLieAlgebra)
  @show_name(io, L)
  @show_special(io, mime, L)
  io = pretty(io)
  println(io, "Abstract Lie algebra")
  println(io, Indent(), "of dimension $(dim(L))", Dedent())
  print(io, "over ")
  print(io, Lowercase(), coefficient_ring(L))
end

function Base.show(io::IO, L::AbstractLieAlgebra)
  @show_name(io, L)
  @show_special(io, L)
  if is_terse(io)
    print(io, "Abstract Lie algebra")
  else
    io = pretty(io)
    print(io, "Abstract Lie algebra over ", Lowercase())
    print(terse(io), coefficient_ring(L))
  end
end

function symbols(L::AbstractLieAlgebra)
  return L.s
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

# no special ones

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function bracket(
  x::AbstractLieAlgebraElem{C}, y::AbstractLieAlgebraElem{C}
) where {C<:FieldElem}
  check_parent(x, y)
  L = parent(x)
  mat = sum(
    cxi * cyj * L.struct_consts[i, j] for (i, cxi) in enumerate(coefficients(x)),
    (j, cyj) in enumerate(coefficients(y));
    init=sparse_row(coefficient_ring(L)),
  )
  return L(mat)
end

###############################################################################
#
#   Properties
#
###############################################################################

function is_abelian(L::AbstractLieAlgebra)
  return all(e -> iszero(length(e)), L.struct_consts)
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    lie_algebra(R::Field, struct_consts::Matrix{SRow{elem_type(R)}}, s::Vector{<:VarName}; check::Bool) -> AbstractLieAlgebra{elem_type(R)}

Construct the Lie algebra over the field `R` with structure constants `struct_consts`
and with basis element names `s`.

The Lie bracket on the newly constructed Lie algebra `L` is determined by the structure
constants in `struct_consts` as follows: let $x_i$ denote the $i$-th standard basis vector
of `L`. Then the entry `struct_consts[i,j][k]` is a scalar $a_{i,j,k}$
such that $[x_i, x_j] = \sum_k a_{i,j,k} x_k$.

* `s`: A vector of basis element names. This is 
  `[Symbol("x_$i") for i in 1:size(struct_consts, 1)]` by default.
* `check`: If `true`, check that the structure constants are anti-symmetric and
  satisfy the Jacobi identity. This is `true` by default.
"""
function lie_algebra(
  R::Field,
  struct_consts::Matrix{SRow{C}},
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:size(struct_consts, 1)];
  check::Bool=true,
) where {C<:FieldElem}
  @req C == elem_type(R) "Invalid coefficient type."
  return AbstractLieAlgebra{elem_type(R)}(R, struct_consts, Symbol.(s); check)
end

@doc raw"""
    lie_algebra(R::Field, struct_consts::Array{elem_type(R),3}, s::Vector{<:VarName}; check::Bool) -> AbstractLieAlgebra{elem_type(R)}

Construct the Lie algebra over the field `R` with structure constants `struct_consts`
and with basis element names `s`.

The Lie bracket on the newly constructed Lie algebra `L` is determined by the structure
constants in `struct_consts` as follows: let $x_i$ denote the $i$-th standard basis vector
of `L`. Then the entry `struct_consts[i,j,k]` is a scalar $a_{i,j,k}$
such that $[x_i, x_j] = \sum_k a_{i,j,k} x_k$.

* `s`: A vector of basis element names. This is
  `[Symbol("x_$i") for i in 1:size(struct_consts, 1)]` by default.
* `check`: If `true`, check that the structure constants are anti-symmetric and
  satisfy the Jacobi identity. This is `true` by default.

# Examples
```jldoctest
julia> struct_consts = zeros(QQ, 3, 3, 3);

julia> struct_consts[1, 2, 3] = QQ(1);

julia> struct_consts[2, 1, 3] = QQ(-1);

julia> struct_consts[3, 1, 1] = QQ(2);

julia> struct_consts[1, 3, 1] = QQ(-2);

julia> struct_consts[3, 2, 2] = QQ(-2);

julia> struct_consts[2, 3, 2] = QQ(2);

julia> sl2 = lie_algebra(QQ, struct_consts, ["e", "f", "h"])
Abstract Lie algebra
  of dimension 3
over rational field

julia> e, f, h = basis(sl2)
3-element Vector{AbstractLieAlgebraElem{QQFieldElem}}:
 e
 f
 h

julia> e * f
h

julia> h * e
2*e

julia> h * f
-2*f
```
"""
function lie_algebra(
  R::Field,
  struct_consts::Array{C,3},
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:size(struct_consts, 1)];
  check::Bool=true,
) where {C<:FieldElem}
  @req C == elem_type(R) "Invalid coefficient type."
  struct_consts2 = Matrix{SRow{elem_type(R)}}(
    undef, size(struct_consts, 1), size(struct_consts, 2)
  )
  for i in axes(struct_consts, 1), j in axes(struct_consts, 2)
    struct_consts2[i, j] = sparse_row(
      R, collect(axes(struct_consts, 3)), struct_consts[i, j, :]
    )
  end

  return AbstractLieAlgebra{elem_type(R)}(R, struct_consts2, Symbol.(s); check)
end

function lie_algebra(
  basis::Vector{AbstractLieAlgebraElem{C}}; check::Bool=true
) where {C<:FieldElem}
  parent_L = parent(basis[1])
  @req all(parent(x) === parent_L for x in basis) "Elements not compatible."
  R = coefficient_ring(parent_L)
  basis_matrix = if length(basis) == 0
    matrix(R, 0, dim(L), C[])
  else
    matrix(R, [coefficients(b) for b in basis])
  end
  struct_consts = Matrix{SRow{elem_type(R)}}(undef, length(basis), length(basis))
  for (i, bi) in enumerate(basis), (j, bj) in enumerate(basis)
    fl, row = can_solve_with_solution(basis_matrix, _matrix(bi * bj); side=:left)
    @req fl "Not closed under the bracket."
    struct_consts[i, j] = sparse_row(row)
  end
  s = map(AbstractAlgebra.obj_to_string, basis)
  return lie_algebra(R, struct_consts, s; check)
end

@doc raw"""
    lie_algebra(R::Field, rs::RootSystem) -> AbstractLieAlgebra{elem_type(R)}

Construct a simple Lie algebra over the field `R` with the root system `rs`.
The internally used basis of this Lie algebra is the Chevalley basis.

The experienced user may supply a boolean vector of length `n_positive_roots(rs) - n_simple_roots(rs)`
via the kwarg `extraspecial_pair_signs::Vector{Bool}` to specify the concrete Lie algebra to be constructed.
If $(\alpha,\beta)$ is the extraspecial pair for the non-simple root `root(rs, i)`,
then $\varepsilon_{\alpha,\beta} = 1$ iff `extraspecial_pair_signs[i - n_simple_roots(rs)] = true`.
For the used notation and the definition of extraspecial pairs, see [CMT04](@cite).
"""
function lie_algebra(
  R::Field,
  rs::RootSystem;
  extraspecial_pair_signs=fill(true, n_positive_roots(rs) - n_simple_roots(rs)),
)
  struct_consts = _struct_consts(R, rs, extraspecial_pair_signs)

  s = [
    [Symbol("x_$i") for i in 1:n_positive_roots(rs)]
    [Symbol("y_$i") for i in 1:n_positive_roots(rs)]
    [Symbol("h_$i") for i in 1:n_simple_roots(rs)]
  ]

  L = lie_algebra(R, struct_consts, s; check=false)
  L.root_system = rs
  return L
end

function _struct_consts(R::Field, rs::RootSystem, extraspecial_pair_signs)
  cm = cartan_matrix(rs)
  @req is_cartan_matrix(cm; generalized=false) "The root system does not induce a finite dimensional Lie algebra."

  nroots = n_roots(rs)
  npos = n_positive_roots(rs)
  nsimp = n_simple_roots(rs)

  n = 2 * npos + nsimp

  N = _N_matrix(rs, extraspecial_pair_signs)

  struct_consts = Matrix{SRow{elem_type(R)}}(undef, n, n)
  for i in 1:nroots, j in i:nroots
    if i == j
      # [e_βi, e_βi] = 0
      struct_consts[i, j] = sparse_row(R)
      continue
    end
    beta_i = root(rs, i)
    beta_j = root(rs, j)
    if iszero(beta_i + beta_j)
      # [e_βi, e_-βi] = h_βi
      struct_consts[i, j] = sparse_row(
        R, collect((nroots + 1):(nroots + nsimp)), _vec(coefficients(coroot(rs, i)))
      )
    elseif ((fl, k) = is_root_with_index(beta_i + beta_j); fl)
      # complicated case
      if i <= npos
        struct_consts[i, j] = sparse_row(R, [k], [N[i, j]])
      elseif j <= npos
        struct_consts[i, j] = sparse_row(R, [k], [-N[i - npos, j + npos]])
      else
        struct_consts[i, j] = sparse_row(R, [k], [-N[i - npos, j - npos]])
      end
    else
      # [e_βi, e_βj] = 0
      struct_consts[i, j] = sparse_row(R)
    end

    # # [e_βj, e_βi] = -[e_βi, e_βj]
    struct_consts[j, i] = -struct_consts[i, j]
  end
  for i in 1:nsimp, j in 1:nroots
    # [h_i, e_βj] = <β_j, α_i> e_βj
    struct_consts[nroots + i, j] = sparse_row(
      R,
      [j],
      [dot(coefficients(root(rs, j)), cm[i, :])],
    )
    # [e_βj, h_i] = -[h_i, e_βj]
    struct_consts[j, nroots + i] = -struct_consts[nroots + i, j]
  end
  for i in 1:nsimp, j in 1:nsimp
    # [h_i, h_j] = 0
    struct_consts[nroots + i, nroots + j] = sparse_row(R)
  end

  return struct_consts
end

function _N_matrix(rs::RootSystem, extraspecial_pair_signs::Vector{Bool})
  # computes the matrix N_αβ from CTM04 Ch. 3 indexed by root indices
  nroots = n_roots(rs)
  npos = n_positive_roots(rs)
  nsimp = n_simple_roots(rs)
  @req length(extraspecial_pair_signs) == npos - nsimp "Invalid extraspecial pair sign vector length."

  N = zeros(Int, npos, nroots)

  # extraspecial pairs
  for (i, alpha_i) in enumerate(simple_roots(rs))
    for (l, beta_l) in enumerate(positive_roots(rs))
      fl, k = is_positive_root_with_index(alpha_i + beta_l)
      fl || continue
      all(
        j -> !is_positive_root_with_index(alpha_i + beta_l - simple_root(rs, j))[1],
        1:(i - 1),
      ) || continue
      p = 0
      while is_root_with_index(beta_l - p * alpha_i)[1]
        p += 1
      end
      N[i, l] = (extraspecial_pair_signs[k - nsimp] ? 1 : -1) * p
      N[l, i] = -N[i, l]
    end
  end

  # special pairs
  for (i, alpha_i) in enumerate(positive_roots(rs))
    for (j, beta_j) in enumerate(positive_roots(rs))
      i < j || continue
      fl = is_positive_root_with_index(alpha_i + beta_j)[1]
      fl || continue
      l = findfirst(
        l -> is_positive_root_with_index(alpha_i + beta_j - simple_root(rs, l))[1], 1:nsimp
      )
      l == i && continue # already extraspecial
      fl, l_comp = is_positive_root_with_index(alpha_i + beta_j - simple_root(rs, l))
      @assert fl
      t1 = 0
      t2 = 0
      if ((fl, m) = is_positive_root_with_index(beta_j - simple_root(rs, l)); fl)
        root_m = positive_root(rs, m)
        t1 = N[l, m] * N[i, m] * dot(root_m, root_m)//dot(beta_j, beta_j)
      end
      if ((fl, m) = is_positive_root_with_index(alpha_i - simple_root(rs, l)); fl)
        root_m = positive_root(rs, m)
        t2 = N[l, m] * N[j, m] * dot(root_m, root_m)//dot(alpha_i, alpha_i)
      end
      @assert t1 - t2 != 0
      p = 0
      while is_root_with_index(beta_j - p * alpha_i)[1]
        p += 1
      end
      N[i, j] = Int(sign(t1 - t2) * sign(N[l, l_comp]) * p) # typo in CMT04
      N[j, i] = -N[i, j]
    end
  end

  # rest
  for (i, alpha_i) in enumerate(positive_roots(rs))
    for (j, beta_j) in enumerate(positive_roots(rs))
      if ((fl, k) = is_positive_root_with_index(alpha_i - beta_j); fl)
        root_k = positive_root(rs, k)
        N[i, npos + j] = Int(N[k, j] * dot(root_k, root_k)//dot(alpha_i, alpha_i))
      end
      if ((fl, k) = is_positive_root_with_index(beta_j - alpha_i); fl)
        root_k = positive_root(rs, k)
        N[i, npos + j] = Int(N[k, i] * dot(root_k, root_k)//dot(beta_j, beta_j))
      end
    end
  end
  return N
end

@doc raw"""
    lie_algebra(R::Field, fam::Symbol, rk::Int) -> AbstractLieAlgebra{elem_type(R)}

Construct a simple Lie algebra over the field `R` with Dynkin type given by `fam` and `rk`.
See `cartan_matrix(fam::Symbol, rk::Int)` for allowed combinations.
The internally used basis of this Lie algebra is the Chevalley basis.
"""
function lie_algebra(R::Field, S::Symbol, n::Int)
  rs = root_system(S, n)
  L = lie_algebra(R, rs)

  characteristic(R) == 0 && set_attribute!(L, :is_simple, true)
  return L
end

function abelian_lie_algebra(::Type{T}, R::Field, n::Int) where {T<:AbstractLieAlgebra}
  @req n >= 0 "Dimension must be non-negative."
  basis = [(b = zero_matrix(R, n, n); b[i, i] = 1; b) for i in 1:n]
  s = ["x_$(i)" for i in 1:n]
  L = lie_algebra(R, n, basis, s; check=false)

  struct_consts = Matrix{SRow{elem_type(R)}}(undef, n, n)
  for i in axes(struct_consts, 1), j in axes(struct_consts, 2)
    struct_consts[i, j] = sparse_row(R)
  end

  L = lie_algebra(R, struct_consts, s)
  set_attribute!(L, :is_abelian => true)
  return L
end
