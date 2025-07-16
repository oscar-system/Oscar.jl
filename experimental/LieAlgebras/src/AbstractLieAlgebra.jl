###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{AbstractLieAlgebraElem{C}}) where {C<:FieldElem} = AbstractLieAlgebra{C}

elem_type(::Type{AbstractLieAlgebra{C}}) where {C<:FieldElem} = AbstractLieAlgebraElem{C}

parent(x::AbstractLieAlgebraElem) = x.parent

coefficient_ring(L::AbstractLieAlgebra{C}) where {C<:FieldElem} = L.R::parent_type(C)

vector_space_dim(L::AbstractLieAlgebra) = L.dim

_struct_consts(L::AbstractLieAlgebra{C}) where {C<:FieldElem} =
  L.struct_consts::Matrix{sparse_row_type(C)}

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
  if has_root_system(L)
    rs = root_system(L)
    if has_root_system_type(rs)
      type, ord = root_system_type_with_ordering(rs)
      print(io, Indent(), "of type ", _root_system_type_string(type))
      if !issorted(ord)
        print(io, " (non-canonical ordering)")
      end
      println(io, Dedent())
    end
  end
  println(io, Indent(), "of dimension ", dim(L), Dedent())
  print(io, "over ", Lowercase(), coefficient_ring(L))
end

function Base.show(io::IO, L::AbstractLieAlgebra)
  @show_name(io, L)
  @show_special(io, L)
  if is_terse(io)
    print(io, "Abstract Lie algebra")
  else
    io = pretty(io)
    print(io, "Abstract Lie algebra")
    if has_root_system(L)
      rs = root_system(L)
      if has_root_system_type(rs)
        type, ord = root_system_type_with_ordering(rs)
        print(io, " of type ", _root_system_type_string(type))
        if !issorted(ord)
          print(io, " (non-canonical ordering)")
        end
      end
    end
    print(terse(io), " over ", Lowercase(), coefficient_ring(L))
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
  vec = sparse_row(coefficient_ring(L))
  for (i, cxi) in enumerate(coefficients(x))
    iszero(cxi) && continue
    for (j, cyj) in enumerate(coefficients(y))
      iszero(cyj) && continue
      vec = addmul!(vec, _struct_consts(L)[i, j], cxi * cyj)
    end
  end
  return L(vec)
end

###############################################################################
#
#   Properties
#
###############################################################################

@attr Bool function is_abelian(L::AbstractLieAlgebra)
  return all(iszero, _struct_consts(L))
end

###############################################################################
#
#   Root system getters
#
###############################################################################

has_root_system(L::AbstractLieAlgebra) = isdefined(L, :root_system)

function root_system(L::AbstractLieAlgebra)
  assure_root_system(L)
  return L.root_system
end

function chevalley_basis(L::AbstractLieAlgebra)
  assure_root_system(L)
  return L.chevalley_basis::NTuple{3,Vector{elem_type(L)}}
end

function set_root_system_and_chevalley_basis!(
  L::AbstractLieAlgebra{C}, R::RootSystem, chev::NTuple{3,Vector{AbstractLieAlgebraElem{C}}}
) where {C<:FieldElem}
  L.root_system = R
  L.chevalley_basis = chev
end

###############################################################################
#
#   Structure constant table and change_base_ring
#
###############################################################################

function change_base_ring(R::Field, L::AbstractLieAlgebra)
  struct_consts = map(e -> change_base_ring(R, e), structure_constant_table(L; copy=false))
  return lie_algebra(R, struct_consts, symbols(L); check=false)
end

function structure_constant_table(L::AbstractLieAlgebra; copy::Bool=true)
  return copy ? deepcopy(_struct_consts(L)) : _struct_consts(L)
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    lie_algebra(R::Field, struct_consts::Matrix{sparse_row_type(R)}, s::Vector{<:VarName}; check::Bool) -> AbstractLieAlgebra{elem_type(R)}

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
  struct_consts::Matrix{<:SRow{C}},
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
  struct_consts2 = Matrix{sparse_row_type(R)}(
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
  @req !isempty(basis) "Basis must not be empty, or provide the Lie algebra as first argument"
  return lie_algebra(parent(basis[1]), basis; check)
end

function lie_algebra(
  L::AbstractLieAlgebra{C}, basis::Vector{AbstractLieAlgebraElem{C}}; check::Bool=true
) where {C<:FieldElem}
  @req all(parent(x) === L for x in basis) "Elements not compatible."
  R = coefficient_ring(L)
  basis_matrix = if length(basis) == 0
    matrix(R, 0, dim(L), C[])
  else
    matrix(R, [coefficients(b) for b in basis])
  end
  struct_consts = Matrix{sparse_row_type(R)}(undef, length(basis), length(basis))
  for (i, bi) in enumerate(basis), (j, bj) in enumerate(basis)
    fl, row = can_solve_with_solution(basis_matrix, _matrix(bi * bj); side=:left)
    @req fl "Not closed under the bracket."
    struct_consts[i, j] = sparse_row(row)
  end
  s = map(AbstractAlgebra.obj_to_string_wrt_times, basis)
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

# Examples
```jldoctest
julia> L = lie_algebra(QQ, root_system(:B, 4))
Abstract Lie algebra
  of type B4
  of dimension 36
over rational field
```
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

  npos = n_positive_roots(rs)
  # set Chevalley basis
  basis_L = basis(L)
  r_plus = basis_L[1:npos]
  r_minus = basis_L[(npos + 1):(2 * npos)]
  h = basis_L[(2 * npos + 1):end]
  chev = (r_plus, r_minus, h)
  set_root_system_and_chevalley_basis!(L, rs, chev)
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

  struct_consts = Matrix{sparse_row_type(R)}(undef, n, n)
  beta_i_plus_beta_j = zero(RootSpaceElem, rs)
  for i in 1:nroots, j in i:nroots
    if i == j
      # [e_βi, e_βi] = 0
      struct_consts[i, j] = sparse_row(R)
      continue
    end
    beta_i = root(rs, i)
    beta_j = root(rs, j)
    beta_i_plus_beta_j = add!(beta_i_plus_beta_j, beta_i, beta_j)
    if iszero(beta_i_plus_beta_j)
      # [e_βi, e_-βi] = h_βi
      struct_consts[i, j] = sparse_row(
        R, collect((nroots + 1):(nroots + nsimp)), _vec(coefficients(coroot(rs, i)));
        sort=false,
      )
    elseif ((fl, k) = is_root_with_index(beta_i_plus_beta_j); fl)
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

    # [e_βj, e_βi] = -[e_βi, e_βj]
    struct_consts[j, i] = -struct_consts[i, j]
  end
  for i in 1:nsimp, j in 1:nroots
    # [h_i, e_βj] = <β_j, α_i> e_βj
    struct_consts[nroots + i, j] = sparse_row(
      R,
      [j],
      # [dot(coefficients(root(rs, j)), view(cm, i, :))], # currently the below is faster
      [only(coefficients(root(rs, j)) * cm[i, :])],
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

# computes the maximum `p` such that `beta - p*alpha` is still a root
# beta is assumed to be a root
function _root_string_length_down(alpha::RootSpaceElem, beta::RootSpaceElem)
  p = 0
  beta_sub_p_alpha = beta - alpha
  while is_root(beta_sub_p_alpha)
    p += 1
    beta_sub_p_alpha = sub!(beta_sub_p_alpha, alpha)
  end
  return p
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
        j -> !is_positive_root(alpha_i + beta_l - simple_root(rs, j)),
        1:(i - 1),
      ) || continue
      p = _root_string_length_down(beta_l, alpha_i) + 1
      N[i, l] = (extraspecial_pair_signs[k - nsimp] ? p : -p)
      N[l, i] = -N[i, l]
    end
  end

  # special pairs
  alpha_i_plus_beta_j = zero(RootSpaceElem, rs)
  for (i, alpha_i) in enumerate(positive_roots(rs))
    for (j, beta_j) in enumerate(positive_roots(rs))
      i < j || continue
      alpha_i_plus_beta_j = add!(alpha_i_plus_beta_j, alpha_i, beta_j)
      is_positive_root(alpha_i_plus_beta_j) || continue
      l = let alpha_i_plus_beta_j = alpha_i_plus_beta_j # avoid closure capture
        findfirst(
          l -> is_positive_root(alpha_i_plus_beta_j - simple_root(rs, l)), 1:nsimp
        )::Int
      end
      l == i && continue # already extraspecial
      fl, l_comp = is_positive_root_with_index(alpha_i_plus_beta_j - simple_root(rs, l))
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
      p = _root_string_length_down(beta_j, alpha_i) + 1
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

# Examples
```jldoctest
julia> L = lie_algebra(QQ, :C, 4)
Abstract Lie algebra
  of type C4
  of dimension 36
over rational field
```
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

  struct_consts = Matrix{sparse_row_type(R)}(undef, n, n)
  for i in axes(struct_consts, 1), j in axes(struct_consts, 2)
    struct_consts[i, j] = sparse_row(R)
  end

  L = lie_algebra(R, struct_consts, s)
  set_attribute!(L, :is_abelian => true)
  return L
end
