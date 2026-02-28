#################################################
# Linear Algebraic Groups
#################################################

@doc raw"""
    linear_algebraic_group(rs::RootSystem, k::Field) -> LinearAlgebraicGroup

Construct the linear algebraic group of the given type.

Only type ``A_n`` is implemented so far.

# Examples
```jldoctest
julia> F, _  = finite_field(5);

julia> rs = root_system(:A, 3);

julia> LAG = linear_algebraic_group(rs, F)
Linear algebraic group of type A3
  over prime field of characteristic 5
```
"""
function linear_algebraic_group(rs::RootSystem, k::Field)
  if !is_finite(k)
    error("Currently only finite fields are supported.")
  end
  rst = root_system_type(rs)
  if length(rst) == 1 && rst[1][1] == :A
    n = rank(rs)
    G = special_linear_group(n + 1, k)
    LAG = LinearAlgebraicGroup(rs, G, k)
  else
    error("Only type A is implemented so far.")
  end
  return LAG
end

@doc raw"""
    linear_algebraic_group(type::Symbol, n::Int, k::Field) -> LinearAlgebraicGroup

Construct the linear algebraic group of the given type.

Only type ``A_n`` is implemented so far.

# Examples
```jldoctest
julia> F, _  = finite_field(5);

julia> LAG = linear_algebraic_group(:A, 3, F)
Linear algebraic group of type A3
  over prime field of characteristic 5
```
"""
function linear_algebraic_group(type::Symbol, n::Int, k::Field)
  return linear_algebraic_group(root_system(type, n), k)
end

function root_system(LAG::LinearAlgebraicGroup)
  return LAG.rs
end

function _underlying_matrix_group(LAG::LinearAlgebraicGroup{C}) where {C<:FieldElem}
  return LAG.G::matrix_group_type(C)
end

function base_ring(LAG::LinearAlgebraicGroup{C}) where {C<:FieldElem}
  return LAG.k::parent_type(C)
end

function base_ring_type(::Type{<:LinearAlgebraicGroup{C}}) where {C<:FieldElem}
  return parent_type(C)
end

function degree(LAG::LinearAlgebraicGroup)
  return degree(_underlying_matrix_group(LAG))
end

#This is here because of https://github.com/oscar-system/Oscar.jl/issues/5661
function gap_likes_the_group(LAG::LinearAlgebraicGroup)
  return is_finite(base_ring(LAG))
end

function has_gens(LAG::LinearAlgebraicGroup)
  gap_likes_the_group(LAG) && return has_gens(_underlying_matrix_group(LAG))
  return false
end

function number_of_generators(LAG::LinearAlgebraicGroup)
  gap_likes_the_group(LAG) && return number_of_generators(_underlying_matrix_group(LAG))
  error("Group is not finitely generated") # as long as field is QQ
end

function gens(LAG::LinearAlgebraicGroup)
  return [gen(LAG, i) for i in 1:ngens(LAG)]
end

function gen(LAG::LinearAlgebraicGroup, i::Int)
  gap_likes_the_group(LAG) &&
    return linear_algebraic_group_elem(
      LAG, gen(_underlying_matrix_group(LAG), i); check=false
    )
  error("Group is not finitely generated") # as long as field is QQ
end

function isfinite(LAG::LinearAlgebraicGroup)
  gap_likes_the_group(LAG) && return isfinite(_underlying_matrix_group(LAG))
  if degree(LAG) == 0 || degree(LAG) == 1 #Should not occur
    return true
  else
    return false
  end
end

function order(::Type{T}, LAG::LinearAlgebraicGroup) where {T}
  is_finite(LAG) || throw(InfiniteOrderError(LAG))
  return order(T, _underlying_matrix_group(LAG))
end

function Base.rand(
  rng::Random.AbstractRNG, rs::Random.SamplerTrivial{<:LinearAlgebraicGroup}
)
  LAG = rs[]
  if gap_likes_the_group(LAG)
    return linear_algebraic_group_elem(
      LAG, rand(_underlying_matrix_group(LAG)); check=false
    )
  else
    error("Random elements not implemented for this field yet.")
  end
end

function Base.eltype(::Type{LinearAlgebraicGroup{C}}) where {C<:FieldElem}
  return LinearAlgebraicGroupElem{C}
end

function elem_type(::Type{LinearAlgebraicGroup{C}}) where {C<:FieldElem}
  return LinearAlgebraicGroupElem{C}
end

function parent_type(::Type{LinearAlgebraicGroupElem{C}}) where {C<:FieldElem}
  return LinearAlgebraicGroup{C}
end

function one(LAG::LinearAlgebraicGroup)
  return linear_algebraic_group_elem(LAG, one(_underlying_matrix_group(LAG)); check=false)
end

@doc raw"""
    is_subgroup(U::MatGroup, LAG::LinearAlgebraicGroup) -> Bool

Return whether the matrix group `U` is a subgroup of `LAG` as a matrix group.
"""
function is_subgroup(U::MatGroup, LAG::LinearAlgebraicGroup)
  return is_subgroup(U, _underlying_matrix_group(LAG))
end

function Base.show(io::IO, ::MIME"text/plain", LAG::LinearAlgebraicGroup)
  io = pretty(io)
  println(
    io,
    "Linear algebraic group of type ",
    Oscar._root_system_type_string(root_system_type(root_system(LAG))),
  )
  print(io, Indent(), "over ", Lowercase(), base_ring(LAG))
  print(io, Dedent())
end

function Base.show(io::IO, LAG::LinearAlgebraicGroup)
  io = pretty(io)
  if is_terse(io)
    print(io, "LAG")
  else
    print(
      io,
      "Linear algebraic group of type ",
      Oscar._root_system_type_string(root_system_type(root_system(LAG))),
    )
    print(terse(io), " over ", Lowercase(), base_ring(LAG))
  end
end

#################################################
# Linear Algebraic Group Elements
#################################################

@doc raw"""
    linear_algebraic_group_elem(LAG::LinearAlgebraicGroup, MGE::MatGroupElem; check::Bool = true) -> LinearAlgebraicGroupElem

Coerce `MGE` into an element of `LAG`.

Setting `check` to `false` disables the check whether the element `MGE` actually lies in the group.

# Examples
```jldoctest
julia> F, _  = finite_field(5);

julia> LAG = linear_algebraic_group(:A, 2, F);

julia> m = matrix(F, [2 1 0; 1 4 3; 0 1 1]);

julia> MGE = MatGroupElem(GL(3,F), m);

julia> linear_algebraic_group_elem(LAG, MGE)
[2   1   0]
[1   4   3]
[0   1   1]

```
"""
function linear_algebraic_group_elem(
  LAG::LinearAlgebraicGroup{C}, MGE::MatGroupElem{C}; check::Bool=true
) where {C<:FieldElem}
  return LinearAlgebraicGroupElem(LAG, MGE; check)
end

@doc raw"""
    linear_algebraic_group_elem(LAG::LinearAlgebraicGroup, m::MatElem{<:FieldElem}; check::Bool = true) -> LinearAlgebraicGroupElem

Coerce `m` into an element of `LAG`.

Setting `check` to `false` disables the check whether the element `m` actually lies in the group.


# Examples
```jldoctest
julia> F, _  = finite_field(5);

julia> LAG = linear_algebraic_group(:A, 2, F);

julia> m = matrix(F, [2 1 0; 1 4 3; 0 1 1]);

julia> linear_algebraic_group_elem(LAG, m)
[2   1   0]
[1   4   3]
[0   1   1]
```
"""
function linear_algebraic_group_elem(
  LAG::LinearAlgebraicGroup{C}, m::MatElem{C}; check::Bool=true
) where {C<:FieldElem}
  MGE = _underlying_matrix_group(LAG)(m; check)
  return LinearAlgebraicGroupElem(LAG, MGE; check=false)
end

function parent(a::LinearAlgebraicGroupElem)
  return a.parent
end

function Base.:(==)(a::LinearAlgebraicGroupElem, b::LinearAlgebraicGroupElem)
  check_parent(a, b)
  return a.mat == b.mat
end

function Base.:(*)(a::LinearAlgebraicGroupElem, b::LinearAlgebraicGroupElem)
  check_parent(a, b)
  return linear_algebraic_group_elem(parent(a), a.mat * b.mat; check=false)
end

function Base.inv(a::LinearAlgebraicGroupElem)
  return linear_algebraic_group_elem(parent(a), inv(a.mat); check=false)
end

function Base.deepcopy_internal(a::LinearAlgebraicGroupElem, dict::IdDict)
  return get!(dict, a) do
    linear_algebraic_group_elem(parent(a), Base.deepcopy_internal(a.mat, dict); check=false)
  end
end

function order(::Type{T}, a::LinearAlgebraicGroupElem) where {T}
  return order(T, a.mat)
end

function Base.hash(a::LinearAlgebraicGroupElem, h::UInt)
  b = 0x1df4d55a7b37db2f % UInt
  h = hash(parent(a), h)
  h = hash(a.mat, h)

  return xor(h, b)
end

Base.show(io::IO, g::LinearAlgebraicGroupElem) = show(io, g.mat)
Base.show(io::IO, mi::MIME"text/plain", g::LinearAlgebraicGroupElem) = show(io, mi, g.mat)

############# Root Subgroups ############################
#internal function to compute action of root alpha, in case :A return the tuple (i, j) for which alpha acts like e_i-e_j
function _compute_action(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
  if root_system_type(root_system(LAG))[1][1] == :A
    c = coefficients(alpha)
    l = number_of_simple_roots(root_system(LAG))
    e = zeros(Int64, l + 1)
    for i in 1:l
      e[i] = e[i] + Int64(c[i])
      e[i + 1] = e[i + 1] - Int64(c[i])
    end
    i = 0
    j = 0
    for k in 1:(l + 1)
      if e[k] == 1
        i = k
      elseif e[k] == -1
        j = k
      end
    end
    return i, j
  else
    error("Only type A is implemented so far.")
  end
end

@doc raw"""
    root_subgroup(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem) -> MatGroup

Return the root subgroup of `LAG` corresponding to the root `alpha`.


# Examples
```jldoctest
julia> F, _  = finite_field(5);

julia> LAG = linear_algebraic_group(:A, 3, F);

julia> alpha = simple_root(root_system(LAG),2);

julia> root_subgroup(LAG, alpha)
Matrix group of degree 4
  over prime field of characteristic 5
```
"""
function root_subgroup(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
  if isdefined(LAG, :U_alphas)
    if haskey(LAG.U_alphas, alpha)
      return LAG.U_alphas[alpha]
    end
  else
    LAG.U_alphas = Dict{WeightLatticeElem,MatGroup}()
  end
  G = _underlying_matrix_group(LAG)
  i, j = _compute_action(LAG, alpha)
  gs = [
    begin
      g = identity_matrix(base_ring(LAG), degree(LAG))
      g[i, j] = lambda
      G(g)
    end for lambda in basis(base_ring(LAG))
  ]
  U, _ = sub(G, gs)
  LAG.U_alphas[alpha] = U
  return U
end

########### Tori ############################
@doc raw"""
    maximal_torus(LAG::LinearAlgebraicGroup) -> MatGroup

Return the standard maximal torus of diagonal elements of `LAG`.


# Examples
```jldoctest
julia> F, _  = finite_field(5);

julia> LAG = linear_algebraic_group(:A, 3, F);

julia> maximal_torus(LAG)
Matrix group of degree 4
  over prime field of characteristic 5
```
"""
function maximal_torus(LAG::LinearAlgebraicGroup{C}) where {C<:FieldElem}
  if !isdefined(LAG, :T)
    G = _underlying_matrix_group(LAG)
    gs = elem_type(G)[]
    t = Hecke.primitive_element(base_ring(LAG))
    it = inv(t)
    for i in 1:(degree(LAG) - 1)
      m = identity_matrix(base_ring(LAG), degree(LAG))
      m[i, i] = t
      m[i + 1, i + 1] = it
      push!(gs, G(m))
    end
    T, _ = sub(G, gs)
    LAG.T = T
  end
  return LAG.T::matrix_group_type(C)
end

@doc raw"""
    torus_element(LAG::LinearAlgebraicGroup, diag::Vector{T}) where {T<:FieldElem} -> LinearAlgebraicGroupElem

Return the root element consisting of `diag` as diagonal entries which is an element of the maximal torus of `LAG`.


# Examples
```jldoctest
julia> F, _  = finite_field(5);

julia> LAG = linear_algebraic_group(:A, 3, F);

julia> torus_element(LAG, [1,2,1,3])
[1   0   0   0]
[0   2   0   0]
[0   0   1   0]
[0   0   0   3]
```
"""
function torus_element(LAG::LinearAlgebraicGroup, diag::Vector)
  @req length(diag) == degree(LAG) "Wrong number of diagonal entries"
  diag_ = [base_ring(LAG)(d) for d in diag]
  @req isone(prod(diag_)) "Deteminant of torus element must be 1"
  m = diagonal_matrix(base_ring(LAG), diag_)
  return linear_algebraic_group_elem(LAG, m; check=false)
end

@doc raw"""
    apply_root_to_torus_element(alpha::RootSpaceElem, t::LinearAlgebraicGroupElem) -> FieldElem

Return the field element obtained by applying the root `alpha` to the torus element `t`.

# Examples
```jldoctest
julia> F, _  = finite_field(5);

julia> LAG = linear_algebraic_group(:A, 3, F);

julia> alpha = simple_root(root_system(LAG),2);

julia> t = torus_element(LAG, [F(1),F(2),F(1),F(3)]);

julia> apply_root_to_torus_element(alpha, t)
2
```
"""
function apply_root_to_torus_element(
  alpha::RootSpaceElem, t::LinearAlgebraicGroupElem
)
  @req is_root(alpha) "The given element is not a root"
  @req in(t.mat, maximal_torus(parent(t))) "The given element is not a torus element"
  i, j = _compute_action(parent(t), alpha)
  return t.mat[i, i] * inv(t.mat[j, j])
end

############### Bruhat decomposition ###########################
@doc raw"""
    representative_of_root_in_group(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem) -> LinearAlgebraicGroupElem

Return the linear algebraic group element corresponding to the root `alpha`.

# Examples
```jldoctest
julia> F, _  = finite_field(4);

julia> LAG = linear_algebraic_group(:A, 4, F);

julia> alpha = simple_root(root_system(LAG),2);

julia> representative_of_root_in_group(LAG, alpha)
[1   0   0   0   0]
[0   0   1   0   0]
[0   1   0   0   0]
[0   0   0   1   0]
[0   0   0   0   1]
```
"""
function representative_of_root_in_group(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
  @req is_root(alpha) "The given element is not a root"
  i, j = _compute_action(LAG, alpha)
  m = identity_matrix(base_ring(LAG), degree(LAG))
  m[i, i] = zero(base_ring(LAG))
  m[i, j] = -one(base_ring(LAG))
  m[j, i] = one(base_ring(LAG))
  m[j, j] = zero(base_ring(LAG))
  return linear_algebraic_group_elem(LAG, m; check=true) # TODO: eventually set check to false
end

@doc raw"""
    borel_subgroup(LAG::LinearAlgebraicGroup) -> MatGroup

Return the standard Borel subgroup of the linear algebraic group `LAG`.

# Examples
```jldoctest
julia> F, _  = finite_field(3);

julia> LAG = linear_algebraic_group(:A, 3, F);

julia> borel_subgroup(LAG)
Matrix group of degree 4
  over prime field of characteristic 3
```
"""
function borel_subgroup(LAG::LinearAlgebraicGroup{C}) where {C<:FieldElem}
  if !isdefined(LAG, :B)
    T = maximal_torus(LAG)
    G = _underlying_matrix_group(LAG)
    gs = elem_type(G)[]
    for t in gens(T)
      push!(gs, t)
    end
    for alpha in simple_roots(root_system(LAG))
      i, j = _compute_action(LAG, alpha)
      for lambda in basis(base_ring(LAG))
        g = identity_matrix(base_ring(LAG), degree(LAG))
        g[i, j] = lambda
        push!(gs, G(g))
      end
    end
    B, _ = sub(G, gs)
    LAG.B = B
  end
  return LAG.B::matrix_group_type(C)
end

@doc raw"""
    bruhat_cell_representative(LAG::LinearAlgebraicGroup, w::WeylGroupElem) -> MatGroupElem

Return the representative of the Bruhat cell corresponding to the Weyl group element `w`.

# Examples
```jldoctest
julia> F, _  = finite_field(5);

julia> LAG = linear_algebraic_group(:A, 3, F);

julia> W = weyl_group(root_system(LAG));

julia> w = W([1,2]);

julia> bruhat_cell_representative(LAG, w)
[0   0   1   0]
[1   0   0   0]
[0   1   0   0]
[0   0   0   1]
```
"""
function bruhat_cell_representative(LAG::LinearAlgebraicGroup, w::WeylGroupElem)
  @req parent(w) == weyl_group(root_system(LAG)) "parent mismatch"
  rep = identity_matrix(base_ring(LAG), degree(LAG))
  for i in word(w)
    alpha = simple_root(root_system(LAG), Int64(i))
    rep = rep * representative_of_root_in_group(LAG, alpha).mat
  end
  return _underlying_matrix_group(LAG)(rep)
end

@doc raw"""
    bruhat_cell(LAG::LinearAlgebraicGroup, w::WeylGroupElem) -> GroupDoubleCoset{MatGroup, MatGroupElem}

Return the Bruhat cell corresponding to the Weyl group element `w`.

# Examples
```jldoctest
julia> F, _  = finite_field(5);

julia> LAG = linear_algebraic_group(:A, 3, F);

julia> W = weyl_group(root_system(LAG));

julia> w = W([1,2]);

julia> bruhat_cell(LAG,w)
Double coset of matrix group of degree 4 over F
  and matrix group of degree 4 over F
  with representative [0 0 1 0; 1 0 0 0; 0 1 0 0; 0 0 0 1]
  in SL(4,5)
```
"""
function bruhat_cell(LAG::LinearAlgebraicGroup, w::WeylGroupElem)
  rep = bruhat_cell_representative(LAG, w)
  B = borel_subgroup(LAG)
  return double_coset(B, rep, B)
end

@doc raw"""
    bruhat_decomposition(LAG::LinearAlgebraicGroup) -> Vector{GroupDoubleCoset{MatGroup, MatGroupElem}}

Return the Bruhat decomposition of the linear algebraic group `LAG`.

# Examples
```jldoctest
julia> F, _  = finite_field(5);

julia> LAG = linear_algebraic_group(:A, 2, F);

julia> bruhat_decomposition(LAG)
6-element Vector{GroupDoubleCoset{MatGroup{FqFieldElem, FqMatrix}, MatGroupElem{FqFieldElem, FqMatrix}}}:
 Double coset of matrix group and matrix group with representative [1 0 0; 0 1 0; 0 0 1]
 Double coset of matrix group and matrix group with representative [0 4 0; 1 0 0; 0 0 1]
 Double coset of matrix group and matrix group with representative [0 4 0; 0 0 4; 1 0 0]
 Double coset of matrix group and matrix group with representative [0 0 1; 0 4 0; 1 0 0]
 Double coset of matrix group and matrix group with representative [1 0 0; 0 0 4; 0 1 0]
 Double coset of matrix group and matrix group with representative [0 0 1; 1 0 0; 0 1 0]
```
"""
function bruhat_decomposition(LAG::LinearAlgebraicGroup)
  return [bruhat_cell(LAG, w) for w in weyl_group(root_system(LAG))]
end
