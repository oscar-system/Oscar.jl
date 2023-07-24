
###############################################################################
#
#  Accessors
#
###############################################################################

@doc raw"""
    lattice(Lf::ZZLatWithIsom) -> ZZLat

Given a lattice with isometry $(L, f)$, return the underlying lattice `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L; neg=true);

julia> lattice(Lf) === L
true
```
"""
lattice(Lf::ZZLatWithIsom) = Lf.Lb

@doc raw"""
    isometry(Lf::ZZLatWithIsom) -> QQMatrix

Given a lattice with isometry $(L, f)$, return the underlying isometry `f`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L; neg=true);

julia> isometry(Lf)
[-1    0    0    0    0]
[ 0   -1    0    0    0]
[ 0    0   -1    0    0]
[ 0    0    0   -1    0]
[ 0    0    0    0   -1]
```
"""
isometry(Lf::ZZLatWithIsom) = Lf.f

@doc raw"""
    ambient_space(Lf::ZZLatWithIsom) -> QuadSpaceWithIsom

Given a lattice with isometry $(L, f)$, return the pair $(V, g)$ where
`V` is the ambient quadratic space of `L` and `g` is an isometry of `V`
inducing `f` on `L`.

Note that `g` might not be unique and we fix such a global isometry
together with `V` into a container type `QuadSpaceWithIsom`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L; neg=true);

julia> Vf = ambient_space(Lf)
Quadratic space of dimension 5
  with isometry of finite order 2
  given by
  [-1    0    0    0    0]
  [ 0   -1    0    0    0]
  [ 0    0   -1    0    0]
  [ 0    0    0   -1    0]
  [ 0    0    0    0   -1]

julia> typeof(Vf)
QuadSpaceWithIsom
```
"""
ambient_space(Lf::ZZLatWithIsom) = Lf.Vf

@doc raw"""
   ambient_isometry(Lf::ZZLatWithIsom) -> QQMatrix

Given a lattice with isometry $(L, f)$, return an isometry of the ambient
space of `L` inducing `f` on `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L; neg=true);

julia> ambient_isometry(Lf)
[-1    0    0    0    0]
[ 0   -1    0    0    0]
[ 0    0   -1    0    0]
[ 0    0    0   -1    0]
[ 0    0    0    0   -1]
```
"""
ambient_isometry(Lf::ZZLatWithIsom) = isometry(ambient_space(Lf))

@doc raw"""
    order_of_isometry(Lf::ZZLatWithIsom) -> IntExt

Given a lattice with isometry $(L, f)$, return the order of the underlying
isometry `f`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L; neg=true);

julia> order_of_isometry(Lf) == 2
true
```
"""
order_of_isometry(Lf::ZZLatWithIsom) = Lf.n

###############################################################################
#
#  Attributes
#
###############################################################################

@doc raw"""
    rank(Lf::ZZLatWithIsom) -> Integer

Given a lattice with isometry $(L, f)$, return the rank of the underlying lattice
`L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L; neg=true);

julia> rank(Lf)
5
```
"""
rank(Lf::ZZLatWithIsom) = rank(lattice(Lf))

@doc raw"""
    characteristic_polynomial(Lf::ZZLatWithIsom) -> QQPolyRingElem

Given a lattice with isometry $(L, f)$, return the characteristic polynomial of the
underlying isometry `f`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L; neg=true);

julia> factor(characteristic_polynomial(Lf))
1 * (x + 1)^5
```
"""
characteristic_polynomial(Lf::ZZLatWithIsom) = characteristic_polynomial(isometry(Lf))

@doc raw"""
    minimal_polynomial(Lf::ZZLatWithIsom) -> QQPolyRingElem

Given a lattice with isometry $(L, f)$, return the minimal polynomial of the
underlying isometry `f`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L; neg=true);

julia> minimal_polynomial(Lf)
x + 1
```
"""
minimal_polynomial(Lf::ZZLatWithIsom) = minimal_polynomial(isometry(Lf))

@doc raw"""
    genus(Lf::ZZLatWithIsom) -> ZZGenus

Given a lattice with isometry $(L, f)$, return the genus of the underlying 
lattice `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L; neg=true);

julia> genus(Lf)
Genus symbol for integer lattices
Signatures: (5, 0, 0)
Local symbols:
  Local genus symbol at 2: 1^-4 2^1_7
  Local genus symbol at 3: 1^-4 3^1
```
"""
genus(Lf::ZZLatWithIsom) = genus(lattice(Lf))

@doc raw"""
    basis_matrix(Lf::ZZLatWithIsom) -> QQMatrix

Given a lattice with isometry $(L, f)$, return the basis matrix of the underlying
lattice `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ,5,5,[ 1  0  0  0  0;
                          -1 -1 -1 -1 -1;
                           0  0  0  0  1;
                           0  0  0  1  0;
                           0  0  1  0  0])
[ 1    0    0    0    0]
[-1   -1   -1   -1   -1]
[ 0    0    0    0    1]
[ 0    0    0    1    0]
[ 0    0    1    0    0]

julia> Lf = integer_lattice_with_isometry(L, f);

julia> invariant_lattice(Lf);

julia> I = invariant_lattice(Lf);

julia> basis_matrix(I)
[1   0   0   0   0]
[0   0   1   0   1]
[0   0   0   1   0]
```
"""
basis_matrix(Lf::ZZLatWithIsom) = basis_matrix(lattice(Lf))

@doc raw"""
    gram_matrix(Lf::ZZLatWithIsom) -> QQMatrix

Given a lattice with isometry $(L, f)$, return the gram matrix of the underlying
lattice `L` with respect to its basis matrix.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> gram_matrix(Lf)
[ 2   -1    0    0    0]
[-1    2   -1    0    0]
[ 0   -1    2   -1    0]
[ 0    0   -1    2   -1]
[ 0    0    0   -1    2]
```
"""
gram_matrix(Lf::ZZLatWithIsom) = gram_matrix(lattice(Lf))

@doc raw"""
    rational_span(Lf::ZZLatWithIsom) -> QuadSpaceWithIsom

Given a lattice with isometry $(L, f)$, return the rational span
$L \otimes \mathbb{Q}$ of the underlying lattice `L` together with the
underlying isometry of `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> Vf = rational_span(Lf)
Quadratic space of dimension 5
  with isometry of finite order 1
  given by
  [1   0   0   0   0]
  [0   1   0   0   0]
  [0   0   1   0   0]
  [0   0   0   1   0]
  [0   0   0   0   1]

julia> typeof(Vf)
QuadSpaceWithIsom
```
"""
rational_span(Lf::ZZLatWithIsom) = quadratic_space_with_isometry(rational_span(lattice(Lf)), isometry(Lf))

@doc raw"""
    det(Lf::ZZLatWithIsom) -> QQFieldElem

Given a lattice with isometry $(L, f)$, return the determinant of the
underlying lattice `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> det(Lf)
6
```
"""
det(Lf::ZZLatWithIsom) = det(lattice(Lf))

@doc raw"""
    scale(Lf::ZZLatWithIsom) -> QQFieldElem

Given a lattice with isometry $(L, f)$, return the scale of the underlying
lattice `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> scale(Lf)
1
```
"""
scale(Lf::ZZLatWithIsom) = scale(lattice(Lf))

@doc raw"""
    norm(Lf::ZZLatWithIsom) -> QQFieldElem

Given a lattice with isometry $(L, f)$, return the norm of the underlying
lattice `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> norm(Lf)
2
```
"""
norm(Lf::ZZLatWithIsom) = norm(lattice(Lf))

@doc raw"""
    is_positive_definite(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the underlying
lattice `L` is positive definite.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> is_positive_definite(Lf)
true
```
"""
is_positive_definite(Lf::ZZLatWithIsom) = is_positive_definite(lattice(Lf))

@doc raw"""
    is_negative_definite(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the underlying
lattice `L` is negative definite.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> is_negative_definite(Lf)
false
```
"""
is_negative_definite(Lf::ZZLatWithIsom) = is_negative_definite(lattice(Lf))

@doc raw"""
    is_definite(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the underlying
lattice `L` is definite.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> is_definite(Lf)
true
```
"""
is_definite(Lf::ZZLatWithIsom) = is_definite(lattice(Lf))

@doc raw"""
    minimum(Lf::ZZLatWithIsom) -> QQFieldElem

Given a positive definite lattice with isometry $(L, f)$, return the minimum
of the underlying lattice `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> minimum(Lf)
2
```
"""
function minimum(Lf::ZZLatWithIsom)
  @req is_positive_definite(Lf) "Underlying lattice must be positive definite"
  return minimum(lattice(Lf))
end

@doc raw"""
    is_integral(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the underlying lattice
`L` is integral.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> is_integral(Lf)
true
```
"""
is_integral(Lf::ZZLatWithIsom) = is_integral(lattice(Lf))

@doc raw"""
    degree(Lf::ZZLatWithIsom) -> Int

Given a lattice with isometry $(L, f)$, return the degree of the underlying
lattice `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> degree(Lf)
5
```
"""
degree(Lf::ZZLatWithIsom) = degree(lattice(Lf))

@doc raw"""
    is_even(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the underlying lattice
`L` is even.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> is_even(Lf)
true
```
"""
is_even(Lf::ZZLatWithIsom) = iseven(lattice(Lf))

@doc raw"""
    discriminant(Lf::ZZLatWithIsom) -> QQFieldElem

Given a lattice with isometry $(L, f)$, return the discriminant of the underlying
lattice `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> discriminant(Lf) == det(Lf) == 6
true
```
"""
discriminant(Lf::ZZLatWithIsom) = discriminant(lattice(Lf))

@doc raw"""
    signature_tuple(Lf::ZZLatWithIsom) -> Tuple{Int, Int, Int}

Given a lattice with isometry $(L, f)$, return the signature tuple of the
underlying lattice `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L);

julia> signature_tuple(Lf)
(5, 0, 0)
```
"""
signature_tuple(Lf::ZZLatWithIsom) = signature_tuple(lattice(Lf))

###############################################################################
#
#  Constructors
#
###############################################################################

@doc raw"""
    integer_lattice_with_isometry(L::ZZLat, f::QQMatrix; check::Bool = true,
                                                ambient_representation = true)
				                                                             -> ZZLatWithIsom

Given a $\mathbb Z$-lattice `L` and a matrix `f`, if `f` defines an isometry
of `L` of order `n`, return the corresponding lattice with isometry pair $(L, f)$.

If `ambient_representation` is set to true, `f` is consider as an isometry of the
ambient space of `L` and the induced isometry on `L` is automatically computed.
Otherwise, an isometry of the ambient space of `L` is constructed, setting the identity
on the complement of the rational span of `L` if it is not of full rank.

# Examples

The way we construct the lattice can have an influence on the isometry of the
ambient space we store. Indeed, if one mentions an isometry of the lattice,
this isometry will be extended by the identity on the orthogonal complement
of the rational span of the lattice. In the following example, `Lf` and `Lf2` are
the same object, but the isometry of their ambient space stored are
different (one has order 2, the other one is the identity).

```jldoctest
julia> B = matrix(QQ, 3, 5, [1 0 0 0 0;
                             0 0 1 0 1;
                             0 0 0 1 0]);

julia> G = matrix(QQ, 5, 5, [ 2 -1  0  0  0;
                             -1  2 -1  0  0;
                              0 -1  2 -1  0;
                              0  0 -1  2 -1;
                              0  0  0 -1  2]);

julia> L = integer_lattice(B; gram = G);

julia> f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0]);

julia> Lf = integer_lattice_with_isometry(L, f)
Integer lattice of rank 3 and degree 5
  with isometry of finite order 1
  given by
  [1   0   0]
  [0   1   0]
  [0   0   1]

julia> ambient_isometry(Lf)
[ 1    0    0    0    0]
[-1   -1   -1   -1   -1]
[ 0    0    0    0    1]
[ 0    0    0    1    0]
[ 0    0    1    0    0]

julia> Lf2 = integer_lattice_with_isometry(L, isometry(Lf); ambient_representation=false)
Integer lattice of rank 3 and degree 5
  with isometry of finite order 1
  given by
  [1   0   0]
  [0   1   0]
  [0   0   1]

julia> ambient_isometry(Lf2)
[1   0   0   0   0]
[0   1   0   0   0]
[0   0   1   0   0]
[0   0   0   1   0]
[0   0   0   0   1]
```
"""
function integer_lattice_with_isometry(L::ZZLat, f::QQMatrix; check::Bool = true,
                                                              ambient_representation::Bool = true)
  if rank(L) == 0
    Vf = quadratic_space_with_isometry(ambient_space(L))
    return ZZLatWithIsom(Vf, L, matrix(QQ,0,0,[]), -1)
  end

  if check
    @req det(f) != 0 "f is not invertible"
  end

  if ambient_representation
    f_ambient = f
    Vf = quadratic_space_with_isometry(ambient_space(L), f_ambient; check=check)
    B = basis_matrix(L)
    ok, f = can_solve_with_solution(B, B*f_ambient; side = :left)
    @req ok "Isometry does not restrict to L"
  else
    V = ambient_space(L)
    B = basis_matrix(L)
    B2 = orthogonal_complement(V, B)
    C = vcat(B, B2)
    f_ambient = block_diagonal_matrix([f, identity_matrix(QQ, nrows(B2))])
    f_ambient = inv(C)*f_ambient*C
    Vf = quadratic_space_with_isometry(V, f_ambient; check=check)
  end

  n = multiplicative_order(f)

  if check
    @req f*gram_matrix(L)*transpose(f) == gram_matrix(L) "f does not define an isometry of L"
    @hassert :ZZLatWithIsom 1 basis_matrix(L)*f_ambient == f*basis_matrix(L)
  end

  return ZZLatWithIsom(Vf, L, f, n)::ZZLatWithIsom
end

@doc raw"""
    integer_lattice_with_isometry(L::ZZLat; neg::Bool = false) -> ZZLatWithIsom

Given a $\mathbb Z$-lattice `L` return the lattice with isometry pair $(L, f)$,
where `f` corresponds to the identity mapping of `L`.

If `neg` is set to `true`, then the isometry `f` is negative the identity of `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> Lf = integer_lattice_with_isometry(L; neg=true)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 2
  given by
  [-1    0    0    0    0]
  [ 0   -1    0    0    0]
  [ 0    0   -1    0    0]
  [ 0    0    0   -1    0]
  [ 0    0    0    0   -1]
```
"""
function integer_lattice_with_isometry(L::ZZLat; neg::Bool = false)
  d = degree(L)
  f = identity_matrix(QQ, d)
  if neg
    f = -f
  end
  return integer_lattice_with_isometry(L, f; check = false, ambient_representation = true)::ZZLatWithIsom
end

@doc raw"""
    lattice(Vf::QuadSpaceWithIsom) -> ZZLatWithIsom

Given a quadratic space with isometry $(V, f)$, return the full rank lattice `L`
in `V` with basis the standard basis, together with the induced action of `f`
on `L`.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, QQ[2 -1; -1 2])
Quadratic space of dimension 2
  over rational field
with gram matrix
[ 2   -1]
[-1    2]

julia> f = matrix(QQ, 2, 2, [1 1; 0 -1])
[1    1]
[0   -1]

julia> Vf = quadratic_space_with_isometry(V, f)
Quadratic space of dimension 2
  with isometry of finite order 2
  given by
  [1    1]
  [0   -1]

julia> Lf = lattice(Vf)
Integer lattice of rank 2 and degree 2
  with isometry of finite order 2
  given by
  [1    1]
  [0   -1]
```
"""
lattice(Vf::QuadSpaceWithIsom) = ZZLatWithIsom(Vf, lattice(space(Vf)), isometry(Vf), order_of_isometry(Vf))::ZZLatWithIsom

@doc raw"""
    lattice(Vf::QuadSpaceWithIsom, B::MatElem{<:RationalUnion};
                                   isbasis::Bool = true, check::Bool = true)
                                                              -> ZZLatWithIsom

Given a quadratic space with isometry $(V, f)$ and a matrix `B` generating a
lattice `L` in `V`, if `L` is preserved under the action of `f`, return the
lattice with isometry $(L, f_L)$ where $f_L$ is induced by the action of `f`
on `L`.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, QQ[ 2 -1  0  0  0;
                                  -1  2 -1  0  0;
                                   0 -1  2 -1  0;
                                   0  0 -1  2 -1;
                                   0  0  0 -1  2]);

julia> f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0]);

julia> Vf = quadratic_space_with_isometry(V, f);

julia> B = matrix(QQ,3,5,[1 0 0 0 0;
                          0 0 1 0 1;
                          0 0 0 1 0])
[1   0   0   0   0]
[0   0   1   0   1]
[0   0   0   1   0]

julia> lattice(Vf, B)
Integer lattice of rank 3 and degree 5
  with isometry of finite order 1
  given by
  [1   0   0]
  [0   1   0]
  [0   0   1]
```
"""
function lattice(Vf::QuadSpaceWithIsom, B::MatElem{<:RationalUnion}; isbasis::Bool = true, check::Bool = true)
  L = lattice(space(Vf), B; isbasis=isbasis, check=check)
  ok, fB = can_solve_with_solution(basis_matrix(L), basis_matrix(L)*isometry(Vf); side = :left)
  n = is_zero(fB) ? -1 : multiplicative_order(fB)
  @req ok "The lattice defined by B is not preserved under the action of the isometry of Vf"
  return ZZLatWithIsom(Vf, L, fB, n)
end

@doc raw"""
    lattice_in_same_ambient_space(L::ZZLatWithIsom, B::MatElem;
                                                    check::Bool = true)
                                                            -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$ and a matrix `B` whose rows define a free
system of vectors in the ambient space `V` of `L`, if the lattice `M` in `V` defined
by `B` is preserved under the fixed isometry `g` of `V` inducing `f` on `L`, return
the lattice with isometry pair $(M, f_M)$ where $f_M$ is induced by the action of
`g` on `M`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> B = matrix(QQ,3,5,[1 0 0 0 0;
                          0 0 1 0 1;
                          0 0 0 1 0])
[1   0   0   0   0]
[0   0   1   0   1]
[0   0   0   1   0]

julia> I = lattice_in_same_ambient_space(Lf, B)
Integer lattice of rank 3 and degree 5
  with isometry of finite order 1
  given by
  [1   0   0]
  [0   1   0]
  [0   0   1]

julia> ambient_space(I) === ambient_space(Lf)
true
```
"""
function lattice_in_same_ambient_space(L::ZZLatWithIsom, B::MatElem; check::Bool = true)
  @req !check || (rank(B) == nrows(B)) "The rows of B must define a free system of vectors"
  Vf = ambient_space(L)
  return lattice(Vf, B; check = check)
end

###############################################################################
#
#  Operations on lattices with isometry
#
###############################################################################

@doc raw"""
    rescale(Lf::ZZLatWithIsom, a::RationalUnion) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$ and a rational number `a`, return the
lattice with isometry $(L(a), f)$.

# Examples
```jldoctest
julia> L = root_lattice(:A,5)
Integer lattice of rank 5 and degree 5
with gram matrix
[ 2   -1    0    0    0]
[-1    2   -1    0    0]
[ 0   -1    2   -1    0]
[ 0    0   -1    2   -1]
[ 0    0    0   -1    2]

julia> f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0]);

julia> Lf = integer_lattice_with_isometry(L, f)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 2
  given by
  [ 1    0    0    0    0]
  [-1   -1   -1   -1   -1]
  [ 0    0    0    0    1]
  [ 0    0    0    1    0]
  [ 0    0    1    0    0]

julia> Lf2 = rescale(Lf, 1//2)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 2
  given by
  [ 1    0    0    0    0]
  [-1   -1   -1   -1   -1]
  [ 0    0    0    0    1]
  [ 0    0    0    1    0]
  [ 0    0    1    0    0]

julia> lattice(Lf2)
Integer lattice of rank 5 and degree 5
with gram matrix
[    1   -1//2       0       0       0]
[-1//2       1   -1//2       0       0]
[    0   -1//2       1   -1//2       0]
[    0       0   -1//2       1   -1//2]
[    0       0       0   -1//2       1]
```
"""
function rescale(Lf::ZZLatWithIsom, a::RationalUnion)
  return lattice(rescale(ambient_space(Lf), a), basis_matrix(Lf); check=false)
end

@doc raw"""
    ^(Lf::ZZLatWithIsom, n::Int) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$ and an integer $n$, return the pair
$(L, f^n)$.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0]);

julia> Lf = integer_lattice_with_isometry(L, f)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 2
  given by
  [ 1    0    0    0    0]
  [-1   -1   -1   -1   -1]
  [ 0    0    0    0    1]
  [ 0    0    0    1    0]
  [ 0    0    1    0    0]

julia> Lf^0
Integer lattice of rank 5 and degree 5
  with isometry of finite order 1
  given by
  [1   0   0   0   0]
  [0   1   0   0   0]
  [0   0   1   0   0]
  [0   0   0   1   0]
  [0   0   0   0   1]
```
"""
function Base.:^(Lf::ZZLatWithIsom, n::Int)
  return lattice(ambient_space(Lf)^n, basis_matrix(Lf))
end

@doc raw"""
    dual(Lf::ZZLatWithIsom) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$ inside the space $(V, \Phi)$, such that
`f` is induced by an isometry `g` of $(V, \Phi)$, return the lattice with
isometry $(L^{\vee}, h)$ where $L^{\vee}$ is the dual of `L` in $(V, \Phi)$
and `h` is induced by `g`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0]);

julia> Lf = integer_lattice_with_isometry(L, f)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 2
  given by
  [ 1    0    0    0    0]
  [-1   -1   -1   -1   -1]
  [ 0    0    0    0    1]
  [ 0    0    0    1    0]
  [ 0    0    1    0    0]

julia> Lfv = dual(Lf)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 2
  given by
  [1   -1   0   0   0]
  [0   -1   0   0   0]
  [0   -1   0   0   1]
  [0   -1   0   1   0]
  [0   -1   1   0   0]

julia> ambient_space(Lfv) == ambient_space(Lf)
true
```
"""
function dual(Lf::ZZLatWithIsom)
  @req is_integral(Lf) "Underlying lattice must be integral"
  return lattice_in_same_ambient_space(Lf, basis_matrix(dual(lattice(Lf))); check = false)
end

@doc raw"""
    lll(Lf::ZZLatWithIsom) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$, return the same lattice with isometry
with a different basis matrix for `L`  obtained by performing an LLL-reduction
on the associated gram matrix of `L`.

Note that matrix representing the action of `f` on `L` changes but the global action
on the ambient space of `L` stays the same.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0]);

julia> Lf = integer_lattice_with_isometry(L, f)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 2
  given by
  [ 1    0    0    0    0]
  [-1   -1   -1   -1   -1]
  [ 0    0    0    0    1]
  [ 0    0    0    1    0]
  [ 0    0    1    0    0]

julia> Lf2 = lll(Lf)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 2
  given by
  [ 1    0    0    0    0]
  [-1    0    0    0   -1]
  [-1    0    0   -1    0]
  [-1    0   -1    0    0]
  [-1   -1    0    0    0]

julia> ambient_space(Lf2) == ambient_space(Lf)
true
```
"""
function lll(Lf::ZZLatWithIsom; same_ambient::Bool = true)
  L2 = lll(lattice(Lf); same_ambient = same_ambient)
  if same_ambient
    return lattice_in_same_ambient_space(Lf, basis_matrix(L2); check=false)
  else
    return integer_lattice_with_isometry(L2, isometry(Lf); check=false, ambient_representation=false)
  end
end

@doc raw"""
    direct_sum(x::Vector{ZZLatWithIsom}) -> ZZLatWithIsom,
                                                       Vector{AbstractSpaceMor}
    direct_sum(x::Vararg{ZZLatWithIsom}) -> ZZLatWithIsom,
                                                       Vector{AbstractSpaceMor}

Given a collection of lattices with isometries $(L_1, f_1) \ldots, (L_n, f_n)$,
return the lattice with isometry $(L, f)$ together with the injections $L_i \to L$,
where `L` is the direct sum $L := L_1 \oplus \ldots \oplus L_n$ and `f` is the
isometry of `L` induced by the diagonal actions of the $f_i$'s.

For objects of type `ZZLatWithIsom`, finite direct sums and finite direct products
agree and they are therefore called biproducts.
If one wants to obtain $(L, f)$ as a direct product with the projections $L \to L_i$,
one should call `direct_product(x)`.
If one wants to obtain $(L, f)$ as a biproduct with the injections $L_i \to L$ and
the projections $L \to L_i$, one should call `biproduct(x)`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0]);

julia> g = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 2
  given by
  [ 1    0    0    0    0]
  [-1   -1   -1   -1   -1]
  [ 0    0    0    0    1]
  [ 0    0    0    1    0]
  [ 0    0    1    0    0]

julia> Lg = integer_lattice_with_isometry(L, g)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 5
  given by
  [1    1    1    1    1]
  [0   -1   -1   -1   -1]
  [0    1    0    0    0]
  [0    0    1    0    0]
  [0    0    0    1    0]

julia> Lh, inj = direct_sum(Lf, Lg)
(Integer lattice with isometry of finite order 10, AbstractSpaceMor[Map with following data
Domain:
=======
Quadratic space of dimension 5
Codomain:
=========
Quadratic space of dimension 10, Map with following data
Domain:
=======
Quadratic space of dimension 5
Codomain:
=========
Quadratic space of dimension 10])

julia> Lh
Integer lattice of rank 10 and degree 10
  with isometry of finite order 10
  given by
  [ 1    0    0    0    0   0    0    0    0    0]
  [-1   -1   -1   -1   -1   0    0    0    0    0]
  [ 0    0    0    0    1   0    0    0    0    0]
  [ 0    0    0    1    0   0    0    0    0    0]
  [ 0    0    1    0    0   0    0    0    0    0]
  [ 0    0    0    0    0   1    1    1    1    1]
  [ 0    0    0    0    0   0   -1   -1   -1   -1]
  [ 0    0    0    0    0   0    1    0    0    0]
  [ 0    0    0    0    0   0    0    1    0    0]
  [ 0    0    0    0    0   0    0    0    1    0]
```
"""
function direct_sum(x::Vector{ZZLatWithIsom})
  Vf, inj = direct_sum(ambient_space.(x))
  Bs = diagonal_matrix(basis_matrix.(x))
  return lattice(Vf, Bs; check=false), inj
end

direct_sum(x::Vararg{ZZLatWithIsom}) = direct_sum(collect(x))

@doc raw"""
    direct_product(x::Vector{ZZLatWithIsom}) -> ZZLatWithIsom,
                                                       Vector{AbstractSpaceMor}
    direct_product(x::Vararg{ZZLatWithIsom}) -> ZZLatWithIsom,
                                                       Vector{AbstractSpaceMor}

Given a collection of lattices with isometries $(L_1, f_1), \ldots, (L_n, f_n)$,
return the lattice with isometry $(L, f)$ together with the projections $L \to L_i$,
where `L` is the direct product $L := L_1 \times \ldots \times L_n$ and `f` is the
isometry of `L` induced by the diagonal actions of the $f_i$'s.

For objects of type `ZZLatWithIsom`, finite direct sums and finite direct products
agree and they are therefore called biproducts.
If one wants to obtain $(L, f)$ as a direct sum with the injections $L_i \to L$,
one should call `direct_sum(x)`.
If one wants to obtain $(L, f)$ as a biproduct with the injections $L_i \to L$ and
the projections $L \to L_i$, one should call `biproduct(x)`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0]);

julia> g = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 2
  given by
  [ 1    0    0    0    0]
  [-1   -1   -1   -1   -1]
  [ 0    0    0    0    1]
  [ 0    0    0    1    0]
  [ 0    0    1    0    0]

julia> Lg = integer_lattice_with_isometry(L, g)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 5
  given by
  [1    1    1    1    1]
  [0   -1   -1   -1   -1]
  [0    1    0    0    0]
  [0    0    1    0    0]
  [0    0    0    1    0]

julia> Lh, proj = direct_product(Lf, Lg)
(Integer lattice with isometry of finite order 10, AbstractSpaceMor[Map with following data
Domain:
=======
Quadratic space of dimension 10
Codomain:
=========
Quadratic space of dimension 5, Map with following data
Domain:
=======
Quadratic space of dimension 10
Codomain:
=========
Quadratic space of dimension 5])

julia> Lh
Integer lattice of rank 10 and degree 10
  with isometry of finite order 10
  given by
  [ 1    0    0    0    0   0    0    0    0    0]
  [-1   -1   -1   -1   -1   0    0    0    0    0]
  [ 0    0    0    0    1   0    0    0    0    0]
  [ 0    0    0    1    0   0    0    0    0    0]
  [ 0    0    1    0    0   0    0    0    0    0]
  [ 0    0    0    0    0   1    1    1    1    1]
  [ 0    0    0    0    0   0   -1   -1   -1   -1]
  [ 0    0    0    0    0   0    1    0    0    0]
  [ 0    0    0    0    0   0    0    1    0    0]
  [ 0    0    0    0    0   0    0    0    1    0]
```
"""
function direct_product(x::Vector{ZZLatWithIsom})
  Vf, proj = direct_product(ambient_space.(x))
  Bs = diagonal_matrix(basis_matrix.(x))
  return lattice(Vf, Bs; check=false), proj
end

direct_product(x::Vararg{ZZLatWithIsom}) = direct_product(collect(x))

@doc raw"""
    biproduct(x::Vector{ZZLatWithIsom}) -> ZZLatWithIsom,
                                                      Vector{AbstractSpaceMor},
                                                      Vector{AbstractSpaceMor}
    biproduct(x::Vararg{ZZLatWithIsom}) -> ZZLatWithIsom,
                                                      Vector{AbstractSpaceMor},
                                                      Vector{AbstractSpaceMor}

Given a collection of lattices with isometries $(L_1, f_1), \ldots, (L_n, f_n)$,
return the lattice with isometry $(L, f)$ together with the injections
$L_i \to L$ and the projections $L \to L_i$, where `L` is the biproduct
$L := L_1 \oplus \ldots \oplus L_n$ and `f` is the isometry of `L` induced by the
diagonal actions of the $f_i$'s.

For objects of type `ZZLatWithIsom`, finite direct sums and finite direct products
agree and they are therefore called biproducts.
If one wants to obtain $(L, f)$ as a direct sum with the injections $L_i \to L$,
one should call `direct_sum(x)`.
If one wants to obtain $(L, f)$ as a direct product with the projections $L \to L_i$,
one should call `direct_product(x)`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0]);

julia> g = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 2
  given by
  [ 1    0    0    0    0]
  [-1   -1   -1   -1   -1]
  [ 0    0    0    0    1]
  [ 0    0    0    1    0]
  [ 0    0    1    0    0]

julia> Lg = integer_lattice_with_isometry(L, g)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 5
  given by
  [1    1    1    1    1]
  [0   -1   -1   -1   -1]
  [0    1    0    0    0]
  [0    0    1    0    0]
  [0    0    0    1    0]

julia> Lh, inj, proj = biproduct(Lf, Lg)
(Integer lattice with isometry of finite order 10, AbstractSpaceMor[Map with following data
Domain:
=======
Quadratic space of dimension 5
Codomain:
=========
Quadratic space of dimension 10, Map with following data
Domain:
=======
Quadratic space of dimension 5
Codomain:
=========
Quadratic space of dimension 10], AbstractSpaceMor[Map with following data
Domain:
=======
Quadratic space of dimension 10
Codomain:
=========
Quadratic space of dimension 5, Map with following data
Domain:
=======
Quadratic space of dimension 10
Codomain:
=========
Quadratic space of dimension 5])

julia> Lh
Integer lattice of rank 10 and degree 10
  with isometry of finite order 10
  given by
  [ 1    0    0    0    0   0    0    0    0    0]
  [-1   -1   -1   -1   -1   0    0    0    0    0]
  [ 0    0    0    0    1   0    0    0    0    0]
  [ 0    0    0    1    0   0    0    0    0    0]
  [ 0    0    1    0    0   0    0    0    0    0]
  [ 0    0    0    0    0   1    1    1    1    1]
  [ 0    0    0    0    0   0   -1   -1   -1   -1]
  [ 0    0    0    0    0   0    1    0    0    0]
  [ 0    0    0    0    0   0    0    1    0    0]
  [ 0    0    0    0    0   0    0    0    1    0]

julia> matrix(compose(inj[1], proj[1]))
[1   0   0   0   0]
[0   1   0   0   0]
[0   0   1   0   0]
[0   0   0   1   0]
[0   0   0   0   1]

julia> matrix(compose(inj[1], proj[2]))
[0   0   0   0   0]
[0   0   0   0   0]
[0   0   0   0   0]
[0   0   0   0   0]
[0   0   0   0   0]
```
"""
function biproduct(x::Vector{ZZLatWithIsom})
  Vf, inj, proj = biproduct(ambient_space.(x))
  Bs = diagonal_matrix(basis_matrix.(x))
  return lattice(Vf, Bs; check=false), inj, proj
end

biproduct(x::Vararg{ZZLatWithIsom}) = biproduct(collect(x))

###############################################################################
#
#  Equality and hash
#
###############################################################################

function Base.:(==)(L1::ZZLatWithIsom, L2::ZZLatWithIsom)
  ambient_space(L1) == ambient_space(L2) || return false
  return lattice(L1) == lattice(L2)
end

function Base.hash(L::ZZLatWithIsom, u::UInt)
  u = Base.hash(ambient_space(L), u)
  return Base.hash(lattice(L), u)
end

###############################################################################
#
#  Hermitian structure
#
###############################################################################

@doc raw"""
    is_of_hermitian_type(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the minimal polynomial of
the underlying isometry `f` is irreducible.

Note that if $(L, f)$ is of hermitian type with `f` of minimal polynomial $\chi$,
then `L` can be seen as a hermitian lattice over the order $\mathbb{Z}[\chi]$.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f)
Integer lattice of rank 5 and degree 5
  with isometry of finite order 5
  given by
  [1    1    1    1    1]
  [0   -1   -1   -1   -1]
  [0    1    0    0    0]
  [0    0    1    0    0]
  [0    0    0    1    0]

julia> is_of_hermitian_type(Lf)
false

julia> is_of_hermitian_type(coinvariant_lattice(Lf))
true
```
"""
function is_of_hermitian_type(Lf::ZZLatWithIsom)
  @req rank(Lf) > 0 "Underlying lattice must have positive rank"
  return is_irreducible(minimal_polynomial(Lf))
end

@doc raw"""
    hermitian_structure(Lf::ZZLatWithIsom) -> HermLat

Given a lattice with isometry $(L, f)$ such that the minimal polynomial of the
underlying isometry `f` is irreducible, return the hermitian structure of the
underlying lattice `L` over the equation order of the minimal polynomial of
`f`.

If it exists, the hermitian structure is stored. For now, we only cover the case
where the equation order is maximal (which is always the case when the order is
finite, for instance, since the minimal polynomial is cyclotomic).

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> M = coinvariant_lattice(Lf)
Integer lattice of rank 4 and degree 5
  with isometry of finite order 5
  given by
  [-1   -1   -1   -1]
  [ 1    0    0    0]
  [ 0    1    0    0]
  [ 0    0    1    0]

julia> H = hermitian_structure(M)
Hermitian lattice of rank 1 and degree 1
  over relative maximal order of Relative number field of degree 2 over maximal real subfield of cyclotomic field of order 5
  with pseudo-basis 
  (1, 1//1 * <1, 1>)
  (z_5, 1//1 * <1, 1>)

julia> res = get_attribute(M, :transfer_data)
Map of change of scalars
  from quadratic space of dimension 4
  to hermitian space of dimension 1

julia> M2, f2 = trace_lattice_with_isometry(H, res)
(Integer lattice of rank 4 and degree 4, [-1 -1 -1 -1; 1 0 0 0; 0 1 0 0; 0 0 1 0])

julia> genus(M) == genus(M2) # One class in this genus, so they are isometric
true

julia> f2 == isometry(M)
true
```
"""
@attr HermLat function hermitian_structure(Lf::ZZLatWithIsom)
  @req is_of_hermitian_type(Lf) "Lf is not of hermitian type"
  f = isometry(Lf)

  H, res = hermitian_structure_with_transfer_data(lattice(Lf), f; ambient_representation = false)

  set_attribute!(Lf, :transfer_data, res)
  return H
end

###############################################################################
#
#  Image centralizer in discriminant group
#
###############################################################################

@doc raw"""
    discriminant_group(Lf::ZZLatWithIsom) -> TorQuadModule, AutomorphismGroupElem

Given an integral lattice with isometry $(L, f)$, return the discriminant group `q`
of the underlying lattice `L` as well as this image of the underlying isometry
`f` inside $O(q)$.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> qL, qf = discriminant_group(Lf)
(Finite quadratic module: Z/6 -> Q/2Z, [1])

julia> qL
Finite quadratic module
  over integer ring
Abelian group: Z/6
Bilinear value module: Q/Z
Quadratic value module: Q/2Z
Gram matrix quadratic form: 
[5//6]

julia> qf
Isometry of Finite quadratic module: Z/6 -> Q/2Z defined by 
[1]

julia> f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> discriminant_group(Lf)[2]
Isometry of Finite quadratic module: Z/6 -> Q/2Z defined by 
[5]
```
"""
function discriminant_group(Lf::ZZLatWithIsom)
  @req is_integral(Lf) "Underlying lattice must be integral"
  L = lattice(Lf)
  f = ambient_isometry(Lf)
  q = discriminant_group(L)
  Oq = orthogonal_group(q)
  return (q, Oq(gens(matrix_group(f))[1]; check = false))::Tuple{TorQuadModule, AutomorphismGroupElem{TorQuadModule}}
end

@doc raw"""
    image_centralizer_in_Oq(Lf::ZZLatWithIsom) -> AutomorphismGroup{TorQuadModule}

Given an integral lattice with isometry $(L, f)$, return the image $G_L$ in
$O(q_L, \bar{f})$ of the centralizer $O(L, f)$ of `f` in $O(L)$. Here $q_L$
denotes the discriminant group of `L` and $\bar{f}$ is the isometry of
$q_L$ induced by `f`.

# Examples
```jldoctest
julia> L = root_lattice(:A,2);

julia> f = matrix(QQ, 2, 2, [1 1; 0 -1]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> G = image_centralizer_in_Oq(Lf)
Group of isometries of Finite quadratic module: Z/3 -> Q/2Z generated by 2 elements

julia> order(G)
2
```
"""
@attr AutomorphismGroup{TorQuadModule} function image_centralizer_in_Oq(Lf::ZZLatWithIsom)
  @req is_integral(Lf) "Underlying lattice must be integral"
  n = order_of_isometry(Lf)
  L = lattice(Lf)
  f = ambient_isometry(Lf)
  if (n in [1, -1]) || (isometry(Lf) == -identity_matrix(QQ, rank(L)))
    GL, _ = image_in_Oq(L)
  elseif is_definite(L) 
    OL = orthogonal_group(L)
    f = OL(f)
    UL = QQMatrix[matrix(OL(s)) for s in gens(centralizer(OL, f)[1])]
    qL = discriminant_group(L)
    UL = ZZMatrix[matrix(hom(qL, qL, elem_type(qL)[qL(lift(t)*g) for t in gens(qL)])) for g in UL]
    unique!(UL)
    GL = Oscar._orthogonal_group(qL, UL; check = false)
  elseif rank(L) == euler_phi(n)
    qL = discriminant_group(L)
    UL = ZZMatrix[matrix(hom(qL, qL, elem_type(qL)[qL(-lift(t)) for t in gens(qL)]))]
    unique!(UL)
    GL = Oscar._orthogonal_group(qL, UL; check = false)
  else
    @req is_of_hermitian_type(Lf) "Not yet implemented for indefinite lattices with isometry which are not of hermitian type"
    dets = Oscar._local_determinants_morphism(Lf)
    GL, _ = kernel(dets)
  end
  return GL::AutomorphismGroup{TorQuadModule}
end

###############################################################################
#
#  Signatures
#
###############################################################################

function _real_kernel_signatures(L::ZZLat, M) 
  C = base_ring(M)
  bL = basis_matrix(L)
  GL = gram_matrix(ambient_space(L))
  bLC = change_base_ring(C, bL)
  GLC = change_base_ring(C, GL)
  k, KC = left_kernel(M)
  newGC = KC*bLC*GLC*transpose(KC*bLC)

  newGC = Hecke._gram_schmidt(newGC, C)[1]
  diagC = diagonal(newGC)

  @hassert :ZZLatWithIsom 1 all(z -> isreal(z), diagC)
  @hassert :ZZLatWithIsom 1 all(z -> !iszero(z), diagC)

  k1 = count(z -> z > 0, diagC)
  k2 = length(diagC) - k1

  return k1, k2
end

@doc raw"""
    signatures(Lf::ZZLatWithIsom) -> Dict{Int, Tuple{Int, Int}}

Given a lattice with isometry $(L, f)$ where the minimal polynomial of `f`
is irreducible cyclotomic, return the signatures of $(L, f)$.

In this context, if we denote $z$ a primitive `n`-th root of unity, where `n`
is the order of `f`, then for each $1 \leq i \leq n/2$ such that $(i, n) = 1$,
the $i$-th signature of $(L, f)$ is given by the signatures of the real quadratic
form $\ker(f + f^{-1} - z^i - z^{-i})$.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> M = coinvariant_lattice(Lf);

julia> signatures(M)
Dict{Integer, Tuple{Int64, Int64}} with 2 entries:
  2 => (2, 0)
  1 => (2, 0)
```
"""
function signatures(Lf::ZZLatWithIsom)
  @req is_cyclotomic_polynomial(minimal_polynomial(Lf)) "Lf must be of finite hermitian type"
  L = lattice(Lf)
  f = isometry(Lf)
  n = order_of_isometry(Lf)
  C = CalciumField()
  eig = eigenvalues(f, QQBar)
  j = findfirst(z -> findfirst(k -> isone(z^k), 1:n) == n, eig)
  lambda = C(eig[j])
  Sq = [i for i in 1:div(n,2) if gcd(i,n) == 1]
  D = Dict{Integer, Tuple{Int64, Int64}}()
  fC = change_base_ring(C, f)
  for i in Sq
    M = fC + inv(fC) - lambda^i - lambda^(-i)
    D[i] = _real_kernel_signatures(L, M)
  end
  return D
end

###############################################################################
#
#  Kernel sublattices
#
###############################################################################

function _divides(k::IntExt, n::Int)
  is_finite(k) && return is_divisible_by(k, n)
  return true
end

@doc raw"""
    kernel_lattice(Lf::ZZLatWithIsom, p::Union{fmpz_poly, QQPolyRingElem})
                                                         -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$ and a polynomial `p` with rational
coefficients, return the sublattice $M := \ker(p(f))$ of the underlying lattice
`L` with isometry `f`, together with the restriction $f_{\mid M}$.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> Zx,x = ZZ["x"]
(Univariate polynomial ring in x over ZZ, x)

julia> mf = minimal_polynomial(Lf)
x^5 - 1

julia> factor(mf)
1 * (x - 1) * (x^4 + x^3 + x^2 + x + 1)

julia> kernel_lattice(Lf, x-1)
Integer lattice of rank 1 and degree 5
  with isometry of finite order 1
  given by
  [1]

julia> kernel_lattice(Lf, cyclotomic_polynomial(5))
Integer lattice of rank 4 and degree 5
  with isometry of finite order 5
  given by
  [-1   -1   -1   -1]
  [ 1    0    0    0]
  [ 0    1    0    0]
  [ 0    0    1    0]
```
"""
kernel_lattice(Lf::ZZLatWithIsom, p::Union{ZZPolyRingElem, QQPolyRingElem})

function kernel_lattice(Lf::ZZLatWithIsom, p::QQPolyRingElem)
  n = order_of_isometry(Lf)
  L = lattice(Lf)
  f = isometry(Lf)
  M = p(f)
  d = denominator(M)
  k, K = left_kernel(change_base_ring(ZZ, d*M))
  return lattice(ambient_space(Lf), K*basis_matrix(L))
end

kernel_lattice(Lf::ZZLatWithIsom, p::ZZPolyRingElem) = kernel_lattice(Lf, change_base_ring(QQ, p))

@doc raw"""
    kernel_lattice(Lf::ZZLatWithIsom, l::Integer) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$ and an integer `l`, return the kernel
lattice of $(L, f)$ associated to the `l`-th cyclotomic polynomial.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> kernel_lattice(Lf, 5)
Integer lattice of rank 4 and degree 5
  with isometry of finite order 5
  given by
  [-1   -1   -1   -1]
  [ 1    0    0    0]
  [ 0    1    0    0]
  [ 0    0    1    0]

julia> kernel_lattice(Lf, 1)
Integer lattice of rank 1 and degree 5
  with isometry of finite order 1
  given by
  [1]
```
"""
function kernel_lattice(Lf::ZZLatWithIsom, l::Integer)
  @req _divides(order_of_isometry(Lf), l)[1] "l must divide the order of the underlying isometry"
  p = cyclotomic_polynomial(l)
  return kernel_lattice(Lf, p)
end

@doc raw"""
    invariant_lattice(Lf::ZZLatWithIsom) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$, return the invariant lattice $L^f$ of
$(L, f)$ together with the restriction of `f` to $L^f$ (which is the identity
in this case).

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> invariant_lattice(Lf)
Integer lattice of rank 1 and degree 5
  with isometry of finite order 1
  given by
  [1]
```
"""
invariant_lattice(Lf::ZZLatWithIsom) = kernel_lattice(Lf, 1)

@doc raw"""
    invariant_lattice(L::ZZLat, G::MatrixGroup;
                                ambient_representation::Bool = true) -> ZZLat

Given an integer lattice `L` and a group `G` of isometries of `L` in matrix,
return the invariant sublattice $L^G$ of `L`.

If `ambient_representation` is set to true, the isometries in `G` are seen as
isometries of the ambient quadratic space of `L` preserving `L`. Otherwise, they
are considered as honnest isometries of `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,2);

julia> G = isometry_group(L);

julia> invariant_lattice(L, G)
Integer lattice of rank 0 and degree 2
with gram matrix
0 by 0 empty matrix
```
"""
function invariant_lattice(L::ZZLat, G::MatrixGroup; ambient_representation::Bool = true)
  return invariant_lattice(L, matrix.(gens(G)); ambient_representation = ambient_representation)
end

@doc raw"""
    coinvariant_lattice(Lf::ZZLatWithIsom) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$, return the coinvariant lattice $L_f$ of
$(L, f)$ together with the restriction of `f` to $L_f$.

The coinvariant lattice $L_f$ of $(L, f)$ is the orthogonal complement in
`L` of the invariant lattice $L_f$.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> coinvariant_lattice(Lf)
Integer lattice of rank 4 and degree 5
  with isometry of finite order 5
  given by
  [-1   -1   -1   -1]
  [ 1    0    0    0]
  [ 0    1    0    0]
  [ 0    0    1    0]
```
"""
function coinvariant_lattice(Lf::ZZLatWithIsom)
  chi = minimal_polynomial(Lf)
  if chi(1) == 0
    R = parent(chi)
    x = gen(R)
    chi = divexact(chi, x-1)
  end
  return kernel_lattice(Lf, chi)
end

@doc raw"""
    coinvariant_lattice(L::ZZLat, G::MatrixGroup;
                                  ambient_representation::Bool = true)
                                                        -> ZZLat, MatrixGroup

Given an integer lattice `L` and a group `G` of isometries of `L` in matrix,
return the coinvariant sublattice $L_G$ of `L`, together with the subgroup `H`
of isometries of $L_G$ induced by the action of $G$.

If `ambient_representation` is set to true, the isometries in `G` and `H` are seen
as isometries of the ambient quadratic space of `L` preserving `L`. Otherwise,
they are considered as honnest isometries of `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,2);

julia> G = isometry_group(L);

julia> L2, G2 = coinvariant_lattice(L, G)
(Integer lattice of rank 2 and degree 2, Matrix group of degree 2 over Rational field)

julia> L == L2
true

julia> G == G2
true
```
"""
function coinvariant_lattice(L::ZZLat, G::MatrixGroup; ambient_representation::Bool = true)
  F = invariant_lattice(L, G; ambient_representation = ambient_representation)
  C = orthogonal_submodule(L, F)
  gene = QQMatrix[]
  for g in gens(G)
    if !ambient_representation
      mL = matrix(g)
      m_amb = solve(basis_matrix(L), mL*basis_matrix(L))
      mC = solve_left(basis_matrix(C), basis_matrix(C)*m_amb)
      push!(gene, mC)
    else
      push!(gene, matrix(g))
    end
  end
  return C, matrix_group(gene)
end

@doc raw"""
    invariant_coinvariant_pair(L::ZZLat, G::MatrixGroup;
                                         ambient_representation::Bool = true)
                                                   -> ZZLat, ZZLat, MatrixGroup

Given an integer lattice `L` and a group `G` of isometries of `L` in matrix,
return the invariant sublattice $L^G$ of `L` and its coinvariant sublattice
$L_G$ together with the subgroup `H` of isometries of $L_G$ induced by the
action of $G$.

If `ambient_representation` is set to true, the isometries in `G` and `H` are seen
as isometries of the ambient quadratic space of `L` preserving `L`. Otherwise,
they are considered as honnest isometries of `L`.

# Examples
```jldoctest
julia> L = root_lattice(:A,2);

julia> G = isometry_group(L);

julia> Gsub, _ = sub(G, [gens(G)[end]]);

julia> F, C, G2 = invariant_coinvariant_pair(L, Gsub)
(Integer lattice of rank 1 and degree 2, Integer lattice of rank 1 and degree 2, Matrix group of degree 2 over Rational field)

julia> F
Integer lattice of rank 1 and degree 2
with gram matrix
[2]

julia> C
Integer lattice of rank 1 and degree 2
with gram matrix
[6]
```
"""
function invariant_coinvariant_pair(L::ZZLat, G::MatrixGroup; ambient_representation::Bool = true)
  F = invariant_lattice(L, G; ambient_representation = ambient_representation)
  C = orthogonal_submodule(L, F)
  gene = QQMatrix[]
  for g in gens(G)
    if !ambient_representation
      mL = matrix(g)
      m_amb = solve(basis_matrix(L), mL*basis_matrix(L))
      mC = solve_left(basis_matrix(C), basis_matrix(C)*m_amb)
      push!(gene, mC)
    else
      push!(gene, matrix(g))
    end
  end
  return F, C, matrix_group(gene)
end

##############################################################################
#
#  Type
#
###############################################################################

@doc raw"""
    type(Lf::ZZLatWithIsom)
                      -> Dict{Int, Tuple{ <: Union{ZZGenus, HermGenus}, ZZGenus}}

Given a lattice with isometry $(L, f)$ with `f` of finite order `n`, return the
type of $(L, f)$.

In this context, the type is defined as follows: for each divisor `k` of `n`,
the `k`-type of $(L, f)$ is the tuple $(H_k, A_K)$ consisting of the genus
$H_k$ of the lattice $\ker(\Phi_k(f))$ viewed as a hermitian $\mathbb{Z}[\zeta_k]$-
lattice (so a $\mathbb{Z}$-lattice for k= 1, 2) and of the genus $A_k$ of the
$\mathbb{Z}$-lattice $\ker(f^k-1)$.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> t = type(Lf);

julia> genus(invariant_lattice(Lf)) == t[1][1]
true
```
"""
@attr function type(Lf::ZZLatWithIsom)
  L = lattice(Lf)
  f = isometry(Lf)
  n = order_of_isometry(Lf)
  @req is_finite(n) "Isometry must be of finite order"
  divs = divisors(n)
  Qx = Hecke.Globals.Qx
  x = gen(Qx)
  t = Dict{Integer, Tuple}()
  for l in divs
    Hl = kernel_lattice(Lf, cyclotomic_polynomial(l))
    if !(order_of_isometry(Hl) in [-1,1,2])
      Hl = Hecke.hermitian_structure(lattice(Hl), isometry(Hl); check=false, ambient_representation=false)
    end
    Al = kernel_lattice(Lf, x^l-1)
    t[l] = (genus(Hl), genus(Al))
  end
  return t
end

@doc raw"""
    is_of_type(Lf::ZZLatWithIsom, t::Dict) -> Bool

Given a lattice with isometry $(L, f)$, return whether $(L, f)$ is of type `t`.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> t = type(Lf);

julia> is_of_type(Lf, t)
true
```
"""
function is_of_type(L::ZZLatWithIsom, t::Dict)
  @req is_finite(order_of_isometry(L)) "Type is defined only for finite order isometries"
  divs = sort(collect(keys(t)))
  x = gen(Hecke.Globals.Qx)
  for l in divs
    Hl = kernel_lattice(L, cyclotomic_polynomial(l))
    if !(order_of_isometry(Hl) in [-1, 1, 2])
      t[l][1] isa Hecke.HermGenus || return false
      Hl = Hecke.hermitian_structure(lattice(Hl), isometry(Hl); check=false, ambient_representation=false, E = base_field(t[l][1]))
    end
    genus(Hl) == t[l][1] || return false
    Al = kernel_lattice(L, x^l-1)
    genus(Al) == t[l][2] || return false
  end
  return true
end

@doc raw"""
    is_of_same_type(Lf::ZZLatWithIsom, Mg::ZZLatWithIsom) -> Bool

Given two lattices with isometry $(L, f)$ and $(M, g)$, return whether they are
of the same type.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> M = coinvariant_lattice(Lf);

julia> is_of_same_type(Lf, M)
false
```
"""
function is_of_same_type(L::ZZLatWithIsom, M::ZZLatWithIsom)
  @req is_finite(order_of_isometry(L)*order_of_isometry(M)) "Type is defined only for finite order isometries"
  order_of_isometry(L) != order_of_isometry(M) && return false
  genus(L) != genus(M) && return false
  return is_of_type(L, type(M))
end

@doc raw"""
    is_hermitian(t::Dict) -> Bool

Given a type `t` of lattices with isometry, return whether `t` is hermitian, i.e.
whether it defines the type of a hermitian lattice with isometry.

# Examples
```jldoctest
julia> L = root_lattice(:A,5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> M = coinvariant_lattice(Lf);

julia> is_hermitian(type(Lf))
false

julia> is_hermitian(type(M))
true
```
"""
function is_hermitian(t::Dict)
  ke = collect(keys(t))
  n = maximum(ke)
  return all(i -> rank(t[i][1]) == rank(t[i][2]) == 0, [i for i in ke if i != n])
end

###############################################################################
#
#  Useful
#
###############################################################################

function to_oscar(io::IO, Lf::ZZLatWithIsom)
  L = lattice(Lf)
  f = ambient_isometry(Lf)
  println(io, "B = matrix(QQ, $(rank(L)), $(degree(L)), " , basis_matrix(L), ");")
  println(io, "G = matrix(QQ, $(degree(L)), $(degree(L)), ", gram_matrix(ambient_space(L)), ");")
  println(io, "L = integer_lattice(B, gram = G);")
  println(io, "f = matrix(QQ, $(degree(L)), $(degree(L)), ", f, ");")
  println(io, "Lf = integer_lattice_with_isometry(L, f);")
end

to_oscar(Lf::ZZLatWithIsom) = to_oscar(stdout, Lf)

