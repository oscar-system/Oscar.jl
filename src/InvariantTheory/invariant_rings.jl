export invariant_ring, affine_algebra
export coefficient_ring, polynomial_ring, action, group
export is_modular
export reynolds_operator, molien_series

################################################################################
#
#  Field access
#
################################################################################

coefficient_ring(I::InvRing) = I.field

polynomial_ring(I::InvRing) = I.poly_ring

action(I::InvRing) = I.action

group(I::InvRing) = I.group

is_modular(I::InvRing) = I.modular

function invariant_ring(M::Vector{<: MatrixElem})
  return invariant_ring(base_ring(M[1]), M)
end

invariant_ring(matrices::MatrixElem{T}...) where {T} = invariant_ring(collect(matrices))

function invariant_ring(K::Field, M::Vector{<: MatrixElem})
  return invariant_ring(matrix_group([change_base_ring(K, g) for g in M]))
end

@doc Markdown.doc"""
    invariant_ring(G::MatrixGroup)
    invariant_ring(K::Field = QQ, G::PermGroup)

Return the invariant ring of the finite matrix group or permutation group `G`.

In the latter case, use the specified field `K` as the coefficient field. 
The default value for `K` is `QQ`.

!!! note
    The creation of invariant rings is lazy in the sense that no explicit computations are done until specifically invoked (for example, by the `primary_invariants` function).

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a");

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0]);

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1]);

julia> G = MatrixGroup(3, K, [M1, M2]);

julia> IRm = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> IRp = invariant_ring(symmetric_group(3))
Invariant ring of
Sym( [ 1 .. 3 ] )
with generators
PermGroupElem[(1,2,3), (1,2)]

julia> coefficient_ring(IRp)
Rational Field
```
"""
function invariant_ring(G::MatrixGroup)
  action = mat_elem_type(typeof(G))[g.elm for g in gens(G)]
  return InvRing(base_ring(G), G, action)
end

invariant_ring(K::Field, G::PermGroup) = InvRing(K, G, gens(G))

invariant_ring(G::PermGroup) = invariant_ring(QQ, G)

function Base.show(io::IO, IR::InvRing)
  print(io, "Invariant ring of\n")
  print(io, group(IR), "\n")
  print(io, "with generators\n")
  print(io, action(IR))
end

# Return a map performing the right action of M on the ring R
function right_action(R::MPolyRing{T}, M::MatrixElem{T}) where T
  @assert nvars(R) == ncols(M)
  @assert nrows(M) == ncols(M)
  n = nvars(R)

  vars = zeros(R, n)
  x = gens(R)
  for i = 1:n
    for j = 1:n
      if iszero(M[i, j])
        continue
      end
      vars[i] = addeq!(vars[i], M[i, j]*x[j])
    end
  end

  right_action_by_M = (f::MPolyElem{T}) -> evaluate(f, vars)

  return MapFromFunc(right_action_by_M, R, R)
end

right_action(R::MPolyRing{T}, M::MatrixGroupElem{T}) where T = right_action(R, M.elm)
right_action(f::MPolyElem{T}, M::MatrixElem{T}) where T = right_action(parent(f), M)(f)
right_action(f::MPolyElem{T}, M::MatrixGroupElem{T}) where T = right_action(f, M.elm)

function right_action(R::MPolyRing{T}, p::PermGroupElem) where T
  n = nvars(R)
  @assert n == degree(parent(p))

  right_action_by_p = (f::MPolyElem{T}) -> on_indeterminates(f, p)

  return MapFromFunc(right_action_by_p, R, R)
end

right_action(f::MPolyElem, p::PermGroupElem) = right_action(parent(f), p)(f)

################################################################################
#
#  Reynolds operator
#
################################################################################

function reynolds_operator(IR::InvRing{FldT, GrpT, PolyElemT}) where {FldT, GrpT, PolyElemT}
  @assert !is_modular(IR)

  if isdefined(IR, :reynolds_operator)
    return nothing
  end

  actions = [ right_action(polynomial_ring(IR), g) for g in group(IR) ]
  function reynolds(f::PolyElemT)
    g = parent(f)()
    for action in actions
      g = addeq!(g, action(f))
    end
    return g*base_ring(f)(1//order(group(IR)))
  end

  IR.reynolds_operator = MapFromFunc(reynolds, polynomial_ring(IR), polynomial_ring(IR))
  return IR.reynolds_operator
end

@doc Markdown.doc"""
     reynolds_operator(IR::InvRing{FldT, GrpT, T}, f::T) where {FldT, GrpT, T <: MPolyElem}

In the non-modular case, return the image of `f` under the Reynolds operator
projecting onto `IR`.

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = MatrixGroup(3, K, [ M1, M2 ])
Matrix group of degree 3 over Cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> R = polynomial_ring(IR)
Multivariate Polynomial Ring in x[1], x[2], x[3] over Cyclotomic field of order 3 graded by
  x[1] -> [1]
  x[2] -> [1]
  x[3] -> [1]

julia> x = gens(R)
3-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 x[1]
 x[2]
 x[3]

julia> f = x[1]^3
x[1]^3

julia> reynolds_operator(IR, f)
1//3*x[1]^3 + 1//3*x[2]^3 + 1//3*x[3]^3

julia> M = matrix(GF(3), [0 1 0; -1 0 0; 0 0 -1])
[0   1   0]
[2   0   0]
[0   0   2]

julia> G = MatrixGroup(3, GF(3), [M])
Matrix group of degree 3 over Galois field with characteristic 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Galois field with characteristic 3
with generators
gfp_mat[[0 1 0; 2 0 0; 0 0 2]]

julia> R = polynomial_ring(IR)
Multivariate Polynomial Ring in x[1], x[2], x[3] over Galois field with characteristic 3 graded by
  x[1] -> [1]
  x[2] -> [1]
  x[3] -> [1]

julia> x = gens(R)
3-element Vector{MPolyElem_dec{gfp_elem, gfp_mpoly}}:
 x[1]
 x[2]
 x[3]

julia> f = x[1]^2
x[1]^2

julia> reynolds_operator(IR, f)
2*x[1]^2 + 2*x[2]^2

julia> f = x[1]^3
x[1]^3

julia> reynolds_operator(IR, f)
0
```
"""
function reynolds_operator(IR::InvRing{FldT, GrpT, T}, f::T) where {FldT, GrpT, T <: MPolyElem}
  @assert !is_modular(IR)
  @assert parent(f) === polynomial_ring(IR)

  if !isdefined(IR, :reynolds_operator)
    reynolds_operator(IR)
  end
  return IR.reynolds_operator(f)
end

function reynolds_operator(IR::InvRing, f::MPolyElem)
  @assert parent(f) === polynomial_ring(IR).R
  return reynolds_operator(IR, polynomial_ring(IR)(f))
end

function reynolds_operator(IR::InvRing{FldT, GrpT, PolyElemT}, chi::GAPGroupClassFunction) where {FldT, GrpT, PolyElemT}
# I expect that this also works in the non-modular case, but haven't found a reference.
  # The only reference for this version of the reynolds operator appears to be [Gat96].
  @assert is_zero(characteristic(coefficient_ring(IR)))
  @assert is_irreducible(chi)

  K = coefficient_ring(IR)

  conjchi = conj(chi)

  actions_and_values = [ (right_action(polynomial_ring(IR), g), K(conjchi(g).data)) for g in group(IR) ]

  function reynolds(f::PolyElemT)
    g = parent(f)()
    for (action, val) in actions_and_values
      g = addeq!(g, action(f)*val)
    end
    return g*base_ring(f)(1//order(group(IR)))
  end

  return MapFromFunc(reynolds, polynomial_ring(IR), polynomial_ring(IR))
end

@doc Markdown.doc"""
     reynolds_operator(IR::InvRing{FldT, GrpT, T}, f::T, chi::GAPGroupClassFunction)
       where {FldT, GrpT, T <: MPolyElem}

In the case of characteristic zero, return the image of `f` under the twisted
Reynolds operator projecting onto the isotypic component of the polynomial ring
with respect to `chi`, that is, the semi-invariants (or relative invariants)
with respect to `chi`, see [Sta79](@cite).
It is assumed that `chi` is an irreducible character.

In case `chi` is a linear character, the returned polynomial, say `h`, fulfils
`h^g = chi(g)h` for all `g` in `group(IR)` (possibly `h` is zero).

!!! note
    It is implicitly assumed that `coefficient_ring(IR)` contains all character values of `chi`.

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a");

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0]);

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1]);

julia> G = MatrixGroup(3, K, [ M1, M2 ]);

julia> IR = invariant_ring(G);

julia> R = polynomial_ring(IR);

julia> x = gens(R);

julia> f = x[1]^3
x[1]^3

julia> reynolds_operator(IR, f, trivial_character(G))
1//3*x[1]^3 + 1//3*x[2]^3 + 1//3*x[3]^3

julia> S2 = symmetric_group(2);

julia> IR = invariant_ring(QQ, S2);

julia> R = polynomial_ring(IR);

julia> x = gens(R);

julia> F = abelian_closure(QQ)[1];

julia> chi = Oscar.group_class_function(S2, [ F(sign(representative(c))) for c in conjugacy_classes(S2) ])
group_class_function(character_table(Sym( [ 1 .. 2 ] )), QQAbElem{nf_elem}[1, -1])

julia> reynolds_operator(IR, x[1], chi)
1//2*x[1] - 1//2*x[2]
```
"""
function reynolds_operator(IR::InvRing{FldT, GrpT, T}, f::T, chi::GAPGroupClassFunction) where {FldT, GrpT, T <: MPolyElem}
  return reynolds_operator(IR, chi)(f)
end

function reynolds_operator(IR::InvRing, f::MPolyElem, chi::GAPGroupClassFunction)
  @assert parent(f) === polynomial_ring(IR).R
  return reynolds_operator(IR, polynomial_ring(IR)(f), chi)
end

################################################################################
#
#  Bases
#
################################################################################

@doc Markdown.doc"""
     basis(IR::InvRing, d::Int, algo::Symbol = :default)

Given an invariant ring `IR` and an integer `d`, return a basis for the invariants in degree `d`.

The optional argument `algo` specifies the algorithm to be used.
If `algo = :reynolds`, the Reynolds operator is utilized (this method is only available in the non-modular case).
Setting `algo = :linear_algebra` means that plain linear algebra is used.
The default option `algo = :default` asks to select the heuristically best algorithm.

See also [`iterate_basis`](@ref).

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = MatrixGroup(3, K, [ M1, M2 ])
Matrix group of degree 3 over Cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> basis(IR, 6)
4-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 x[1]^2*x[2]^2*x[3]^2
 x[1]^4*x[2]*x[3] + x[1]*x[2]^4*x[3] + x[1]*x[2]*x[3]^4
 x[1]^3*x[2]^3 + x[1]^3*x[3]^3 + x[2]^3*x[3]^3
 x[1]^6 + x[2]^6 + x[3]^6

julia> M = matrix(GF(3), [0 1 0; -1 0 0; 0 0 -1])
[0   1   0]
[2   0   0]
[0   0   2]

julia> G = MatrixGroup(3, GF(3), [M])
Matrix group of degree 3 over Galois field with characteristic 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Galois field with characteristic 3
with generators
gfp_mat[[0 1 0; 2 0 0; 0 0 2]]

julia> basis(IR, 2)
2-element Vector{MPolyElem_dec{gfp_elem, gfp_mpoly}}:
 x[1]^2 + x[2]^2
 x[3]^2

julia> basis(IR, 3)
2-element Vector{MPolyElem_dec{gfp_elem, gfp_mpoly}}:
 x[1]*x[2]*x[3]
 x[1]^2*x[3] + 2*x[2]^2*x[3]
```
"""
basis(IR::InvRing, d::Int, algo = :default) = collect(iterate_basis(IR, d, algo))

@doc Markdown.doc"""
    basis(IR::InvRing, d::Int, chi::GAPGroupClassFunction)

Given an invariant ring `IR`, an integer `d` and an irreducible character `chi`,
return a basis for the semi-invariants (or relative invariants) in degree `d`
with respect to `chi`.

This function is only implemented in the case of characteristic zero.

!!! note
    It is implicitly assumed that `coefficient_ring(IR)` contains all character values of `chi`.

See also [`iterate_basis`](@ref).

# Examples
```
julia> K, a = CyclotomicField(3, "a");

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0]);

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1]);

julia> G = MatrixGroup(3, K, [ M1, M2 ]);

julia> IR = invariant_ring(G);

julia> basis(IR, 6, trivial_character(G))
4-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 x[1]^6 + x[2]^6 + x[3]^6
 x[1]^4*x[2]*x[3] + x[1]*x[2]^4*x[3] + x[1]*x[2]*x[3]^4
 x[1]^3*x[2]^3 + x[1]^3*x[3]^3 + x[2]^3*x[3]^3
 x[1]^2*x[2]^2*x[3]^2

julia> S2 = symmetric_group(2);

julia> R = invariant_ring(QQ, S2);

julia> F = abelian_closure(QQ)[1];

julia> chi = Oscar.group_class_function(S2, [ F(sign(representative(c))) for c in conjugacy_classes(S2) ])
group_class_function(character_table(Sym( [ 1 .. 2 ] )), QQAbElem{nf_elem}[1, -1])

julia> basis(R, 3, chi)
2-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x[1]^3 - x[2]^3
 x[1]^2*x[2] - x[1]*x[2]^2

```
"""
basis(IR::InvRing, d::Int, chi::GAPGroupClassFunction) = collect(iterate_basis(IR, d, chi))

################################################################################
#
#  Primary invariants
#
################################################################################

# See primary_invariants.jl

################################################################################
#
#  Secondary invariants
#
################################################################################

# See secondary_invariants.jl

################################################################################
#
#  Molien series
#
################################################################################

function _molien_series_char0(S::PolyRing, I::InvRing)
  G = group(I)
  n = degree(G)
  if G isa MatrixGroup{T, T1} where T1<:MatElem{T} where T<:Union{fmpq, fmpz, nf_elem}
    Gp, GtoGp = isomorphic_group_over_finite_field(G)
  else
    Gp, GtoGp = (G, id_hom(G))
  end
  K = coefficient_ring(I)
  Kt, _ = PolynomialRing(K, "t", cached = false)
  C = conjugacy_classes(Gp)
  res = zero(FractionField(Kt))
  for c in C
    g = (GtoGp\(representative(c)))::elem_type(G)
    if g isa MatrixGroupElem
      f = charpoly(Kt, g.elm)
    elseif g isa PermGroupElem
      f = charpoly(Kt, permutation_matrix(K, g))
    else
      error("problem to compute the char. pol. of $g")
    end
    res = res + length(c)::fmpz * 1//reverse(f)
  end
  res = divexact(res, order(Gp)::fmpz)
  num = change_coefficient_ring(coefficient_ring(S),
                                numerator(res), parent = S)
  den = change_coefficient_ring(coefficient_ring(S),
                                denominator(res), parent = S)
  return num//den
end

function _molien_series_nonmodular_via_gap(S::PolyRing, I::InvRing, chi::Union{GAPGroupClassFunction, Nothing} = nothing)
  @assert !is_modular(I)
  G = group(I)
  @assert G isa MatrixGroup || G isa PermGroup
  t = GAP.Globals.CharacterTable(G.X)
  if G isa MatrixGroup
    if is_zero(characteristic(coefficient_ring(I)))
      psi = natural_character(G).values
    else
      psi = [GAP.Globals.BrauerCharacterValue(GAP.Globals.Representative(c))
             for c in GAP.Globals.ConjugacyClasses(t)]
    end
  else
    deg = GAP.Obj(degree(G))
    psi = [deg - GAP.Globals.NrMovedPoints(GAP.Globals.Representative(c))
           for c in GAP.Globals.ConjugacyClasses(t)]
  end
  if chi === nothing
    info = GAP.Globals.MolienSeriesInfo(GAP.Globals.MolienSeries(t,
                                                                 GAP.GapObj(psi)))
  else
    info = GAP.Globals.MolienSeriesInfo(GAP.Globals.MolienSeries(t,
                                                                 GAP.GapObj(psi), chi.values))
  end
  num = S(Vector{fmpz}(GAP.Globals.CoefficientsOfUnivariatePolynomial(info.numer))::Vector{fmpz})
  den = S(Vector{fmpz}(GAP.Globals.CoefficientsOfUnivariatePolynomial(info.denom))::Vector{fmpz})
  return num//den
end

@doc Markdown.doc"""
    molien_series([S::PolyRing], I::InvRing, [chi::GAPGroupClassFunction])

In the non-modular case, return the Molien series of `I` as a rational function.

If a univariate polynomial ring with rational coefficients is specified by the
optional argument `S::PolyRing`, then return the Molien series as an element
of the fraction field of that ring.

If a character `chi` is specified, the series relative to `chi` is returned.
This is the Molien series of the module of semi-invariants (or relative invariants)
with respect to `chi`, see [Sta79](@cite).

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a");

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0]);

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1]);

julia> G = MatrixGroup(3, K, [M1, M2]);

julia> IR = invariant_ring(G);

julia> MS = molien_series(IR)
(-t^6 + t^3 - 1)//(t^9 - 3*t^6 + 3*t^3 - 1)

julia> parent(MS)
Fraction field of Univariate Polynomial Ring in t over Rational Field

julia> S2 = symmetric_group(2);

julia> IR = invariant_ring(QQ, S2);

julia> F = abelian_closure(QQ)[1];

julia> chi = Oscar.group_class_function(S2, [ F(sign(representative(c))) for c in conjugacy_classes(S2) ])
group_class_function(character_table(Sym( [ 1 .. 2 ] )), QQAbElem{nf_elem}[1, -1])

julia> molien_series(IR)
1//(t^3 - t^2 - t + 1)

julia> molien_series(IR, chi)
t//(t^3 - t^2 - t + 1)
```
"""
function molien_series(S::PolyRing, I::InvRing, chi::Union{GAPGroupClassFunction, Nothing} = nothing)
  if isdefined(I, :molien_series) && chi === nothing
    if parent(I.molien_series) === S
      return I.molien_series
    end
    num = change_coefficient_ring(coefficient_ring(S), numerator(I.molien_series), parent = S)
    den = change_coefficient_ring(coefficient_ring(S), denominator(I.molien_series), parent = S)
    return num//den
  end

  if characteristic(coefficient_ring(I)) == 0 && chi === nothing
    return _molien_series_char0(S, I)
  else
    if !is_modular(I)
      return _molien_series_nonmodular_via_gap(S, I, chi)
    else
      throw(Hecke.NotImplemented())
    end
  end
end

function molien_series(I::InvRing)
  if !isdefined(I, :molien_series)
    S, t = PolynomialRing(QQ, "t", cached = false)
    I.molien_series = molien_series(S, I)
  end
  return I.molien_series
end

function molien_series(I::InvRing, chi::GAPGroupClassFunction)
  S, t = PolynomialRing(QQ, "t", cached = false)
  return molien_series(S, I, chi)
end

# There are some situations where one needs to know whether one can ask for the
# Molien series without throwing an error.
# And maybe some day we can also compute Molien series in some modular cases.
is_molien_series_implemented(I::InvRing) = !is_modular(I)

################################################################################
#
#  Fundamental invariants
#
################################################################################

# See fundamental_invariants.jl

################################################################################
#
#  Presentation as affine algebra
#
################################################################################

# See affine_algebra.jl
