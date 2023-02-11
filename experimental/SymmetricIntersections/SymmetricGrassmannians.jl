import AbstractAlgebra: is_empty

import Oscar: defining_ideal, describe, irreducible_components

export character_grassmannian,
       determinant_grassmannian,
       invariant_grassmannian,
       irreducible_component,
       isotypical_factor,
       isotypical_factors,
       isotypical_grassmannian,
       module_representation,
       parametrization_data,
       projective_dimension,
       standard_element,
       submodule_character,
       submodule_determinant_character,
       submodule_dimension

###############################################################################
#
# Accessors/Attributes
#
###############################################################################

@doc Markdown.doc"""
    is_irreducible(M::T) where T <: SymmetricGrassmannian -> Bool

Given a symmetric Grassmannian, return whether it is irreducible.
"""
is_irreducible(M::T) where T <: SymmetricGrassmannian

@doc Markdown.doc"""
    is_empty(M::T) where T <: SymmetricGrassmannian -> Bool

Given a symmetric Grassmannian, return whether it is empty.
"""
is_empty(M::T) where T <: SymmetricGrassmannian

@doc Markdown.doc"""
    projective_dimension(M::T) where T <: SymmetricGrassmannian -> Int

Given a symmetric Grassmannian, return its dimension as projective
variety.
"""
projective_dimension(M::T) where T <: SymmetricGrassmannian

@doc Markdown.doc"""
    submodule_dimension(M::T) where T <: SymmetricGrassmannian -> Int

Given a symmetric Grassmannian, return the dimension of the submodules
it parametrizes.
"""
submodule_dimension(M::T) where T <: SymmetricGrassmannian

@doc Markdown.doc"""
    module_representation(M::T) where T <: SymmetricGrassmannian -> LinRep

Given a symmetric Grassmannian parametrizing submodules of a module `V`,
return the representation mapping associated to `V`.
"""
module_representation(M::T) where T <: SymmetricGrassmannian

@doc Markdown.doc"""
    defining_ideal(M::T) where T <: SymmetricGrassmannian -> MPolyIdeal

Given a symmetric Grassmannian parametrizing submodules of a module `V`,
return the ideal defining its structure as projective variety.
"""
defining_ideal(M::T) where T <: SymmetricGrassmannian

@doc Markdown.doc"""
    describe(M::T) where T <: Union{IsotypicalGrassmannian,
                                    CharacterGrassmannian}  -> nothing
                           
Given a symmetric Grassmannian parametrizing submodules of a module `V`,
return a description of M in terms of Grassmannian varieties.
"""
describe(M::T) where T <: Union{IsotypicalGrassmannian, CharacterGrassmannian}

###############################################################################

### Isotypical Grassmannians

@doc Markdown.doc"""
    submodule_character(M::IsotypicalGrassmannian) -> Oscar.GAPGroupClassFunction

Given the isotypical Grassmannian variety `M` parametrizing isotypical submodules
with character `chi` of a module `V`, return `chi`.
"""
submodule_character(M::IsotypicalGrassmannian) = M.chi

submodule_dimension(M::IsotypicalGrassmannian) = Int(degree(M.chi))

module_representation(M::IsotypicalGrassmannian) = M.rep_mod

is_irreducible(M::IsotypicalGrassmannian) = true

is_empty(M::IsotypicalGrassmannian) = false

projective_dimension(M::IsotypicalGrassmannian) = M.t

function defining_ideal(M::IsotypicalGrassmannian)
  rep = module_representation(M)
  F = base_field(representation_ring(rep))
  chi = submodule_character(M)
  cd = character_decomposition(chi)[1]
  n = Int(scalar_product(character_representation(rep), cd[2]))
  t = cd[1]
  if t in [1, n-1, n]
    return defining_ideal(projective_space(F, binomial(n,t)-1))
  end
  t = cd[1]
  S, _ = GradedPolynomialRing(F, ["x[$j]" for j in 0:binomial(n, t)-1], [1 for i in 1:binomial(n,t)])
  return grassmann_pluecker_ideal(S, t, n)
end

function describe(M::IsotypicalGrassmannian)
  rep = module_representation(M)
  F = base_field(representation_ring(rep))
  chi = submodule_character(M)
  cd = character_decomposition(chi)[1]
  t = cd[1]
  n = Int(scalar_product(character_representation(rep), cd[2]))
  if n == t
    return "A point"
  elseif t == 1
    return "$(n-1)-dimensional projective space over F = $F"
  else
    return "Grassmannian Gr($t, F^$n) where F = $F"
  end
end

@doc Markdown.doc"""
    parametrization_data(M::IsotypicalGrassmannian)
                                        -> Tuple{Vector{MatSpaceElem}, Int}

Given a symmetric Grassmannian parametrizing submodules of a module `V` with
given isotypical character `chi`, return a tuple `(B, n)` where `B` is a basis
of the group `Hom(N, V)` of homomorphisms from `N` to `V`, where `N` is the simple
module of the representation ring `V` affording the simple character associated
to `chi` and `n` is the size of the submodules of `Hom(N, V)` parametrized by `M`.

The elements of `B` are given by matrices in the standard coordinates of `V` and `N`.
"""
function parametrization_data(M::IsotypicalGrassmannian{S, T, U}) where {S, T, U}
  f = M.vs_struct
  chi = submodule_character(M)
  n, _ = character_decomposition(chi)[1]
  B = f.(gens(domain(f)))::Vector{dense_matrix_type(U)}
  return B, n
end

@doc Markdown.doc"""
    standard_element(M::IsotypicalGrassmannian)
                                             -> Vector{MatSpaceElem}

Given a symmetric Grassmannian parametrizing submodules of a module `V` with
given isotypical character, return a specific elements parametrized by `M`
made from all possible elements in `M`.
"""
function standard_element(M::IsotypicalGrassmannian)
  B, n = parametrization_data(M)
  std_el = eltype(B)[]
  for i in 1:n
    push!(std_el, sum([B[j] for j in i:length(B)]))
  end
  return std_el
end

###############################################################################

### Character Grassmannians

@doc Markdown.doc"""
    isotypical_factors(M::CharacterGrassmannian{S, T, U})
                    where {S, T, U} -> Vector{IsotypicalGrassmannian{S, T, U}}

Given the character Grassmannian variety `M` parametrizing `t`-dimensional
submodules of a module `V` with given character, return the isotypical
Grassmannian which define `M` as product variety.
"""
isotypical_factors(M::CharacterGrassmannian) = M.constituents

@doc Markdown.doc"""
    isotypical_factor(M::CharacterGrassmannian{S, T, U},
                      chi::Oscar.GAPGroupClassFunction)
                           where {S, T, U} -> IsotypicalGrassmannian{S, T, U}

Given the character Grassmannian variety `M` parametrizing `t`-dimensional
submodules of a module `V` with given character, return the isotypical
Grassmannian in the factors of `M` corresponding to the isotypical character
`chi`.
"""
function isotypical_factor(M::CharacterGrassmannian{S, T, U},
                              chi::Oscar.GAPGroupClassFunction) where {S, T, U}
  @req is_isotypical(chi) "chi must be an isotypical character"
  @req chi in character_decomposition(M) "chi does not define the character of a factor of M"
  j = findfirst(N -> submodule_character(N) == chi, isotypical_factors(M))
  @assert j !== nothing
  return isotypical_factors(M)[j]
end

@doc Markdown.doc"""
    character_decomposition(M::CharacterGrassmannian)
                                     -> Vector{Oscar.GAPGroupClassFunction}

Given the character Grassmannian variety `M` parametrizing `t`-dimensional
submodules of a module `V` with given character, return the isotypical
characters associated to its isotypical factors.
"""
character_decomposition(M::CharacterGrassmannian) = M.dec

@doc Markdown.doc"""
    submodule_character(M::CharacterGrassmannian) -> Oscar.GAPGroupClassFunction

Given the character Grassmannian variety `M` parametrizing submodules
with character `chi` of a module `V`, return `chi`.
"""
submodule_character(M::CharacterGrassmannian) = sum(M.dec)

submodule_dimension(M::CharacterGrassmannian) = Int(degree(sum(M.dec)))

module_representation(M::CharacterGrassmannian) = module_representation(isotypical_factors(M)[1])

is_irreducible(M::CharacterGrassmannian) = true

is_empty(M::CharacterGrassmannian) = false

projective_dimension(M::CharacterGrassmannian) = sum(projective_dimension.(isotypical_factors(M)))

defining_ideal(M::CharacterGrassmannian) = error("Not yet implemented")

function describe(M::CharacterGrassmannian)
  ifs = isotypical_factors(M)
  length(ifs) == 1 && return describe(ifs[1])
  str = ""
  rep = module_representation(M)
  F = base_field(representation_ring(rep))
  chi = submodule_character(M)
  cds = character_decomposition(chi)
  t = cds[1][1]
  n = Int(scalar_product(character_representation(rep), cds[1][2]))
  if n == t
    str *= "{pt}"
  elseif t == 1
    str *= "PP^$(n-1)_F"
  else
    str *= "Gr($t, F^$n)"
  end
  for j in 2:length(ifs)
    str *= " x "
    t = cds[j][1]
    n = Int(scalar_product(character_representation(rep), cds[j][2]))
    if n == 1
      str *= "{pt}"
    elseif t == 1
      str *= "PP^$(n-1)_F"
    else
      str *= "Gr($t, F^$n)"
    end
  end
  if projective_dimension(M) != 0
    str *= " where F = $F"
  end
  return str
end

@doc Markdown.doc"""
    parametrization_data(M::CharacterGrassmannian)
                                        -> Vector{Tuple{Vector{MatSpaceElem}, Int}}

Given a symmetric Grassmannian parametrizing submodules of a module `V` with
given character `chi`, return a list of tuples `(B, n)` where `B` is a basis
of the group `Hom(N, V)` of homomorphisms from `N` to `V`, where `N` is the simple
module of the representation ring `V` affording the simple character associated
to an isotypical component of `chi` and `n` is the size of the submodules of
`Hom(N, V)` parametrized by `M`.

The elements of `B` are given by matrices in the standard coordinates of `V` and `N`.
"""
parametrization_data(M::CharacterGrassmannian) = parametrization_data.(isotypical_factors(M))

@doc Markdown.doc"""
    standard_element(M::CharacterGrassmannian)
                                             -> Vector{Vector{MatSpaceElem}}

Given a symmetric Grassmannian parametrizing submodules of a module `V` with
given character, return a list whose entries are specific elements parametrized
by the isotypical factors of`M`.
"""
standard_element(M::CharacterGrassmannian) = standard_element.(isotypical_factors(M))

###############################################################################

### Determinant Grassmannians

submodule_dimension(M::DeterminantGrassmannian) = M.d

module_representation(M::DeterminantGrassmannian) = M.rep

@doc Markdown.doc"""
    submodule_determinant_character(M::DeterminantGrassmannian)
                                                 -> Oscar.GAPGroupClassFunction

Given the determinant Grassmannian `M` parametrizing `d`-dimensional submodule
of a module `V` with determinant character `chi`, return `chi`.
"""
submodule_determinant_character(M::DeterminantGrassmannian) = M.det_char

@doc Markdown.doc"""
    irreducible_components(M::DeterminantGrassmannian{S, T, U}) where {S, T, U}
                                        -> Vector{CharacterGrassmannian{S, T, U}}

Given the determinant Grassmannian `M` parametrizing `d`-dimensional submodule
of a module `V` with determinantal character `chi`, return the irreducible
components of `M`. These are the character Grassmannians parametrizing submodules
of `V` with a given character `mu` such that `mu` has dimension `d` and determinant `chi`.
"""
irreducible_components(M::DeterminantGrassmannian) = M.irr_comp

is_irreducible(M::DeterminantGrassmannian) = length(irreducible_components(M)) == 1

is_empty(M::DeterminantGrassmannian) = length(irreducible_components(M)) == 0

projective_dimension(M::DeterminantGrassmannian) = maximum(projective_dimension.(M.irr_comp); init = -1)

function defining_ideal(M::DeterminantGrassmannian)
  @req !is_empty(M) "M is empty"
  r = module_representation(M)
  chi = submodule_determinant_character(M)
  t = submodule_dimension(M)
  return _defining_ideal_determinant_grassmannian(r, chi, t)
end

###############################################################################

### Invariant Grassmannians

submodule_dimension(M::InvariantGrassmannian) = M.d

module_representation(M::InvariantGrassmannian) = M.rep

@doc Markdown.doc"""
    irreducible_components(M::InvariantGrassmannian{S, T, U}) where {S, T, U}
                                        -> Vector{CharacterGrassmannian{S, T, U}}

Given the invariant Grassmannian `M` parametrizing `d`-dimensional submodule
of a module `V`, return the irreducible components of `M`. These are the character
Grassmannians parametrizing submodules of `V` with a given character `mu` such that
`mu` has dimension `d`.
"""
irreducible_components(M::InvariantGrassmannian) = M.irr_comp

function irreducible_component(M::InvariantGrassmannian, chi::Oscar.GAPGroupClassFunction)
  irr = irreducible_components(M)
  chis = submodule_character.(irr)
  j = findfirst(j -> chis[j] == chi, 1:length(chis))
  @req j !== nothing "chi is not a character associated to an irreducible component of M"
  return irr[j]
end

is_irreducible(M::InvariantGrassmannian) = length(M.irr_comp) == 1

is_empty(M::InvariantGrassmannian) = length(M.irr_comp) == 0

projective_dimension(M::InvariantGrassmannian) = maximum(projective_dimension.(M.irr_comp); init = -1)

function defining_ideal(M::InvariantGrassmannian)
  @req !is_empty(M) "M is empty"
  r = module_representation(M)
  t = submodule_dimension(M)
  return _defining_ideal_invariant_grassmannian(r, t)
end

###############################################################################
#
# I/O Printing
#
###############################################################################

### Isotypical Grassmannian

function Base.show(io::IO, ::MIME"text/plain", M::IsotypicalGrassmannian)
  chi = submodule_character(M)
  println(io, "Symmetric Grassmannian of $(degree(chi))-dimensional submodules of")
  println(io, module_representation(M))
  println(io, "with isotypical character")
  print(io, chi)
end

function Base.show(io::IO, M::IsotypicalGrassmannian)
  print(io, "Isotypical Grassmannian of dimension $(projective_dimension(M))")
end

###############################################################################

### Character Grassmannian

function Base.show(io::IO, ::MIME"text/plain", M::CharacterGrassmannian)
  chi = submodule_character(M)
  println(io, "Symmetric Grassmannian of $(degree(chi))-dimensional submodules of")
  println(io, module_representation(M))
  println(io, "with character")
  print(io, chi)
end

function Base.show(io::IO, M::CharacterGrassmannian)
  print(io, "Character Grassmannian of dimension $(projective_dimension(M))")
end

###############################################################################

### Determinant Grassmannian

function Base.show(io::IO, ::MIME"text/plain", M::DeterminantGrassmannian)
  chi = submodule_determinant_character(M)
  println(io, "Symmetric Grassmannian of $(M.d)-dimensional submodules of")
  println(io, module_representation(M))
  println(io, "with determinant character")
  print(io, chi)
end

function Base.show(io::IO, M::DeterminantGrassmannian)
  print(io, "Determinant Grassmannian of dimension $(projective_dimension(M))")
end

###############################################################################

### Invariant Grassmannian

function Base.show(io::IO, ::MIME"text/plain", M::InvariantGrassmannian)
  println(io, "Symmetric Grassmannian of $(M.d)-dimensional submodules of")
  print(io, module_representation(M))
end

function Base.show(io::IO, M::InvariantGrassmannian)
  print(io, "Invariant Grassmannian of dimension $(projective_dimension(M))")
end

###############################################################################
#
# Constructors
#
###############################################################################

function _submodules_space_isotypical_as_vs(rep::LinRep{S, T, U},
                                            chi::Oscar.GAPGroupClassFunction) where {S, T, U}
  @assert is_isotypical(chi)
  RR = representation_ring(rep)
  cd = character_decomposition(chi)
  alpha, chis = cd[1]
  B = basis_isotypical_component(rep, chis)
  n = Int(degree(chis))
  d = length(B)
  @assert alpha <= d
  F = base_field(RR)
  V = VectorSpace(F, d)

  function _basis_parametrisation(v)
    return sum([v[i]*B[i] for i in 1:length(B)])
  end

  return MapFromFunc(_basis_parametrisation, V, parent(B[1]))
end

@doc Markdown.doc"""
    isotypical_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction) where {S, T, U}
    isotypical_grassmannian(prep::ProjRep{S, T, U, V}, chi::Oscar.GAPGroupClassFunction) where {S, T, U, V}
                                                                                                          -> IsotypicalGrassmannian{S, T, U}

Given a linear representation `rep` of a finite group `E` on a finite dimensional
`F`-vector space `V`  and an isotypical `F`-character `chi` of `E` where `F` is a 
splitting field of `E` of characteristic zero, return the symmetric Grassmannian
parametrizing isotypical $FE$-submodules of $(V, rep)$ affording the character `chi`.

In the case of we are given a projective representation `prep`, it is constructed
with respect to a stored linear lift `prep`.
"""
function isotypical_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction) where {S, T, U}
  nu = character_representation(rep)
  @req is_isotypical(chi) "chi is not isotypical"
  @req is_constituent(nu, chi) "chi is not a constituent of the character of rep"
  f = _submodules_space_isotypical_as_vs(rep, chi)
  d = Int(scalar_product(chi, nu-chi))
  return IsotypicalGrassmannian(chi, rep, d, f)
end

isotypical_grassmannian(prep::ProjRep, chi::Oscar.GAPGroupClassFunction) = isotypical_grassmannian(linear_lift(prep), chi)

@doc Markdown.doc"""
    character_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction) where {S, T, U}
    character_grassmannian(prep::ProjRep{S, T, U, V}, chi::Oscar.GAPGroupClassFunction) where {S, T, U, V}
                                                                                            -> CharacterGrassmannian{S, T, U}

Given a linear representation `rep` of a finite group `E` on a finite dimensional
`F`-vector space `V`  and a `F`-character `chi` of `E` where `F` is a splitting
field of `E` of characteristic zero, return the symmetric Grassmannian parametrizing
$FE$-submodules of $(V, rep)$ affording the character `chi`.

In the case of we are given a projective representation `prep`, it is constructed
with respect to a stored linear lift `prep`.
"""
function character_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction) where {S, T, U}
  nu = character_representation(rep)
  @req is_constituent(nu, chi) "chi is not a constituent of the character of prep"
  dec = character_decomposition(chi)
  constituents = IsotypicalGrassmannian{S, T, U}[]
  for d in  dec
    M = isotypical_grassmannian(rep, d[1]*d[2])
    push!(constituents, M)
  end

  modi = CharacterGrassmannian(constituents)
  @assert projective_dimension(modi) == scalar_product(chi, nu-chi)
  return modi
end

character_grassmannian(prep::ProjRep, chi::Oscar.GAPGroupClassFunction) = character_grassmannian(linear_lift(prep), chi)

@doc Markdown.doc"""
    invariant_grassmannian(rep::LinRep{S, T, U}, t::Int) where {S, T, U}
    invariant_grassmannian(prep::ProjRep{S, T, U, V}, t::Int) where {S, T, U, V}
                                                                     -> InvariantGrassmannian{S, T, U}

Given a linear representation `rep` of a finite group `E` on a finite dimensional
`F`-vector space `V` where `F` is a splitting field of `E` of characteristic zero,
and an integer `t`, return the symmetric Grassmannian parametrizing `t`-dimensional
$FE$-submodules of $(V, rep)$.

In the case of we are given a projective representation `prep`, it is constructed
with respect to a stored linear lift `prep`.
"""
function invariant_grassmannian(rep::LinRep{S, T, U}, t::Int) where {S, T, U}
  @req 1 <= t < dimension_representation(rep) "t must be positive and (strictly) smaller than the dimension of rep"
  chis = constituents(character_representation(rep), t)
  M = CharacterGrassmannian{S, T, U}[]
  for chi in chis
    N = character_grassmannian(rep, chi)
    push!(M, N)
  end

  return InvariantGrassmannian(rep, M, t)
end

invariant_grassmannian(prep::ProjRep, t::Int) = invariant_grassmannian(linear_lift(prep), t)

@doc Markdown.doc"""
    determinant_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction, t::Int) where {S, T, U}
    determinant_grassmannian(prep::ProjRep{S, T, U, V}, chi::Oscar.GAPGroupClassFunction, t::Int) where {S, T, U, V}
                                                                                                                      -> DeterminantalGrassmannian{S, T, U}

Given a linear representation `rep` of a finite group `E` on a finite dimensional
`F`-vector space `V`, a `F`-character `chi` of `E` where `F` is a splitting field
of `E` of characteristic zero, return the symmetric Grassmannian parametrizing
`t`-dimensional $FE$-submodules of $(V, rep)$ with determinant character `chi`.

In the case of we are given a projective representation `prep`, it is constructed
with respect to a stored linear lift `prep`.
"""
function determinant_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction, t::Int) where {S, T, U}
  @req Int(degree(chi)) == 1 "chi must be a degree 1 character"
  chis = constituents(character_representation(rep), t)
  filter!(nu -> determinant(nu) == chi, chis)
  irr_comp = CharacterGrassmannian{S, T, U}[]
  for nu in chis
    push!(irr_comp, character_grassmannian(rep, nu))
  end
  return DeterminantGrassmannian(rep, chi, irr_comp, t)
end

determinant_grassmannian(prep::ProjRep, chi::Oscar.GAPGroupClassFunction, t::Int) = determinant_grassmannian(linear_lift(prep), chi, t)

###############################################################################
#
# Defining ideals
#
###############################################################################

###  Intersections with Grassmannians

function _intersection_with_grassmannian(V::Vector{T}, n::Int, t::Int; S = nothing) where T
  F = base_ring(V[1])
  
  if S === nothing
    S, _ = grade(PolynomialRing(F, "x" => 0:binomial(n, t)-1)[1])
  end
  
  X = ProjectiveScheme(S)
  if t == 1
    ideal_Gr = ideal(S, [S(0)])
  else
    ideal_Gr = grassmann_pluecker_ideal(S, t, n)
  end
  Grtn = subscheme(X, ideal_Gr)
  B = reduce(vcat, V)
  _, K = right_kernel(B)
  
  if ncols(K) == 0
    return ideal_Gr
  end
  
  ideal_PV = ideal(S, vec(collect(transpose(matrix(gens(S)))*K)))
  PV = subscheme(X, ideal_PV)
  J = modulus(OO(intersect(affine_cone(Grtn), affine_cone(PV))))
  J = ideal(S, [map_coefficients(x -> F(x), p, parent = S) for p in gens(J)])
  J = saturation(J, ideal(S, gens(S)))
  return J
end

### For invariant and determinant grassmannians

function _defining_ideal_determinant_grassmannian(r::LinRep, chi::Oscar.GAPGroupClassFunction, t::Int)
  @req degree(chi) == 1 "chi must be a linear character"
  @req chi in irreducible_characters_underlying_group(representation_ring(r)) "chi is not a character of the undrlying group of r"
  rt = t == 1 ? r : exterior_power_representation(r, t)
  F = base_field(representation_ring(r))
  k = dimension_representation(rt)
  S, _ = grade(PolynomialRing(F, "x" => 0:k-1)[1])
  bas = basis_isotypical_component(rt, chi)
  if length(bas) == 0
    return ideal(S, [S(1)])
  end
  return _intersection_with_grassmannian(bas, dimension_representation(r), t, S = S)::ideal_type(S)
end

function _defining_ideal_invariant_grassmannian(r::LinRep, t::Int)
  rt = t == 1 ? r : exterior_power_representation(r, t)
  F = base_field(representation_ring(r))
  cds = character_decomposition(rt)
  chis = [cd[2] for cd in cds if Int(degree(cd[2])) == 1]
  k = dimension_representation(rt)
  S, _ = grade(PolynomialRing(F, "x" => 0:k-1)[1])
  irre = ideal(S, gens(S))
  I = ideal(S, [S(1)])
  for chi in chis
    bas = basis_isotypical_component(rt, chi)
    J = _intersection_with_grassmannian(bas, dimension_representation(r), t, S=S)
    I = saturation(I*J, irre)
  end
  return I
end

