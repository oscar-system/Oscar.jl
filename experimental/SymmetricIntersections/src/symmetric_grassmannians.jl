###############################################################################
#
# Accessors/Attributes
#
###############################################################################

@doc raw"""
    is_irreducible(M::T) where T <: SymGrass -> Bool

Given a symmetric Grassmannian, return whether it is irreducible.
"""
is_irreducible(M::T) where T <: SymGrass

@doc raw"""
    is_empty(M::T) where T <: SymGrass -> Bool

Given a symmetric Grassmannian, return whether it is empty.
"""
is_empty(M::T) where T <: SymGrass

@doc raw"""
    projective_dimension(M::T) where T <: SymGrass -> Int

Given a symmetric Grassmannian, return its dimension as projective
variety.
"""
projective_dimension(M::T) where T <: SymGrass

@doc raw"""
    submodule_dimension(M::T) where T <: SymGrass -> Int

Given a symmetric Grassmannian, return the dimension of the submodules
it parametrizes.
"""
submodule_dimension(M::T) where T <: SymGrass

@doc raw"""
    module_representation(M::T) where T <: SymGrass -> LinRep

Given a symmetric Grassmannian parametrizing submodules of a module `V`,
return the representation mapping associated to `V`.
"""
module_representation(M::T) where T <: SymGrass

@doc raw"""
    defining_ideal(M::T) where T <: SymGrass -> MPolyIdeal

Given a symmetric Grassmannian parametrizing submodules of a module `V`,
return the ideal defining its structure as projective variety.
"""
defining_ideal(M::T) where T <: SymGrass

@doc raw"""
    describe(M::T) where T <: Union{IsotGrass,
                                    CharGrass}  -> nothing
                           
Given a symmetric Grassmannian parametrizing submodules of a module `V`,
return a description of M in terms of Grassmannian varieties.
"""
describe(M::T) where T <: Union{IsotGrass, CharGrass}

###############################################################################

### Isotypical Grassmannians

@doc raw"""
    submodule_character(M::IsotGrass) -> Oscar.GAPGroupClassFunction

Given the isotypical Grassmannian variety `M` parametrizing isotypical submodules
with character `chi` of a module `V`, return `chi`.
"""
submodule_character(M::IsotGrass) = M.chi

submodule_dimension(M::IsotGrass) = degree(Int, M.chi)

module_representation(M::IsotGrass) = M.rep_mod

is_irreducible(M::IsotGrass) = true

is_empty(M::IsotGrass) = false

projective_dimension(M::IsotGrass) = M.t

function defining_ideal(M::IsotGrass)
  rep = module_representation(M)
  F = base_field(representation_ring(rep))
  chi = submodule_character(M)
  cd = character_decomposition(chi)[1]
  n = scalar_product(Int, character_representation(rep), cd[2])
  t = cd[1]
  if t in Int[1, n-1, n]
    return defining_ideal(projective_space(F, binomial(n, t)-1))
  end
  t = cd[1]
  S, _ = graded_polynomial_ring(F, :x => 0:binomial(n, t)-1)
  return grassmann_pluecker_ideal(S, t, n)
end

function describe(M::IsotGrass)
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

@doc raw"""
    parametrization_data(M::IsotGrass) -> Vector{MatSpaceElem}, Int

Given a symmetric Grassmannian parametrizing submodules of a module `V` with
given isotypical character `chi`, return a tuple `(B, n)` where `B` is a basis
of the group `Hom(N, V)` of homomorphisms from `N` to `V`, where `N` is the simple
module of the representation ring `V` affording the simple character associated
to `chi` and `n` is the size of the submodules of `Hom(N, V)` parametrized by `M`.

The elements of `B` are given by matrices in the standard coordinates of `V` and `N`.
"""
function parametrization_data(M::IsotGrass{S, T, U}) where {S, T, U}
  f = M.vs_struct
  chi = submodule_character(M)
  n, _ = character_decomposition(chi)[1]
  B = f.(gens(domain(f)))::Vector{dense_matrix_type(U)}
  return B, n
end

@doc raw"""
    standard_element(M::IsotGrass) -> Vector{MatSpaceElem}

Given a symmetric Grassmannian parametrizing submodules of a module `V` with
given isotypical character, return a specific elements parametrized by `M`
made from all possible elements in `M`.
"""
function standard_element(M::IsotGrass)
  B, n = parametrization_data(M)
  std_el = eltype(B)[]
  for i in 1:n
    push!(std_el, sum(B[j] for j in i:length(B)))
  end
  return std_el
end

###############################################################################

### Character Grassmannians

@doc raw"""
    isotypical_factors(M::CharGrass{S, T, U})
                    where {S, T, U} -> Vector{IsotGrass{S, T, U}}

Given the character Grassmannian variety `M` parametrizing `t`-dimensional
submodules of a module `V` with given character, return the isotypical
Grassmannian which define `M` as product variety.
"""
isotypical_factors(M::CharGrass) = M.constituents

@doc raw"""
    isotypical_factor(M::CharGrass{S, T, U},
                      chi::Oscar.GAPGroupClassFunction)
                           where {S, T, U} -> IsotGrass{S, T, U}

Given the character Grassmannian variety `M` parametrizing `t`-dimensional
submodules of a module `V` with given character, return the isotypical
Grassmannian in the factors of `M` corresponding to the isotypical character
`chi`.
"""
function isotypical_factor(M::CharGrass{S, T, U},
                              chi::Oscar.GAPGroupClassFunction) where {S, T, U}
  @req is_isotypical(chi) "chi must be an isotypical character"
  @req chi in character_decomposition(M) "chi does not define the character of a factor of M"
  j = findfirst(N -> submodule_character(N) == chi, isotypical_factors(M))
  @assert j !== nothing
  return isotypical_factors(M)[j]
end

@doc raw"""
    character_decomposition(M::CharGrass)
                                     -> Vector{Oscar.GAPGroupClassFunction}

Given the character Grassmannian variety `M` parametrizing `t`-dimensional
submodules of a module `V` with given character, return the isotypical
characters associated to its isotypical factors.
"""
character_decomposition(M::CharGrass) = M.dec

@doc raw"""
    submodule_character(M::CharGrass) -> Oscar.GAPGroupClassFunction

Given the character Grassmannian variety `M` parametrizing submodules
with character `chi` of a module `V`, return `chi`.
"""
submodule_character(M::CharGrass) = sum(M.dec)

submodule_dimension(M::CharGrass) = degree(Int, sum(M.dec))

module_representation(M::CharGrass) = module_representation(isotypical_factors(M)[1])

is_irreducible(M::CharGrass) = true

is_empty(M::CharGrass) = false

projective_dimension(M::CharGrass) = sum(projective_dimension.(isotypical_factors(M)))

defining_ideal(M::CharGrass) = error("Not yet implemented")

function describe(M::CharGrass)
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

@doc raw"""
    parametrization_data(M::CharGrass)
                                        -> Vector{Tuple{Vector{MatSpaceElem}, Int}}

Given a symmetric Grassmannian parametrizing submodules of a module `V` with
given character `chi`, return a list of tuples `(B, n)` where `B` is a basis
of the group `Hom(N, V)` of homomorphisms from `N` to `V`, where `N` is the simple
module of the representation ring `V` affording the simple character associated
to an isotypical component of `chi` and `n` is the size of the submodules of
`Hom(N, V)` parametrized by `M`.

The elements of `B` are given by matrices in the standard coordinates of `V` and `N`.
"""
parametrization_data(M::CharGrass) = parametrization_data.(isotypical_factors(M))

@doc raw"""
    standard_element(M::CharGrass) -> Vector{Vector{MatSpaceElem}}

Given a symmetric Grassmannian parametrizing submodules of a module `V` with
given character, return a list whose entries are specific elements parametrized
by the isotypical factors of`M`.
"""
standard_element(M::CharGrass) = standard_element.(isotypical_factors(M))

###############################################################################

### Determinant Grassmannians

submodule_dimension(M::DetGrass) = M.d

module_representation(M::DetGrass) = M.rep

@doc raw"""
    submodule_determinant_character(M::DetGrass)
                                                 -> GAPGroupClassFunction

Given the determinant Grassmannian `M` parametrizing `d`-dimensional submodule
of a module `V` with determinant character `chi`, return `chi`.
"""
submodule_determinant_character(M::DetGrass) = M.det_char

@doc raw"""
    irreducible_components(M::DetGrass{S, T, U}) where {S, T, U}
                                        -> Vector{CharGrass{S, T, U}}

Given the determinant Grassmannian `M` parametrizing `d`-dimensional submodule
of a module `V` with determinantal character `chi`, return the irreducible
components of `M`. These are the character Grassmannians parametrizing submodules
of `V` with a given character `mu` such that `mu` has dimension `d` and determinant `chi`.
"""
irreducible_components(M::DetGrass) = M.irr_comp

is_irreducible(M::DetGrass) = length(irreducible_components(M)) == 1

is_empty(M::DetGrass) = length(irreducible_components(M)) == 0

projective_dimension(M::DetGrass) = maximum(projective_dimension.(M.irr_comp); init = -1)

function defining_ideal(M::DetGrass)
  @req !is_empty(M) "M is empty"
  r = module_representation(M)
  chi = submodule_determinant_character(M)
  t = submodule_dimension(M)
  return _defining_ideal_determinant_grassmannian(r, chi, t)
end

###############################################################################

### Invariant Grassmannians

submodule_dimension(M::InvGrass) = M.d

module_representation(M::InvGrass) = M.rep

@doc raw"""
    irreducible_components(M::InvGrass{S, T, U}) where {S, T, U}
                                        -> Vector{CharGrass{S, T, U}}

Given the invariant Grassmannian `M` parametrizing `d`-dimensional submodule
of a module `V`, return the irreducible components of `M`. These are the character
Grassmannians parametrizing submodules of `V` with a given character `mu` such that
`mu` has dimension `d`.
"""
irreducible_components(M::InvGrass) = M.irr_comp

function irreducible_component(M::InvGrass, chi::Oscar.GAPGroupClassFunction)
  irr = irreducible_components(M)
  chis = submodule_character.(irr)
  j = findfirst(j -> chis[j] == chi, 1:length(chis))
  @req j !== nothing "chi is not a character associated to an irreducible component of M"
  return irr[j]
end

is_irreducible(M::InvGrass) = length(M.irr_comp) == 1

is_empty(M::InvGrass) = length(M.irr_comp) == 0

projective_dimension(M::InvGrass) = maximum(projective_dimension.(M.irr_comp); init = -1)

function defining_ideal(M::InvGrass)
  @req !is_empty(M) "M is empty"
  r = module_representation(M)
  t = submodule_dimension(M)
  return _defining_ideal_invariant_grassmannian(r, t)
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
  V = vector_space(F, d)

  function _basis_parametrisation(v)
    return sum(v[i]*B[i] for i in 1:length(B))
  end

  return MapFromFunc(V, parent(B[1]), _basis_parametrisation)
end

@doc raw"""
    isotypical_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction) where {S, T, U}
    isotypical_grassmannian(prep::ProjRep{S, T, U, V}, chi::Oscar.GAPGroupClassFunction) where {S, T, U, V}
                                                                                                          -> IsotGrass{S, T, U}

Given a linear representation `rep` of a finite group `E` on a finite dimensional
`F`-vector space `V`  and an isotypical `F`-character `chi` of `E` where `F` is a 
splitting field of `E` of characteristic zero, return the symmetric Grassmannian
parametrizing isotypical $FE$-submodules of $(V, rep)$ affording the character `chi`.

In the case of we are given a projective representation `prep`, it is constructed
with respect to a stored linear lift `rep`.
"""
function isotypical_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction) where {S, T, U}
  nu = character_representation(rep)
  @req is_isotypical(chi) "chi is not isotypical"
  @req is_constituent(nu, chi) "chi is not a constituent of the character of rep"
  f = _submodules_space_isotypical_as_vs(rep, chi)
  d = Int(scalar_product(chi, nu-chi))
  return IsotGrass(chi, rep, d, f)
end

isotypical_grassmannian(prep::ProjRep, chi::Oscar.GAPGroupClassFunction) = isotypical_grassmannian(linear_lift(prep), chi)

@doc raw"""
    character_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction) where {S, T, U}
    character_grassmannian(prep::ProjRep{S, T, U, V}, chi::Oscar.GAPGroupClassFunction) where {S, T, U, V}
                                                                                            -> CharGrass{S, T, U}

Given a linear representation `rep` of a finite group `E` on a finite dimensional
`F`-vector space `V`  and a `F`-character `chi` of `E` where `F` is a splitting
field of `E` of characteristic zero, return the symmetric Grassmannian parametrizing
$FE$-submodules of $(V, rep)$ affording the character `chi`.

In the case of we are given a projective representation `prep`, it is constructed
with respect to a stored linear lift `rep`.
"""
function character_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction) where {S, T, U}
  nu = character_representation(rep)
  @req is_constituent(nu, chi) "chi is not a constituent of the character of prep"
  dec = character_decomposition(chi)
  constituents = IsotGrass{S, T, U}[]
  for d in  dec
    M = isotypical_grassmannian(rep, d[1]*d[2])
    push!(constituents, M)
  end

  modi = CharGrass(constituents)
  @assert projective_dimension(modi) == scalar_product(chi, nu-chi)
  return modi
end

character_grassmannian(prep::ProjRep, chi::Oscar.GAPGroupClassFunction) = character_grassmannian(linear_lift(prep), chi)

@doc raw"""
    invariant_grassmannian(rep::LinRep{S, T, U}, t::Int) where {S, T, U}
    invariant_grassmannian(prep::ProjRep{S, T, U, V}, t::Int) where {S, T, U, V}
                                                                     -> InvGrass{S, T, U}

Given a linear representation `rep` of a finite group `E` on a finite dimensional
`F`-vector space `V` where `F` is a splitting field of `E` of characteristic zero,
and an integer `t`, return the symmetric Grassmannian parametrizing `t`-dimensional
$FE$-submodules of $(V, rep)$.

In the case of we are given a projective representation `prep`, it is constructed
with respect to a stored linear lift `rep`.
"""
function invariant_grassmannian(rep::LinRep{S, T, U}, t::Int) where {S, T, U}
  @req 1 <= t < dimension_representation(rep) "t must be positive and (strictly) smaller than the dimension of rep"
  chis = constituents(character_representation(rep), t)
  M = CharGrass{S, T, U}[]
  for chi in chis
    N = character_grassmannian(rep, chi)
    push!(M, N)
  end

  return InvGrass(rep, M, t)
end

invariant_grassmannian(prep::ProjRep, t::Int) = invariant_grassmannian(linear_lift(prep), t)

@doc raw"""
    determinant_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction, t::Int) where {S, T, U}
    determinant_grassmannian(prep::ProjRep{S, T, U, V}, chi::Oscar.GAPGroupClassFunction, t::Int) where {S, T, U, V}
                                                                                                                      -> DeterminantalGrassmannian{S, T, U}

Given a linear representation `rep` of a finite group `E` on a finite dimensional
`F`-vector space `V`, a `F`-character `chi` of `E` where `F` is a splitting field
of `E` of characteristic zero, return the symmetric Grassmannian parametrizing
`t`-dimensional $FE$-submodules of $(V, rep)$ with determinant character `chi`.

In the case of we are given a projective representation `prep`, it is constructed
with respect to a stored linear lift `rep`.
"""
function determinant_grassmannian(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction, t::Int) where {S, T, U}
  @req Int(degree(chi)) == 1 "chi must be a degree 1 character"
  chis = constituents(character_representation(rep), t)
  filter!(nu -> det(nu) == chi, chis)
  irr_comp = CharGrass{S, T, U}[]
  for nu in chis
    push!(irr_comp, character_grassmannian(rep, nu))
  end
  return DetGrass(rep, chi, irr_comp, t)
end

determinant_grassmannian(prep::ProjRep, chi::Oscar.GAPGroupClassFunction, t::Int) = determinant_grassmannian(linear_lift(prep), chi, t)

###############################################################################
#
# Defining ideals
#
###############################################################################

###  Intersections with Grassmannians

function _intersection_with_grassmannian(V::Vector{T}, n::Int, t::Int;
                                         S::Union{MPolyDecRing, Nothing} = nothing) where T
  F = base_ring(V[1])
  
  if S === nothing
    S, _ = graded_polynomial_ring(F, :x => 0:binomial(n, t)-1)
  end
  
  X = proj(S)
  if t == 1
    ideal_Gr = ideal(S, elem_type(S)[zero(S)])
  else
    ideal_Gr = grassmann_pluecker_ideal(S, t, n)
  end
  Grtn = subscheme(X, ideal_Gr)
  B = reduce(vcat, V)
  K = kernel(B; side = :right)
 
  if ncols(K) == 0
    return ideal_Gr
  end
  
  ideal_PV = ideal(S, vec(collect(matrix(S, 1, nvars(S), gens(S))*K)))
  PV = subscheme(X, ideal_PV)
  _J = modulus(OO(intersect(affine_cone(Grtn)[1], affine_cone(PV)[1])))
  J = ideal(S, elem_type(S)[map_coefficients(F, p; parent = S) for p in gens(_J)])
  J = saturation(J)
  return J::ideal_type(S)
end

### For invariant and determinant grassmannians

function _defining_ideal_determinant_grassmannian(r::LinRep, chi::Oscar.GAPGroupClassFunction, t::Int)
  @req degree(chi) == 1 "chi must be a linear character"
  @req chi in irreducible_characters_underlying_group(representation_ring(r)) "chi is not a character of the undrlying group of r"
  rt = t == 1 ? r : exterior_power_representation(r, t)
  F = base_field(representation_ring(r))
  k = dimension_representation(rt)
  S, _ = graded_polynomial_ring(F, :x => 0:k-1)
  bas = basis_isotypical_component(rt, chi)
  if length(bas) == 0
    return ideal(S, elem_type(S)[one(S)])
  end
  return _intersection_with_grassmannian(bas, dimension_representation(r), t; S)::ideal_type(S)
end

function _defining_ideal_invariant_grassmannian(r::LinRep, t::Int)
  rt = t == 1 ? r : exterior_power_representation(r, t)
  F = base_field(representation_ring(r))
  cds = character_decomposition(rt)
  chis = [cd[2] for cd in cds if Int(degree(cd[2])) == 1]
  k = dimension_representation(rt)
  S, _ = graded_polynomial_ring(F, :x => 0:k-1)
  I = ideal(S, elem_type(S)[S(1)])
  for chi in chis
    bas = basis_isotypical_component(rt, chi)
    J = _intersection_with_grassmannian(bas, dimension_representation(r), t; S)
    I = saturation(I*J)
  end
  return I
end

