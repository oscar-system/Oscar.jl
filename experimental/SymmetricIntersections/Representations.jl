GG = GAP.Globals

import Oscar: is_irreducible, base_field, is_submodule, is_equivalent, is_projective

export action_on_submodule,
       affording_representation,
       all_characters,
       all_irreducible_representations,
       associated_schur_cover,
       basis_exterior_power,
       basis_isotypical_component,
       center_of_character,
       character_decomposition,
       character_linear_lift,
       character_representation,
       character_table_underlying_group,
       complement_submodule,
       constituents,
       determinant,
       dimension_linear_lift,
       dimension_representation,
       direct_sum_representation,
       dual_representation,
       exterior_power_representation,
       faithful_projective_representations,
       generators_underlying_group,
       homogeneous_polynomial_representation,
       irreducible_characters_underlying_group,
       irreducible_affording_representation,
       is_character,
       is_constituent,
       is_faithful,
       is_isotypical,
       is_isotypical_component,
       is_similar,
       is_subrepresentation,
       isotypical_components,
       linear_lift,
       linear_representation,
       matrix_representation,
       matrix_representation_linear_lift,
       projective_representation,
       quotient_representation,
       representation,
       representation_mapping,
       representation_mapping_linear_lift,
       representation_ring,
       representation_ring_linear_lift,
       symmetric_power_representation,
       to_equivalent_block_representation,
       underlying_group

"""
Here are some tools for representation theory and projective representations.
We use representation rings as parent object and we represent projective
representations of a group `G` by a linear lift of a Schur cover of `G`.
Here we use by default as splitting field the cyclotomic field of degree the
exponent of the asscoiated gorup, but a priori functions should work if one
allow representation ring to be constructed over another splitting field.

Note that all this code is well supported for finite groups and splitting
fields of characteristic zero, i.e. `QQAb` or appropriate cyclotomic fields
for instance. Outside of this context, I don't guarantee that everything is
working, or return the expected results. However, one should be able to adapt
part of my codes to extend a bit more the notion of representations in Oscar.
"""

############################################################################
#
#  Accessors
#
############################################################################

### Representation ring

@doc Markdown.doc"""
    base_field(RR::RepRing{S, T}) where {S, T} -> S

Given a representation ring `RR`, return the (cached) splitting field of the
associated finite group.
"""
base_field(RR::RepRing) = RR.field

@doc Markdown.doc"""
    underlying_group(RR::RepRing{S, T}) where {S, T} -> T

Return the underlying group associated to the representation ring `RR`.
"""
underlying_group(RR::RepRing) = RR.group

@doc Markdown.doc"""
    generators_underlying_group(RR::RepRing) -> Vector

Return a (cached) finite set of generators for the underlying group of `RR`.
"""
generators_underlying_group(RR::RepRing) = RR.gens

@doc Markdown.doc"""
    character_table_underlying_group(RR::RepRing)
                                           -> Oscar.GAPGroupCharacterTable

Return the `QQAb`-character table of the underlying group of `RR`, where
`QQAb` is an abelian closure of the rationals.
"""
character_table_underlying_group(RR::RepRing) = RR.ct

@doc Markdown.doc"""
    irreducible_characters_underlying_group(RR::RepRing)
                                     -> Vector{Oscar.GAPGroupClassFunction}

Return the list of irreducible `QQAb`-characters of the underlying group of `RR`,
where `QQAb` is an abelian closure of the rationals.
"""
irreducible_characters_underlying_group(RR::RepRing) = RR.irr

###############################################################################

### Linear representations of finite group

@doc Markdown.doc"""
    representation_ring(LR::LinRep{S, T, U}) where {S, T, U}
                                                   -> RepRing{S, T}

Return the representation ring which `LR` belongs to.
"""
representation_ring(LR::LinRep) = LR.rep_ring

@doc Markdown.doc"""
    representation_mapping(LR::LinRep) -> GAPGroupHomomorphism

Return the underlying representation mapping associated to `LR`.
"""
representation_mapping(LR::LinRep) = LR.f

@doc Markdown.doc"""
    character_representation(LR::LinRep) -> Oscar.GAPGroupClassFunction

Return the corresponding character to the linear representation `LR`.
"""
character_representation(LR::LinRep) = LR.char

@doc Markdown.doc"""
    underlying_group(LR::LinRep{S, T, U}) -> T

Return the underlying group of `LR`.
"""
underlying_group(LR::LinRep) = underlying_group(representation_ring(LR))

@doc Markdown.doc"""
    dimension_representation(LR::LinRep) -> Int

Return the dimension of the underlying module of `LR`.
"""
dimension_representation(LR::LinRep) = Int(degree(LR.char))

###############################################################################

### Projective representations of finite group

@doc Markdown.doc"""
    linear_lift(PR::ProjRep{S, T, U, V} where {S, T, U, V} -> LinRep{S, T, U}

Return a (cached) linear lift of `PR` to a Schur cover of the underlying group
of `PR`.
"""
linear_lift(PR::ProjRep) = PR.LR

@doc Markdown.doc"""
    underlying_group(PR::ProjRep{S, T, U, V}) where {S, T, U, V} -> V

Return the underlying group of `PR`.
"""
underlying_group(PR::ProjRep) = codomain(PR.p)

@doc Markdown.doc"""
    associated_schur_cover(PR::ProjRep) -> GAPGroupHomomorphism

Return the Schur cover along which the projective representation is
lifted.
"""
associated_schur_cover(PR::ProjRep) = PR.p

@doc Markdown.doc"""
    representation_ring_linear_lift(PR::ProjRep{S, T, U, V}) where {S, T, U, V}
                                                                -> RepRing{S, T}

Return the representation ring of the linear lifts of `PR`.
"""
representation_ring_linear_lift(PR::ProjRep) = representation_ring(PR.LR)

@doc Markdown.doc"""
    representation_mapping_linear_lift(PR::ProjRep) -> GAPGroupHomomorphism

Return the representation map associated to a linear lift of `PR`.
"""
representation_mapping_linear_lift(PR::ProjRep) = representation_mapping(linear_lift(PR))

@doc Markdown.doc"""
    character_linear_lift(PR::ProjRep) -> Oscar.GAPGroupClassFunction

Return the character of the linear lift of `PR` which is cached.
"""
character_linear_lift(PR::ProjRep) = character_representation(PR.LR)

@doc Markdown.doc"""
    dimension_linear_lift(PR::ProjRep) -> Int

Return the dimension of any linear lift of `PR`.
"""
dimension_linear_lift(PR::ProjRep) = dimension_representation(PR.LR)

###############################################################################
#
#  Character operations (not yet available in src)
#
###############################################################################

@doc Markdown.doc"""
    character_representation(RR::RepRing, f::GAPGroupHomomorphism)
                                                 -> Oscar.GAPGroupClassFunction

Given a representation mapping `f` of the underlying group of `RR`, return the
character afforded by `f`.
"""
function character_representation(RR::RepRing{S, T}, f::GAPGroupHomomorphism{T, MatrixGroup{U, AbstractAlgebra.Generic.MatSpaceElem{U}}}) where {S, T, U}
  @req domain(f) === underlying_group(RR) "f does not define a representation in the given representation ring"
  @req codomain(f) isa MatrixGroup "f should take value in a matrix group"
  Q = abelian_closure(QQ)[1]
  coll = elem_type(Q)[]
  H = conjugacy_classes(underlying_group(RR))
  # Any representation rho affords the character e -> trace(rho(e)) 
  for h in H
    b = f(representative(h))
    push!(coll, Q(trace(b)))
  end
  c = Oscar.group_class_function(character_table_underlying_group(RR), coll)
  return c
end

@doc Markdown.doc"""
    character_decompositon(rep::LinRep)
                           -> Vector{Tuple{Int, Oscar.GAPGroupClassFunction}}

Return the decomposition of the character afforded by `rep`.
"""
character_decomposition(rep::LinRep) = character_decomposition(character_representation(rep))

@doc Markdown.doc"""
    character_decomposition(prep::ProjRep)
                          -> Vector{Tuple{Int, Oscar.GAPGroupClassFunction}}

Return the decomposition of the character afforded by the cached linear lift
of `prep`.
"""
character_decomposition(prep::ProjRep) = character_decomposition(linear_lift(prep))

# We decompose characters as sum of positive multiple of some irreducible characters.
# Actually, we could decompose any class functions in this way, one will obtain some
# negative multiplicities then.
@doc Markdown.doc"""
    character_decomposition(char::Oscar.GAPGroupClassFunction)
                              -> Vector{Tuple{Int, Oscar.GAPGroupClassFunction}}

Given a character `char`, return a set of tuple `(alpha, chi)` where `chi` is
irreducible in the character table of `char` and $alpha = scalar_product(char, chi)$
is the multiplicity of chi in char. If `alpha = 0`, it is not added.
"""
function character_decomposition(char::Oscar.GAPGroupClassFunction)
  ct = char.table
  decomp = Tuple{Int, typeof(char)}[]
  for chi in ct
    alpha = Int(scalar_product(chi, char))
    alpha == 0 && continue
    push!(decomp, (alpha, chi))
  end
  sort!(decomp, by = c -> Int(degree(c[2])))
  return decomp
end

# This is just an import of a GAP function
@doc Markdown.doc"""
    is_character(chi::Oscar.GAPGroupClassFunction) -> Bool

Given a class function `chi` on a group `G`, return whether `chi` defines a
character of `G` (over its codomain).
"""
is_character(chi::Oscar.GAPGroupClassFunction) = GG.IsCharacter(chi.table.GAPTable, chi.values)::Bool

@doc Markdown.doc"""
    is_constituent(chi::T, nu::T) where T <: Oscar.GAPGroupClassFunction -> Bool

Given two characters `chi` and `nu` of a group `G`, return whether `nu` is
constituent of `chi`.
"""
is_constituent(chi::T, nu::T) where T <: Oscar.GAPGroupClassFunction = chi == nu || is_character(chi-nu)

@doc Markdown.doc"""
    is_isotypical(char::Oscar.GAPGroupClassFunction) -> Bool

Given a character `char` of a group `G`, return whether it is isotypical.

A character is said to be isotypical if it is a positive multiple of an
irreducible character.
"""
is_isotypical(char::Oscar.GAPGroupClassFunction) = length(character_decomposition(char)) == 1

# We compute all consistuents of chi of degree t
@doc Markdown.doc"""
    constituents(chi::Oscar.GAPGroupClassFunction, t::Int)
                                     -> Vector{Oscar.GAPGroupClassFunction)

Given a character `chi` and an integer `t`, return all the consistuents of
`chi` of degree `t`.
"""
function constituents(chi::Oscar.GAPGroupClassFunction, t::Int)
  @req is_character(chi) "chi should be a character"
  @req t >= 0 "t must be non negative"
  
  # sort of fallback when t = degree(chi)
  if t == 0
    return [zero(chi)]
  end
  
  # constituents of degree t in chi corresponds bijectively to their
  # complement which are of degree degree(chi)-t. So to do the smallest
  # computations, we always reduce to the case t <= floor(degree(chi)/2)
  if t > div(Int(degree(chi)), 2)
    els = constituents(chi, Int(degree(chi)) - t)
    return [chi-nu for nu in els]
  end
  
  cd = character_decomposition(chi)
  L = [c[2] for c in cd]
  ubs = [c[1] for c in cd]
  
  function _deg(x::Oscar.GAPGroupClassFunction)
    return ZZ(degree(x))
  end
  
  # _deg(L) represents the degree of all irreducible constituents of chi,
  # and ubs keeps track of their multiplicities in chi. So that we return
  # all the t-elevations of (L, _deg) where we impose some upper bound
  # on their of time we can take each index (otherwise the resulting
  # character won't be a constituent of chi
  el = elevator(L, _deg, t, ubs = ubs)
  (number_of_elevations(el) == 0) && return Oscar.GAPGroupClassFunction[]
  return Oscar.GAPGroupClassFunction[sum(L[l]) for l in el]
end

@doc Markdown.doc"""
    all_characters(RR::RepRing, t::Int) -> Vector{Oscar.GAPGroupClassFunction}

Given a representation ring `RR` associated to a finite group `E` and an integer
`t`, return all characters of `E` of degree `t`.
"""
function all_characters(RR::RepRing, t::Int)
  chis = irreducible_characters_underlying_group(RR)
  el = elevator(chis, x -> ZZ(degree(x)), t)
  return Oscar.GAPGroupClassFunction[sum(chis[l]) for l in el]
end

@doc Markdown.doc"""
    determinant(chi::Oscar.GAPGroupClassFunction) -> Oscar.GAPGroupClassFunction

Given a character `chi` of degree `n`, return its determinant character,
which corresponds to the `n`-th antisymmetric parts of `chi`.
"""
determinant(chi::Oscar.GAPGroupClassFunction) = exterior_power(chi, Int(degree(chi)))

@doc Markdown.doc"""
    center_of_character(chi::Oscar.GAPGroupClassFunction) -> GAPGroup, GAPGroupHomomorphism

Return the center of the character `chi` which corresponds to the subgroup
of the underlying group `E` of `chi` consisting on the elements `e` of `E`
such that $chi(e)$ is a primitive root of unity times the degree of `chi`.
"""
function center_of_character(chi::Oscar.GAPGroupClassFunction)
  E = chi.table.GAPGroup
  _H = GG.CenterOfCharacter(chi.values)::GAP.GAP_jll.GapObj
  H = typeof(E)(_H)
  ok, j = is_subgroup(E, H)
  @assert ok
  return H, j
end

###########################################################################
#
# I/O printing
#
###########################################################################

## RR

function Base.show(io::IO, ::MIME"text/plain", RR::RepRing)
  println(io, "Representation ring of")
  println(io, "$(underlying_group(RR))")
  println(io, "over")
  print(io, "$(base_field(RR))")
end

function Base.show(io::IO, RR::RepRing)
  print(io, "Representation ring of finite group over a field of characteristic 0")
end

## LR

function Base.show(io::IO, ::MIME"text/plain", LR::LinRep)
  println(io, "Linear $(dimension_representation(LR))-dimensional representation of")
  println(io, "$(underlying_group(LR))")
  println(io, "over")
  print(io, "$(base_field(representation_ring(LR)))")
end 

function Base.show(io::IO, LR::LinRep)
  print(io, "Linear representation of finite group of dimension $(dimension_representation(LR))")
end

## PR

function Base.show(io::IO, ::MIME"text/plain", PR::ProjRep)
  println(io, "Linear lift of a $(dimension_linear_lift(PR))-dimensional projective representation of")
  println(io, "$(underlying_group(PR))")
  println(io, "over")
  print(io, "$(base_field(representation_ring_linear_lift(PR)))")
end

function Base.show(io::IO, PR::ProjRep)
    print(io, "Linear lift of a projective representation of finite group of dimension $(dimension_linear_lift(PR))")
end

##############################################################################
#
# Constructors
#
##############################################################################

@doc Markdown.doc"""
    representation_ring(E::T) where T <: Oscar.GAPGroup -> RepRing

Given a finite group `E`, return the representation ring over `E` over the
`e`-th cyclotomic field, where `e` is the exponent of `E`.
"""
function representation_ring(E::T) where T <: Oscar.GAPGroup
  e = exponent(E)
  F, _ = CyclotomicField(Int(e), cached = false)
  return RepRing{typeof(F), T}(F, E)
end

# Old functions that I used before turning everything into GAPGroupHomomorphism.
# Still useful if one knows what one does, so that we avoid to always reconstruct
# independently in each functions the homomorphism: we just feed this/these
# internal(s) with an appropriate list of matrices.
"""
    _linear_representation(RR::RepRing{S, T}, mr::Vector{V},
                          chi::Oscar.GAPGroupClassFunction)
                                where {S, T, U, V <: MatElem{U}} -> LinRep{S, T, U}

Return the linear representation in `RR` associated to the matrix representation
`mr`. `chi` must be correspond to the character afforded by this representation.
If not provided, `chi` is automatically computed.

Note: here matrices in `mr` should correspond to the cached set of
generators of the underlying group of `RR` (order matters.)
"""
function _linear_representation(RR::RepRing{S, T}, mr::Vector{V}, chi::Oscar.GAPGroupClassFunction) where {S, T, U, V <: MatElem{U}}
  @req chi.table === character_table_underlying_group(RR) "Character should belong to the character table attached to the given representation ring"
  G = underlying_group(RR)
  mg = MatrixGroup(nrows(mr[1]), base_ring(mr[1]), mr)
  f = hom(G, mg, generators_underlying_group(RR), gens(mg), check = false)
  return linear_representation(RR, f, chi)
end

function _linear_representation(RR::RepRing{S ,T}, mr::Vector{V}) where {S ,T ,U, V <: MatElem{U}}
  G = underlying_group(RR)
  mg = MatrixGroup(nrows(mr[1]), base_ring(mr[1]), mr)
  f = hom(G, mg, generators_underlying_group(RR), gens(mg), check = false)
  return linear_representation(RR, f)
end

# The associated constructors to the corresponding types.
@doc Markdown.doc"""
    linear_representation(RR::RepRing, f::GAPGroupHomomorphism,
                                       chi::Oscar.GAPGroupClassFunction) -> LinRep

Return the linear representation in `RR` associated to the mapping `f`, if the
domain of `f` corresponds to the underlying group of `RR`.
`chi` must be correspond to the character afforded by this representation.
If not provided, `chi` is automatically computed.
"""
function linear_representation(RR::RepRing{S, T}, f::GAPGroupHomomorphism, chi::Oscar.GAPGroupClassFunction) where {S, T}
  @req chi.table === character_table_underlying_group(RR) "Character should belong to the character table attached to the given representation ring"
  @req domain(f) === underlying_group(RR) "f should be defined over the underlying group of the given representation ring"  
  @req codomain(f) isa MatrixGroup "f should take value in a matrix group"
  return LinRep{S, T, elem_type(S)}(RR, f, chi)
end

function linear_representation(RR::RepRing{S ,T}, f::GAPGroupHomomorphism) where {S ,T}
  @req domain(f) === underlying_group(RR) "f should be defined over the underlying group of the given representation ring"
  @req codomain(f) isa MatrixGroup "f should take value in a matrix group"
  chi = character_representation(RR, f)
  @assert chi.table === character_table_underlying_group(RR)
  return LinRep{S, T, elem_type(S)}(RR, f, chi)
end

@doc Markdown.doc"""
    is_projective(rep::LinRep, p::GAPGroupHomomorphism) -> Bool
    is_projective(chi::GAPGroupClassFunction, p::GAPGroupHomomorphism) -> Bool

Given a linear representation `rep` of the domain `E` of the cover `p` affording the
character `chi`, return whether `rep` is `p`-projective, i.e. `rep` reduces to a
projective representation of the codomain of `p`.

This is equivalent to ask that the center of `chi` contains the kernel of `p`.
"""
function is_projective(chi::Oscar.GAPGroupClassFunction, p::GAPGroupHomomorphism)
  @req chi.table.GAPGroup === domain(p) "Incompatible representation ring of rep and domain of the cover p"
  return is_subgroup(center_of_character(chi)[1], kernel(p)[1])[1]
end

is_projective(rep::LinRep, p::GAPGroupHomomorphism) = is_projective(character_representation(rep), p)

@doc Markdown.doc"""
    projective_representation(RR::RepRing, f::GAPGroupHomomorphism,
                                           p::GAPGroupHomomorphism) -> ProjRep

Return a linear representation in `RR` and associating mapping `f` representation
a lift along the cover `p` of a projective representation of $codomain(p)$.

If `check=true`, `f` is check to be `p`-projective.
"""
function projective_representation(RR::RepRing{S, T}, f::GAPGroupHomomorphism{T, MatrixGroup{U, AbstractAlgebra.Generic.MatSpaceElem{U}}}, p::V; check::Bool = true) where {S, T, U, V}
  lr = linear_representation(RR, f)
  if check
    @req underlying_group(RR) === domain(p) "Incompatible representation ring and cover"
    @req is_projective(lr, p) "f is not p-projective"
  end
  return ProjRep{S, T, U, V}(lr, p)
end

@doc Markdown.doc"""
    is_irreducible(rep::LinRep) -> Bool

Return whether `rep` is an irreducible linear representation.
"""
function is_irreducible(rep::LinRep)
  return Oscar.is_irreducible(character_representation(rep))
end

@doc Markdown.doc"""
    is_irreducible(prep::ProjRep) -> Bool

Return whether the underlying linear lift of `prep` is an irreducible
linear representation.
"""
function is_irreducible(prep::ProjRep)
  return is_irreducible(linear_lift(prep))
end

@doc Markdown.doc"""
    is_isotypical(rep::LinRep) -> Bool

Return whether `rep` is a positive multiple of an irreducible linear
representation.
"""
function is_isotypical(rep::LinRep)
  return is_isotypical(character_representation(rep))
end

@doc Markdown.doc"""
    is_isotypical(prep::ProjRep) -> Bool

Return whether the underlying linear lift of `prep` is a positive
multiple of an irreducible linear representation.
"""
function is_isotypical(prep::ProjRep)
  return is_isotypical(linear_lift(prep))
end

# This is an import of the `repsn` GAP package: while they do construct the
# homomorphism without context, here we insist on the fact that the representation
# "belongs to" the given representation ring, whenever it makes sense. Then,
# since we have a fixed underlying group and an parent object for the entries of
# the matrices in the image, both stored in `RR`, we translate the homomorphism
# constructed in GAP to match properly the data here. We also add an extra caching
# in `RR` to avoid having to reconstruct irreducible representations already computed.
@doc Markdown.doc"""
    irreducible_affording_representation(RR::RepRing,
                                         chi::Oscar.GAPGroupClassFunction)
                                                                      -> LinRep

Given a representation ring `RR` and an irreducible character `chi` of the
underlying group of `RR`, return a linear representation in `RR` affording
`chi`.

Note: we store already computed irreducible representations in `RR` to avoid
redundacy.
"""
function irreducible_affording_representation(RR::RepRing{S, T}, chi::Oscar.GAPGroupClassFunction) where {S, T}
  @req chi in irreducible_characters_underlying_group(RR) "chi is not an irreducible character of RR"
  # we have first to inspect the caching: if one already computed
  # an affording representation, we are done
  if has_attribute(RR, :irr_rep)
    irr_rep = get_attribute(RR, :irr_rep)::Dict{Oscar.GAPGroupClassFunction, LinRep{S, T, elem_type(S)}}
  else
    irr_rep = Dict{Oscar.GAPGroupClassFunction, LinRep{S, T, elem_type(S)}}()
  end
  haskey(irr_rep, chi) && return irr_rep[chi] 
  F = base_field(RR)

  # we use F and H to make the translation between the representation
  # mapping on GAP and the one we aim to create on Oscar
  H = generators_underlying_group(RR)

  # the GAP function which we rely on
  rep = GG.IrreducibleAffordingRepresentation(chi.values)::GAP.GAP_jll.GapObj

  # we compute certain images than we will convert then.
  # we use them then to construct the mapping in 
  # `_linear_representation`
  _Mat = GAP.GAP_jll.GapObj[GG.Image(rep, h.X) for h in H]
  Mat = AbstractAlgebra.Generic.MatSpaceElem{elem_type(F)}[]
  for _M in _Mat
    rM = eltype(Mat)[transpose(matrix(F.(m))) for m in _M]
    M = reduce(vcat, rM)::eltype(Mat)
    push!(Mat, M)
  end
  lin = _linear_representation(RR, Mat, chi)

  # and we cache the representation for later
  irr_rep[chi] = lin
  set_attribute!(RR, :irr_rep, irr_rep)

  return lin::elem_type(RR)
end

@doc Markdown.doc"""
    affording_representation(RR::RepRing, chi::Oscar.GAPGroupClassFunction)
                                                                      -> LinRep

Given a representation ring `RR` and a character `chi` of the underlying
group of `RR`, return a linear representation in `RR` affording `chi`.
"""
function affording_representation(RR::RepRing{S, T}, chi::Oscar.GAPGroupClassFunction) where {S, T}
  chi in irreducible_characters_underlying_group(RR) && return irreducible_affording_representation(RR, chi)
  cd = character_decomposition(chi)

  # The blocks for the block diagonal representation. The
  # easiest way for us is to start by looking at the 
  # irreducible constituent of chi and their multiplicities
  # and then we take the direct sum of the corresponding
  # irreducible affording representations with the given multiplicities
  blocks = [c[2] for c in cd for i in 1:c[1]]
  rep = irreducible_affording_representation(RR, blocks[1])
  for j in 2:length(blocks)
    rep = direct_sum_representation(rep, irreducible_affording_representation(RR, blocks[j]))
  end

  return rep::elem_type(RR)
end

# This can be expensive if not necessary. Though, to setup a complex algorithm,
# it could allow to save some time by already precomputing all irreducible
# representations in `RR` and storing them.
@doc Markdown.doc"""
    all_irreducible_representations(RR::RepRing) -> Vector{LinRep}
    all_irreducible_representation(E::GAPGroup) -> Vector{LinRep}

Given a representation ring `RR` with underlying group `E` over a
field `F`, or more generally a finite group `E` with `F = QQAb`, return
affording linear `F`-representations of `E` for all irreducible `F`-characters
of `E`.
"""
function all_irreducible_representations(RR::RepRing{S, T}) where {S, T}
  chis = irreducible_characters_underlying_group(RR)
  return LinRep{S, T, elem_type(S)}[irreducible_affording_representation(RR, chi) for chi in chis]
end

function all_irreducible_representations(E::T) where T <: Oscar.GAPGroup
  RR = representation_ring(E)
  coll = all_irreducible_representations(RR)
  return RR, coll
end

###############################################################################
#
# Operations on representation
#
###############################################################################

### Matrix representation

# In most of our "induced action function", we reconstruct the representations
# from a list of matrices coming from the images of H by the original representation
# where H is the set of generators cached in the corresponding representation ring
@doc Markdown.doc"""
    matrix_representation(LR::LinRep{S, T, U}) -> Vector{MatElem{U}}

Return the matrix representation of the corresponding images of the
cached list of generators of the underlying group of `LR`.
"""
function matrix_representation(LR::LinRep{S, T, U}) where {S, T, U}
  f = representation_mapping(LR)
  return matrix.(gens(codomain(f)))::Vector{dense_matrix_type(U)}
end

@doc Markdown.doc"""
    matrix_representation_linear_lift(PR::ProjRep{S, T, U, V})
                                       where {S, T, U, V} -> Vector{MatElem{U}}

Return representatives (modulo scalars) for the matrix representation of the
underlying actions encoded by `PR`.
"""
matrix_representation_linear_lift(PR::ProjRep{S, T, U, V}) where {S, T, U, V} = matrix_representation(linear_lift(PR))

###############################################################################

### Faithfulness

@doc Markdown.doc"""
    is_faithful(rep::LinRep) -> Bool

Return whether the underlying representation mapping of `rep` is injective.
"""
function is_faithful(rep::LinRep)
  f = representation_mapping(rep)
  return is_injective(f)
end

@doc Markdown.doc"""
    is_projectively_faithful(rep::LinRep, p::GAPGroupHomomorphism) -> Bool
    is_projectively_faithful(chi::Oscar.GAPGroupClassFunction, p::GAPGroupHomomorphism) -> Bool

Given a linear representation `rep` of the domain `E` of the cover `p` affording the
character `chi`, return whether `rep` is `p`-faithful, i.e. `rep` reduces to a faithful
projective representation of the codomain of `p`.

This is equivalent to ask that the center of `chi` coincides with the kernel of `p`.
"""
function is_faithful(chi::Oscar.GAPGroupClassFunction, p::GAPGroupHomomorphism{T, V}) where {T, V}
  E = chi.table.GAPGroup::T
  @req E === domain(p) "Incompatible underlying group of chi and domain of the cover p"
  @req is_projective(chi, p) "chi is not afforded by a p-projective representation"
  Z = center_of_character(chi)[1]::T
  Q = kernel(p)[1]
  return Q.X == Z.X
end

is_faithful(rep::LinRep, p::GAPGroupHomomorphism) = is_faithful(character_representation(rep), p)

@doc Markdown.doc"""
    is_faithful(prep::ProjRep) -> Bool

Given a projective representation `prep` lifting to a linear representation `rep` along
a cover `p`, return whether `rep` is `p`-faithful.
"""
is_faithful(pr::ProjRep) = is_faithful(character_linear_lift(pr), associated_schur_cover(pr))

###############################################################################

### Comparison

@doc Markdown.doc"""
    is_equivalent(rep1::T, rep2::T) where T <: LinRep -> Bool

Given two representations `rep1` and `rep2` of the same representation
ring `RR`, return whether `rep1` and `rep2` are equivalent.

Two linear representations are equivalent if and only if the afford the
same character.
"""
function is_equivalent(rep1::T, rep2::T) where T <: LinRep
  @req representation_ring(rep1) === representation_ring(rep2) "Representations must be in the same representation ring"
  return character_representation(rep1) == character_representation(rep2)
end

@doc Markdown.doc"""
    is_similar(prep1::T, prep2::T) where T <: ProjRep -> Bool

Given two projective representations `prep1` and `prep2` of the same group
with linear lifts in the same representation ring, return whether they are
similar.

Two projective representations are similar if and only if their
respective linear lifts are equivalent modulo a 1-dimensional character.
"""
function is_similar(prep1::T, prep2::T) where T <: ProjRep
  @req representation_ring_linear_lift(prep1) === representation_ring_linear_lift(prep2) "Linear lifts must be in the same representation ring"
  @req underlying_group(prep1) == underlying_group(prep2) "Underlying groups mismatch"
  RR = representation_ring_linear_lift(prep1)
  chi1 = character_linear_lift(prep1)
  chi2 = character_linear_lift(prep2)
  Irr = irreducible_characters_underlying_group(RR)
  return any(irr -> irr*chi1 == chi2, Irr)
end

# This works nicely here because we only have conditions on the action of
# the group on a hypothetical sub-module, without involving any basis. So if the
# comparison of the character return true, we can find a submodule of the module
# defined by `rep` on which the group acts via `rep2`.
@doc Markdown.doc"""
    is_subrepresentation(rep::T, rep2::T) where T <: LinRep -> Bool
    is_subrepresentation(prep::T, prep2::T) where T <: ProjRep -> Bool

Given two linear representations `rep` and `rep2` into the same representation
ring `RR` of the finite group `E` over the field `F`, return
whether the abstract `FE`-module represented by `rep2` can be embedded,
`E`-equivariantly, into the `FE`-module represented by `rep`.

In the case of projective representations, this comparison is done at the level
of linear lifts.
"""
function is_subrepresentation(rep::T, rep2::T) where T <: LinRep
  @req representation_ring(rep) === representation_ring(rep2) "Representation must be in the same representation ring"
  return is_constituent(character_representation(rep), character_representation(rep2))
end

function is_subrepresentation(prep::T, prep2::T) where T <: ProjRep
  return is_subrepresentation(linear_lift(prep), linear_lift(prep2))
end

@doc Markdown.doc"""
    is_isotypical_component(rep::T, rep2::T) where T <: LinRep -> Bool
    is_isotypical_component(prep::T, prep2::T) where T <: ProjRep -> Bool
    
Given two linear representations `rep` and `rep2` into the same representation
ring `RR` of the finite group `E` over the field `F`, return whether
the abstract `FE`-module represented by `rep2` is an isotypical component of
the `FE`-module represented by `rep`.

In the case of projective representations, this comparison is done at the level
of linear lifts.
"""
function is_isotypical_component(rep::T, rep2::T) where T <: LinRep
  @req is_isotypical(rep2) && is_subrepresentation(rep, rep2) "rep2 is not an isotypical subrepresentation of rep"
  chi = character_representation(rep)
  chi2 = character_representation(rep2)
  irr = character_decomposition(chi2)[1][2]
  return scalar_product(irr, chi-chi2) == 0
end

function is_isotypical_component(prep::T, prep2::T) where T <: ProjRep
  return is_isotypical_component(linear_lift(prep), linear_lift(prep2))
end

###############################################################################

### Dual representation

@doc Markdown.doc"""
    dual_representation(rep::T) where T <: LinRep -> T
    dual_representation(prep::T) where T <: ProjRep -> T

Given a representation of a finite group `E` on a finite dimensional `F`-vector
space `V`, where `F` is a splitting field of `E`, return the induced action of
`E` on the dual space of `V`.

In the case of a projective representation `prep`, return the reduction of the
action on the dual of `V` of any linear lift of `prep`.
"""
function dual_representation(rep::LinRep{S, T, U}) where {S, T, U}
  RR = representation_ring(rep)
  mr = matrix_representation(rep)::Vector{dense_matrix_type(U)}
  chi = character_representation(rep)
  mr_dual = dense_matrix_type(U)[transpose(inv(m)) for m in mr]
  chi_dual = conj(chi)
  return _linear_representation(RR, mr_dual, chi_dual)
end

function dual_representation(prep::ProjRep{S, T, U, V}) where {S, T, U, V}
    rep = dual_representation(linear_lift(prep))::LinRep{S, T, U}
    return ProjRep{S, T, U, V}(rep, associated_schur_cover(prep))
end

###############################################################################

### Direct sum in representation ring

@doc Markdown.doc"""
    direct_sum_representation(rep1::T, rep2::T) where T <: LinRep -> T
    direct_sum_representation(prep1::T, prep2::T) where T <: ProjRep -> T

Given two linear representation `rep1` and `rep2` of a finite group `E` on
finite dimensional a `F`-vector spaces `V1` and `V2`, where `F` is a splitting
field of `E`, return the corresponding linear representation on $V1 \oplus V2$.

In the case of projective representations, return the reduction of the
corresponding direct sum representation of the linear lifts.
"""
function direct_sum_representation(rep1::W, rep2::W) where {S, T, U, W <: LinRep{S, T, U}}
  @req representation_ring(rep1) === representation_ring(rep2) "Representations must be in the same representation ring"
  mr1 = matrix_representation(rep1)::Vector{dense_matrix_type(U)}
  mr2 = matrix_representation(rep2)::Vector{dense_matrix_type(U)}
  mr = eltype(mr1)[block_diagonal_matrix([mr1[i], mr2[i]]) for i in 1:length(mr1)]
  chi = character_representation(rep1) + character_representation(rep2)
  return _linear_representation(representation_ring(rep1), mr, chi)
end

function direct_sum_representation(prep1::W, prep2::W) where {S, T, U, V, W <: ProjRep{S, T, U, V}}
  @req representation_ring_linear_lift(prep1) === representation_ring_linear_lift(prep2) "Linear lifts must be in the same representation ring"
  @req associated_schur_cover(prep1) === associated_schur_cover(prep2) "Associated covers mismatch"
  lin = direct_sum_representation(linear_lift(prep1), linear_lift(prep2))
  return ProjRep{S, T, U, V}(lin, associated_schur_cover(prep1))
end

Base.:(+)(rep1::T, rep2::T) where T <: Union{LinRep, ProjRep} = direct_sum_representation(rep1, rep2)

###############################################################################

### Some bases operations

"""
Standars basis of `F^n`
"""
function _standard_basis(F::U, n::Int) where U
  v = zero_matrix(F, n, 1)
  std_bas = typeof(v)[]
  for j=1:n
    vv =deepcopy(v)
    vv[j,1] = one(F)
    push!(std_bas, vv)
  end
  return std_bas
end


@doc Markdown.doc"""
    basis_exterior_power(B::Vector{T}, t::Int) where T -> Vecotr{Vector{T}}

Given a basis `B` of a finite dimensional vector space `V`, return a basis for
$\bigwedge^t V$, ordered accordingly to the usual Plucker coordinates in the
elements of `B`.
"""
function basis_exterior_power(B::Vector{T}, t::Int64) where T
  l = length(B)
  L = [1 for i in 1:l]
  el = elevator(L, t, ubs=L)
  return Vector{T}[B[lis] for lis in el]
end

# I am not sure whether we have somewhere a fast algorithm to compute
# exterior powers of matrices, so I made it up myself. For this I need
# some operations on "tensors", symmetric or alternating, seen as lists
# of lists of the components of their terms.
"""
A bunch of operations on symmetric/anti-symmetric tensors.
"""
function _change_basis(w, basis::Vector{Vector{T}}) where T 
  @req length(w) == length(basis) "Change of basis impossible" 
  gene = [(w[i],basis[i]) for i =1:length(w) if w[i] != parent(w[i])(0)]
  return gene 
end

function _decompose_in_standard_basis(v)
  std_bas = _standard_basis(base_ring(v), nrows(v))
  return std_bas, [[v[i], std_bas[i]] for i=1:nrows(v)]
end

function _same_support(v::Vector{T}, w::Vector{T}) where T
  @req length(v) == length(w) "Tensor must have the same number of compoenents"
  return MSet(v).dict == MSet(w).dict
end

function _div(v::Vector{T}, w::Vector{T}; symmetric = false) where T
  @req _same_support(v,w) "Tensor must be the same support"
  if symmetric
    return 1
  else
    return sign(perm([findfirst(vv -> vv == ww, v) for ww in w]))
  end
end

function _in_basis(v::Vector{T}, basis::Vector{Vector{T}}; symmetric::Bool = false) where T
  @assert length(v) == length(basis[1])
  i = findfirst(w -> _same_support(v,w),basis)
  i === nothing ? (return false, 0) : (return basis[i], _div(v,basis[i], symmetric = symmetric))
end

function _tensor_in_basis(w, bas)
  @assert length(w[1][2]) == length(bas[1])
  v = zeros(parent(w[1][1]), 1, length(bas))
  for ww in w
    j = indexin([ww[2]], bas)[1]
    v[j] = ww[1]
  end
  return matrix(v)
end

###############################################################################

### Symmetric power representation

# here we could have used symmetric tensors, but we can actually use homogeneous
# polynomials which turns out to be faster.
"""
If `mg` represents the action of a group `E` on a finite dimensional vector
space `V`, return an induced action on Sym^d V.
"""
function _action_symmetric_power(mr::Vector{AbstractAlgebra.Generic.MatSpaceElem{S}}, d::Int) where S
  n = ncols(mr[1])
  F = base_ring(mr[1])
  R, _ = grade(PolynomialRing(F, "x"=>1:n, cached=false)[1])
  R1, R1toR = homogeneous_component(R, 1)
  Rd, RdtoR = homogeneous_component(R, d)
  bsp = reverse(RdtoR.(gens(Rd)))
  coll = eltype(mr)[]
  for m in mr
    poly_m = [R1toR(R1(reverse(vec(collect(m[i,:]))))) for i in 1:nrows(m)]
    @assert length(poly_m) == n
    m2 = zero_matrix(F, length(bsp), length(bsp))
    for i = 1:length(bsp)
      b = bsp[i]
      w = evaluate(b, poly_m)
      C = collect(terms(w))
      for c in C
        m2[i, indexin([collect(monomials(c))[1]], bsp)[1]] += collect(coeffs(c))[1]
      end
    end
    push!(coll, m2)
  end
  return coll
end

@doc Markdown.doc"""
    symmetric_power_representation(rep::T, d::Int) where T <: LinRep -> T
    symmetric_power_representation(prep::T, d::Int) where T <: ProjRep -> T

Given a `F`-representation of a finite group `E` on a finite dimensional
`F`-vector space`V`, where `F` is a splitting field of `E`, return the induced
`F`-representation of `E` on $\text{Sym}^dV$.

In the case of a projective representation `prep`, return the reduction of the
action on the `d`-symmetric power of `V` of any linear lift of `prep`.
"""
function symmetric_power_representation(rep::LinRep{S, T, U}, d::Int) where {S, T, U}
  @req d > 0 "d must be a positive integer"
  mr = matrix_representation(rep)::Vector{dense_matrix_type(U)}
  mrd = _action_symmetric_power(mr, d)::Vector{dense_matrix_type(U)}
  chid = symmetric_power(character_representation(rep), d)
  return _linear_representation(representation_ring(rep), mrd, chid)
end

function symmetric_power_representation(prep::ProjRep{S, T, U, V}, d::Int) where {S, T, U, V}
  lin = symmetric_power_representation(linear_lift(prep), d)
  return ProjRep{S, T, U ,V}(lin, associated_schur_cover(prep))
end

@doc Markdown.doc"""
    homogeneous_polynomial_representation(rep::T, d::Int)
                                                         where T <: LinRep -> T
    homogeneous_polynonial_representation(prep::T, d::Int)
                                                        where T <: ProjRep -> T

Given a `F`-representation of a finite group `E` on a finite dimensional
`F`-vector space`V`, where `F` is a splitting field of `E`, return the
induced `F`-representation of `E` on $\text{Sym}^dV^{\vee}$.

In the case of a projective representation `prep`, return the reduction
of the action on the `d`-homogeneous part of the polynomial algebra
associated to`V` of any linear lift of `prep`.
"""
homogeneous_polynomial_representation(rep::LinRep, d::Integer) =
                   symmetric_power_representation(dual_representation(rep), d)

function homogeneous_polynomial_representation(prep::ProjRep{S, T, U, V}, d::Integer) where {S, T, U, V}
  lin = homogeneous_polynomial_representation(linear_lift(prep), d)
  return ProjRep{S, T, U, V}(lin, associated_schur_cover(prep))
end

###############################################################################

### Exterior power representation

"""
Compute wedge product between vectors of a same vector space, used to compute
induced actions on exterior power.
"""
function _wedge_product(V::Vector{T}) where T
  if length(V) == 1
    v = V[1]
    _, dec = _decompose_in_standard_basis(v)
    w = [(l[1], [l[2]]) for l in dec if l[1] != 0]
    return w
  end
  v = popfirst!(V)
  wp = _wedge_product(copy(V))
  std_bas, dec = _decompose_in_standard_basis(v)
  bas = basis_exterior_power(std_bas, length(V)+1)
  coeffs = zeros(parent(dec[1][1]), length(bas))
  for l in dec
    for w in wp
      if l[2] in w[2]
        continue
      end
      w_co = copy(w[2])
      lw_co = pushfirst!(w_co, l[2])
      ba, mult = _in_basis(lw_co, bas)
      coeffs[indexin([ba], bas)[1]] += mult*w[1]*l[1]
    end
  end
  return _change_basis(coeffs, bas)
end

"""
If `mg` represents an action of a group `E` on a vector space `V`, return
an induced action of \bidwedge^t V according to the basis given by
`basis_exterior_power`.
"""
function _action_exterior_power(mr::Vector{AbstractAlgebra.Generic.MatSpaceElem{S}}, t::Int) where S
  n = ncols(mr[1])
  F = base_ring(mr[1])
  std_bas = _standard_basis(F, n)
  bsp = basis_exterior_power(std_bas, t)
  coll = eltype(mr)[]
  for m in mr
    m2 = zero_matrix(F, length(bsp), length(bsp))
    for i = 1:length(bsp)
      b = bsp[i]
      mb = [transpose(m)*bb for bb in b]
      sp = _wedge_product(mb)
      for bbb in sp
        m2[indexin([bbb[2]], bsp)[1], i] += bbb[1]
      end
    end
    push!(coll, transpose(m2))
  end
  return coll
end

@doc Markdown.doc"""
    exterior_power_representation(rep::T, t::Int) where T <: LinRep -> T
    exterior_power_representation(prep::T, t::Int) where T <: ProjRep -> T

Given a `F`-representation of a finite group `E` on a finite dimensional
`F`-vector space`V`, where `F` is a splitting field of `E`, return the induced
`F`-representation of `E` on $\bigwedge^tV$, where the choice of the basis on
this exterior power is given by `basis_exterior_power`.

In the case of a projective representation `prep`, return the reduction of the
action on the `t`_exterior power of `V` of any linear lift of `prep`.
"""
function exterior_power_representation(rep::LinRep{S, T, U}, t::Int) where {S, T, U}
  @req 1 <= t <= dimension_representation(rep) "t must be an integer between 1 and the dimension of the representation"
  mr = matrix_representation(rep)::Vector{dense_matrix_type(U)}
  mrt = _action_exterior_power(mr, t)::Vector{dense_matrix_type(U)}
  chit = exterior_power(character_representation(rep), t)
  return _linear_representation(representation_ring(rep), mrt, chit)
end

function exterior_power_representation(prep::ProjRep{S, T, U, V}, t::Int) where {S, T, U, V}
  @req 1 <= t <= dimension_linear_lift(prep) "t must be an integer between 1 and the dimension of the linear lift"
  lin = exterior_power_representation(linear_lift(prep), t)
  return ProjRep{S, T, U, V}(lin, associated_schur_cover(prep))
end

###############################################################################
#
# Projectively faithful representations
#
###############################################################################

"""
  Check whether a small group `G` admits faithful projective representations
  of dimension `dim`. If yes, it returns the necessary tools to create the
  projectively faithful representations of a Schur cover associated to those
  faithful projective representations
"""
function _has_pfr(G::T, dim::Int) where T <: Oscar.GAPGroup
  # we start by computing a Schur cover and we turn it into an Oscar object
  G_gap = G.X
  f_gap = GG.EpimorphismSchurCover(G_gap)
  H_gap = GG.Source(f_gap)
  n, p = ispower(GG.Size(H_gap))
  if isprime(p)
    fff_gap = GG.EpimorphismPGroup(H_gap, p)
    E_gap = fff_gap(H_gap)
  else
    fff_gap = GG.IsomorphismPermGroup(H_gap)
    E_gap = fff_gap(H_gap)
  end
  E = Oscar._get_type(E_gap)(E_gap)
  H = Oscar._get_type(H_gap)(H_gap)
  fff = inv(GAPGroupHomomorphism(H, E, fff_gap))
  f = GAPGroupHomomorphism(H, G, f_gap)
  pschur = compose(fff, f)
  @assert is_surjective(f)

  # now we need to enumerate of characters of E of degree dim, modulo
  # multiplication by a linear character and up to the bounds
  # fixed by `setup`
  RR = representation_ring(E)
  ct = character_table_underlying_group(RR)
  _Irr = irreducible_characters_underlying_group(RR)
  _L = Int[Int(degree(irr)) for irr in _Irr]
  ps = sortperm(_L)
  L, Irr = _L[ps], _Irr[ps]
  RR.irr = Irr
  L = Int[i for i in L if i <= dim]
  lins = Oscar.GAPGroupClassFunction[irr for irr in Irr if irr[1] == 1]
  keep_index = Int[]
  sum_index = Vector{Int}[]
  chars = Oscar.GAPGroupClassFunction[]
  local El::ElevCtx
  El = elevator(L, dim)
  for l in El
    chi = sum(Irr[l])
    if !is_projective(chi, pschur) || !is_faithful(chi, pschur)
      continue
    end
    if any(lin -> lin*chi in chars, lins) # afforded by representations with similar reductions
      continue
    end
    push!(chars, chi)
    keep_index = union(keep_index, l)
    push!(sum_index, l)
  end
  bool = length(keep_index) > 0
  return bool, RR, sum_index, pschur
end

@doc Markdown.doc"""
    faithful_projective_representations(o::Int, n::Int, dim::Int)
    faithful_projective_representations(G::Oscar.GAPGroup, dim::Int)
                                                             -> Vector{ProjRep}

Given a small group `G` of ID `[o,n]`, return a set of representatives
of equivalences classes of faithful projective representations of `G` of complex
dimension `dim`.
"""
function faithful_projective_representations(G::Oscar.GAPGroup, dim::Int)
  bool, RR, sum_index, p = _has_pfr(G, dim)
  F = base_field(RR)
  E = underlying_group(RR)
  pfr = ProjRep{typeof(F), typeof(E), elem_type(F), typeof(p)}[]
  if !bool
    return pfr
  end
  Irr = irreducible_characters_underlying_group(RR)
  for l in sum_index
    chara = sum(Irr[l])
    @assert chara[1] == dim
    lr = affording_representation(RR, chara)
    pr = ProjRep{typeof(F), typeof(E), elem_type(F), typeof(p)}(lr, p)
    push!(pfr, pr)
  end
  return pfr
end

faithful_projective_representations(o::Int, n::Int, d::Int) = faithful_projective_representations(small_group(o, n), d)

###############################################################################
#
# Action on submodules, complements and quotients
#
###############################################################################

# All this part could be made cleaner with a proper structure of module which
# endows a vector space and a linear representation of a group on it

@doc Markdown.doc"""
    basis_isotypical_component(rep::LinRep{S, T, U},
                               chi::Oscar.GAPGroupClassFunction)
                                          where {S, T, U} -> Vector{MatElem{U}}
                                                      
Given an irreducible character `chi` which is constituent of the character of
`rep`, return a basis for the isotypical component of `rep` associated to `chi`.
The output is a set of bases each of them corresponding to a copy of a module
affording `chi` in the module associated to `rep`.

It is given in term of coordinates in the standard basis of the underlying
vector space of `rep`.
"""
function basis_isotypical_component(rep::LinRep{S, T, U}, chi::Oscar.GAPGroupClassFunction) where {S, T, U}
  @req Oscar.is_irreducible(chi) "chi must be an irreducible character"
  !is_constituent(character_representation(rep), chi) && return dense_matrix_type(U)[]
  RR = representation_ring(rep)
  irr = irreducible_affording_representation(RR, chi)
  As = matrix_representation(rep)
  Bs = matrix_representation(irr)
  return Hecke._basis_of_commutator_algebra(As, Bs)::Vector{dense_matrix_type(U)}
end

@doc Markdown.doc"""
    to_equivalent_block_representation(rep::LinRep{S, T, U})
                                 where {S, T, U} -> LinRep{S, T, U}, MatElem{U}
    to_equivalent_block_representation(rep::ProjRep{S, T, U, V})
                          where {S, T, U, V} -> ProjRep(S, T, U, V}, MatElem{U}
                                                         
Given a linear representation `rep` of a group, return an equivalent
representation `rep2` whose matrix representation is given by block in the
irreducible subrepresentations of `rep`, together with the matrix of change of
coordinates $rep \to rep2$.

In the case of a projective representation `prep`, return the reduction of an
equivalent block representation of any linear lift of `prep`.

Note that this can be expensive if the representations have irreducible
subrepresentations of high degree.
"""
function to_equivalent_block_representation(rep::LinRep{S, T, U}) where {S, T, U}
  RR = representation_ring(rep)
  cd = character_decomposition(character_representation(rep))
  F = base_field(representation_ring(rep))
  B = zero_matrix(F, 0, dimension_representation(rep))
  for (d, chi) in cd
    Bchi = basis_isotypical_component(rep, chi)
    @assert length(Bchi) == d
    Bd = reduce(vcat, Bchi)
    a = Int(degree(chi))
    @assert size(Bd) == (a*d, dimension_representation(rep))
    B = vcat(B, Bd)
  end
  @assert det(B) != 0
  mr = matrix_representation(rep)
  mr2 = dense_matrix_type(U)[B*m*inv(B) for m in mr]
  rep2 = _linear_representation(RR, mr2)
  return rep2, inv(B)
end

function to_equivalent_block_representation(prep::ProjRep{S, T, U, V}) where {S, T, U, V}
  lin, B = to_equivalent_block_representation(linear_lift(prep))
  return ProjRep{S, T, U, V}(lin,  associated_schur_cover(prep)), B
end

@doc Markdown.doc"""
    is_submodule(rep::LinRep{S, T, U}, M::W)
                                       where {S, T, U, W <: MatElem{U}} -> Bool

Given a matrix `M`, return whether the rows of `M` define a submodule of the
underlying module defined by the representation `rep`. If this is the case, one
can see the rows of `M` as the coordinates of the basis vectors of its underlying
vector space, with respect to the standard basis of the underlying vector space of
`rep`.
"""
function is_submodule(rep::LinRep{S, T, U}, M::W) where {S, T, U, W <: MatElem{U}}
  @req ncols(M) == dimension_representation(rep) "Incompatible size for the basis"
  mr = matrix_representation(rep)
  for m in mr
    ok = can_solve(M, M*m, side=:left)
    ok || return false
  end
  return true
end

@doc Markdown.doc"""
    action_on_submodule(rep::LinRep{S, T, U}, M::W)
                            where {S, T, U, W <: MatElem{U}} -> LinRep{S, T, U}
                                                                             
Given a matrix `M` whose rows define the coordinates of a submodule of `rep` in
the standard basis of the underlying vector space of `rep`, return the induced
action on the submodule defined by `M`, as a linear representation in the
representation ring of `rep`.
"""
function action_on_submodule(rep::LinRep{S, T, U}, M::W) where {S, T, U, W <: MatElem{U}}
  @req ncols(M) == dimension_representation(rep) "Incompatible size for the basis"
  mr = matrix_representation(rep)
  coll = eltype(mr)[]
  for m in mr
    ok, mm = can_solve_with_solution(M, M*m, side=:left)
    @req ok "Matrix does not define a submodule"
    push!(coll, mm)
  end
  return _linear_representation(representation_ring(rep), coll)
end

@doc Markdown.doc"""
    complement_submodule(rep::LinRep{S, T, U}, M::W)
                                     where {S, T, U, W <: MatElem{U}} -> W

Given a matrix `M` whose rows define the coordinates of a submodule of `rep` in
the standard basis of the underlying vector space of `rep`, return a matrix whose
rows span a complement of `M` inside the underlying module of `rep`. This is done
by looking complement of each isotypical components of the modules defined by `M`
inside the corresponding isotypical components of `rep`.
"""
function complement_submodule(rep::LinRep{S, T, U}, M::W) where {S, T, U, W <: MatElem{U}}
  rep2 = action_on_submodule(rep, M)
  cd = character_decomposition(rep)
  F = base_field(representation_ring(rep))
  bas = zero_matrix(F, 0, dimension_representation(rep))
  for (_, chi) in cd
    B1 = basis_isotypical_component(rep, chi)
    if !is_constituent(character_representation(rep2), chi)
      bas = vcat(bas, reduce(vcat, B1))
    else
      d = Int(degree(chi))
      B2 = basis_isotypical_component(rep2, chi)
      _K = zero_matrix(F, d*length(B1), d*length(B2))
      for j in 1:length(B2)
        B2j = B2[j]
        B2jM = (B2j*M)[1,:]
        B1u = reduce(vcat, [BB[1,:] for BB in B1])
        _Kj = solve_left(B1u, B2jM)
        for i in 1:d, k in 1:length(B1)
          _K[(k-1)*d + i, i + (j-1)*d] = _Kj[1, k]
        end
      end
      a, K2 = left_kernel(_K)
      Borth = K2*reduce(vcat, B1)
      @assert is_submodule(rep, Borth)
      bas = vcat(bas, Borth)
    end
  end
  return bas
end

@doc Markdown.doc"""
    quotient_representation(rep::LinRep{S, T, U}, M::W)
                  where {S, T, U, W <: MatElem{U} -> LinRep{S, T, U}, MatElem{U}

Given a matrix `M` whose rows define the coordinates of a submodule of `rep` in
the standard basis of the underlying vector space of `rep`, return the induced
action of the underlying group on the quotient modules `Q`, given as a linear
representation in the representation ring of `rep`. It is given with a matrix
defining the projection `rep -> Q`, in the standard basis of the underlying
vector space of `rep`.
"""
function quotient_representation(rep::LinRep{S, T, U}, M::W) where {S, T, U, W <: MatElem{U}}
  @req is_submodule(rep, M) "M is not invariant"
  F = base_field(representation_ring(rep))
  V = VectorSpace(F, dimension_representation(rep))
  sub_gene = [V(vec(collect(M[i,:]))) for i in 1:nrows(M)]
  Vsub, _ = sub(V, sub_gene)
  _, p = quo(V, Vsub)
  proj = p.matrix::W
  mr = matrix_representation(rep)
  coll = eltype(mr)[]
  for m in mr
    ok, mm = can_solve_with_solution(proj, m*proj, side=:right)
    @assert ok
    push!(coll, mm)
  end
  repQ = _linear_representation(representation_ring(rep), coll)
  @assert is_constituent(character_representation(rep), character_representation(repQ))
  return repQ, proj
end

@doc Markdown.doc"""
    isotypical_components(rep::LinRep{S, T, U}) where {S, T, U}
                         -> Dict{Oscar.GAPGroupClassFunction,
                                 Tuple{MatElem{U}, MatElem{U}}}

Given a linear representation `rep` in the representation ring `RR` associated
to a finite group `E` over a field `F`, return a dictionnary whose keys are the
irreducible `F`-characters of `E` which are components of the `F`-character of `E`
afforded by `rep`. To each key character `chi` is associated a matrix representing
the projection map corresponding to the isotypical component of `rep` associated
to the character `chi`. The matrices are given in the standard bases of the
underlying vector spaces of `rep` and of the corresponding isotypical component.
"""
function isotypical_components(rep::LinRep{S, T, U}) where {S, T, U}
  ic = Dict{Oscar.GAPGroupClassFunction, Tuple{dense_matrix_type(U), dense_matrix_type(U)}}()
  cd = character_decomposition(rep)
  F = base_field(representation_ring(rep))
  V = VectorSpace(F, dimension_representation(rep))
  for c in cd
    _B = basis_isotypical_component(rep ,c[2])
    @assert length(_B) == c[1]
    B = reduce(vcat, _B)
    B2 = complement_submodule(rep, B)
    _, pc = quo(V, sub(V, [V(vec(collect(B2[i,:]))) for i in 1:nrows(B2)])[1])
    ic[c[2]] = (B, pc.matrix)
  end
  return ic
end

