# this file is meant to contain the structs for the collections in the database

################################################################################
# Leech Pair
@doc raw"""
    LeechPair

Container type for the entries in the Leech pair database, computed by G. Höhn
and G. Mason [HM16](@ref).
See also https://arxiv.org/abs/1505.06420v3 and the ancillary files.

Each `LeechPair` consists of a pair $(L, G)$ where $L$ is a negative definite
$\mathbb{Z}$-lattice isometric to the Leech lattice, and $G$ is a saturated
finite subgroup of the Conway group $Co_0 = O(L)$. Here ``saturated'' means
that $G$ is the pointwize stabilizer of the invariant lattice

```math
L^G := \{v \in L | g(v) = v, \forall g \in G\}.
```
"""
struct LeechPair
  leech::ZZLat
  G::MatrixGroup{QQFieldElem, QQMatrix}
  _id::Int
  rank::Int 
  order::ZZRingElem
  alpha::Int
  icoinvbar::Int
  iinvbar::Int
  ind::Int
  hinv::Int
  N::Int
  type::String

  LeechPair(args...) = new(args...)
end

###############################################################################
#
# Accessors and methods
#
###############################################################################

@doc raw"""
    lattice(LG::LeechPair) -> ZZLat

Return the underlying Leech lattice ``L``.
"""
Oscar.lattice(LG::LeechPair) = LG.leech

@doc raw"""
    group(LG::LeechPair) -> MatrixGroup{QQFieldElem, QQMatrix}

Return the underlying matrix group ``G``.
"""
Oscar.group(LG::LeechPair) = LG.G

@doc raw"""
    group_description(LG::LeechPair) -> String

Return a description of the underlying matrix group ``G``. See also
[`describe`](@ref).
"""
group_description(LG::LeechPair) = describe(LG.G)

@doc raw"""
    number(LG::LeechPair) -> Int

Return the number, in Höhn-Mason list, of the associated invariant lattice
$L^G$.

Refer to Table 1 of [HM16](@cite).
"""
number(LG::LeechPair) = LG._id

@doc raw"""
    rank_invariant_lattice(LG::LeechPair) -> Int

Return the rank of the associated invariant lattice $L^G$.

Refer to Table 1 of [HM16](@cite).
"""
rank_invariant_lattice(LG::LeechPair) = LG.rank

@doc raw"""
    group_order(LG::LeechPair) -> ZZRingElem

Return the order of the underlying matrix group ``G``.

Refer to Table 1 of [HM16](@cite).
"""
group_order(LG::LeechPair) = LG.order

@doc raw"""
    alpha(LG::LeechPair) -> Int

Return the difference between the rank of the associated invariant lattice
$L^G$ and the length of its discriminant group $D_{L^G}$.

Refer to Table 1 of [HM16](@cite).
"""
alpha(LG::LeechPair) = LG.alpha

@doc raw"""
    index_image_discriminant_representation_coinvariant(LG::LeechPair) -> Int

Return the index of the image of the map $O(L_G) \to O(D_{L_G})$ inside
$O(D_{L_G})$, where $L_G$ is the associated coinvariant lattice.

Refer to Table 1 of [HM16](@cite).
"""
index_image_discriminant_representation_coinvariant(LG::LeechPair) = LG.icoinvbar

@doc raw"""
    index_image_discriminant_representation_invariant(LG::LeechPair) -> Int

Return the index of the image of the map $O(L^G) \to O(D_{L^G})$ inside
$O(D_{L^G})$, where $L^G$ is the associated coinvariant lattice.

Refer to Table 1 of [HM16](@cite).
"""
index_image_discriminant_representation_invariant(LG::LeechPair) = LG.iinvbar

@doc raw"""
    index_normalizer_modulo_group(LG::LeechPair) -> Int

Return the index of the quotient $N_{O(L)}(G)/G$ inside $O(L^G)$ where $L^G$
is the associated invariant lattice and $N_{O(L)}(G)$ is the normalizer of $G$
in the Conway group $O(L) \cong Co_0$.

Refer to Table 1 of [HM16](@cite).
"""
index_normalizer_modulo_group(LG::LeechPair) = LG.ind

@doc raw"""
    class_number_invariant_lattice(LG::LeechPair) -> Int

Return the number of isometry classes in the genus of the associated invariant
lattice $L^G$.

Refer to Table 1 of [HM16](@cite).
"""
class_number_invariant_lattice(LG::LeechPair) = LG.hinv

@doc raw"""
    number_of_niemeier_embeddings(LG::LeechPair) -> Int

Return the number of (isometry classes of) Niemeier lattices with $(-2)$-roots
in which the associated coinvariant sublattice $L_G$ embeds.

Refer to Table 1 of [HM16](@cite).
"""
number_of_niemeier_embeddings(LG::LeechPair) = LG.N

@doc raw"""
    case_type(LG::LeechPair) -> String

Return the case type of the Leech pair, according to the last column in Table 1
of [HM16](@cite).
"""
case_type(LG::LeechPair) = LG.type

@doc raw"""
    invariant_lattice(LG::LeechPair) -> ZZLat

Return the associated invariant lattice $L^G$.
"""
Oscar.invariant_lattice(LG::LeechPair) = invariant_lattice(lattice(LG), group(LG))

@doc raw"""
    coinvariant_lattice(LG::LeechPair) -> ZZLat

Return the associated coinvariant lattice $L_G$.
"""
Oscar.coinvariant_lattice(LG::LeechPair) = first(coinvariant_lattice(lattice(LG), group(LG)))

@doc raw"""
    invariant_coinvariant_pair(LG::LeechPair) -> ZZLat, ZZLat

Return the pair $(L^G, L_G)$ consisting of the associated invariant and
coinvariant lattices respectively.
"""
Oscar.invariant_coinvariant_pair(LG::LeechPair) = invariant_coinvariant_pair(lattice(LG), group(LG))[1:2]

###############################################################################
#
# Printing
#
###############################################################################

function Base.show(io::IO, LG::LeechPair)
  print(io, "Leech pair no. $(number(LG))")
end
