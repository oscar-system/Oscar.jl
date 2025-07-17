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
end

###############################################################################
#
# Accessors and methods
#
###############################################################################

Oscar.lattice(LG::LeechPair) = LG.leech

Oscar.group(LG::LeechPair) = LG.G

group_description(LG::LeechPair) = describe(LG.G)

number(LG::LeechPair) = LG._id

rank_invariant_lattice(LG::LeechPair) = LG.rank

order_group(LG::LeechPair) = LG.order

alpha(LG::LeechPair) = LG.alpha

index_image_in_Oq_coinvariant(LG::LeechPair) = LG.icoinvbar

index_image_in_Oq_invariant(LG::LeechPair) = LG.iinvbar

index_normaliser_modulo_group(LG::LeechPair) = LG.ind

class_number_invariant_lattice(LG::LeechPair) = LG.hinv

number_of_niemeier_embeddings(LG::LeechPair) = LG.N

Oscar.type(LG::LeechPair) = LG.type

Oscar.invariant_lattice(LG::LeechPair) = invariant_lattice(lattice(LG), group(LG))

Oscar.coinvariant_lattice(LG::LeechPair) = coinvariant_lattice(lattice(LG), group(LG))

Oscar.invariant_coinvariant_pair(LG::LeechPair) = invariant_coinvariant_pair(lattice(LG), group(LG))

###############################################################################
#
# Printing
#
###############################################################################

function Base.show(io::IO, LG::LeechPair)
  print(io, "Leech pair no. $(number(LG))")
end

function Base.show(io::IO, ::MIME"text/plain", LG::LeechPair)
  if get(io, :supercompact, false)
    print(io, "Leech pair")
  else
    println(io, "Leech pair no. $(number(LG))")
    println(io, "Rank:  $(rank_invariant_lattice(LG))")
    println(io, "Order: $(order_group(LG))")
    print(io, "Type:  $(type(LG))")
  end
end

@register_serialization_type LeechPair
type_params(x::LeechPair) = TypeParams(LeechPair, group(x))

function save_object(s::SerializerState, LG::LeechPair)
  save_data_dict(s) do
    save_object(s, group(LG), :group)
    save_object(s, number(LG), :_id)
    save_object(s, rank_invariant_lattice(LG), :rank)
    save_object(s, index_image_in_Oq_coinvariant(LG), :icoinvbar)
    save_object(s, index_image_in_Oq_invariant(LG), :iinvbar)
    save_object(s, index_normaliser_modulo_group(LG), :hinv)
    save_object(s, number_of_niemeier_embeddings(LG), :N)
    # want to avoid using type key here to avoid confusion with file format
    save_object(s, number_of_niemeier_embeddings(LG), :leech_pair_type) 
  end
end

function load_object(s::DeserializerState, ::Type::LeechPair, G::MatrixGroup)
  leech = 
  LeechPair()
end
