@register_serialization_type LeechPair
type_params(x::LeechPair) = TypeParams(LeechPair, group(x))

function save_object(s::SerializerState, LG::LeechPair)
  save_data_dict(s) do
    save_object(s, Oscar.number(LG), :id)
    save_object(s, Oscar.rank_invariant_lattice(LG), :rank)
    save_object(s, Oscar.order_group(LG), :order)
    save_object(s, Oscar.alpha(LG), :alpha)
    save_object(s, Oscar.index_image_in_Oq_coinvariant(LG), :icoinvbar)
    save_object(s, Oscar.index_image_in_Oq_invariant(LG), :iinvbar)
    save_object(s, Oscar.index_normaliser_modulo_group(LG), :ind)
    save_object(s, Oscar.class_number_invariant_lattice(LG), :hinv)
    save_object(s, Oscar.number_of_niemeier_embeddings(LG), :N)
    # want to avoid using type key here to avoid confusion with file format
    save_object(s, LG.type, :leech_pair_type) 
  end
end

function load_object(s::DeserializerState, ::Type{LeechPair}, G::MatrixGroup)
  # dev until we move collection to offical db
  db = Oscar.OscarDB.get_db(;dev=true)
  leech = Oscar.OscarDB.find_one(db["zzlattices"], Dict("_id" => "leech"))

  LeechPair(
    leech,
    G,
    load_object(s, Int, :id),
    load_object(s, Int, :rank),
    load_object(s, ZZRingElem, :order),
    load_object(s, Int, :alpha),
    load_object(s, Int, :icoinvbar),
    load_object(s, Int, :iinvbar),
    load_object(s, Int, :ind),
    load_object(s, Int, :hinv),
    load_object(s, Int, :N), 
    load_object(s, String, :leech_pair_type),
  )
end

@register_serialization_type TransitiveSimplicialComplex

type_params(tsc::TransitiveSimplicialComplex) = TypeParams(
  TransitiveSimplicialComplex,
  automorphism_group(tsc)
)

function save_object(s::SerializerState, tsc::TransitiveSimplicialComplex)
  save_data_dict(s) do 
    save_object(s, simplicial_complex(tsc), :simplicial_complex)
    save_object(s, topological_type(tsc), :topological_type)
  end
end

function load_object(s::DeserializerState, ::Type{TransitiveSimplicialComplex},
                     params::PermGroup)
  load_node(s) do _
    TransitiveSimplicialComplex(
      load_object(s, SimplicialComplex, :simplicial_complex),
      params,
      load_object(s, String, :topological_type)
    )
  end
end
