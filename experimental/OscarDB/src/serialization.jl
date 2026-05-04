function load_object(s::DeserializerState, ::Type{Union{Nothing, T}}, params::Any, key::Symbol) where T
  set_key(s, key)
  isnothing(s.obj) && return nothing
  return load_object(s, T, params)
end

function load_object(s::DeserializerState, ::Type{Union{Nothing, T}}, key::Symbol) where T
  set_key(s, key)
  isnothing(s.obj) && return nothing
  return load_object(s, T)
end

@register_serialization_type LeechPair
type_params(x::LeechPair) = TypeParams(LeechPair, group(x))

function save_object(s::SerializerState, LG::LeechPair)
  save_data_dict(s) do
    save_object(s, Oscar.number(LG), :id)
    save_object(s, Oscar.rank_invariant_lattice(LG), :rank)
    save_object(s, Oscar.group_order(LG), :order)
    save_object(s, Oscar.alpha(LG), :alpha)
    save_object(s, Oscar.index_image_discriminant_representation_coinvariant(LG), :icoinvbar)
    save_object(s, Oscar.index_image_discriminant_representation_invariant(LG), :iinvbar)
    save_object(s, Oscar.index_normalizer_modulo_group(LG), :ind)
    save_object(s, Oscar.class_number_invariant_lattice(LG), :hinv)
    save_object(s, Oscar.number_of_niemeier_embeddings(LG), :N)
    # want to avoid using type key here to avoid confusion with file format
    save_object(s, Oscar.case_type(LG), :leech_pair_type) 
  end
end

function load_object(s::DeserializerState, ::Type{LeechPair}, G::MatGroup)
  db = Oscar.OscarDB.get_db()
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
    save_object(s, tsc._id, :name)
    save_object(s, simplicial_complex(tsc), :complex)
    save_object(s, dim(tsc), :dim)
    save_object(s, n_vertices(tsc), :n_vertices)
    save_object(s, f_vector(tsc), :f_vector)
    save_object(s, betti_numbers(tsc), :betti_numbers)
    save_object(s, automorphism_group(tsc), :aut_group)
    save_object(s, topological_type(tsc), :topological_type)
  end
end

function load_object(s::DeserializerState, ::Type{TransitiveSimplicialComplex},
                     params::PermGroup)
  load_node(s) do _
    TransitiveSimplicialComplex(
      load_object(s, String, :name),
      load_object(s, SimplicialComplex, :complex),
      load_object(s, Int, :dim),
      load_object(s, Int, :n_vertices),
      load_object(s, Vector{Int}, :f_vector),
      load_object(s, Vector{Int}, :betti_numbers),
      params,
      load_object(s, String, :topological_type)
    )
  end
end

################################################################################
# Phylogenetic collection structs
@register_serialization_type SmallTreeModel

type_params(stm::SmallTreeModel) = TypeParams(
  SmallTreeModel,
  group_based_phylogenetic_model(stm)
)

function save_object(s::SerializerState, stm::SmallTreeModel)
  save_data_dict(s) do
    save_object(s, stm._id, :model_encoding)
    save_object(s, stm.model_type, :model_type)
  end
end

function load_object(s::DeserializerState, ::Type{SmallTreeModel}, GBM::GroupBasedPhylogeneticModel)
  return SmallTreeModel(
    load_object(s, String, :model_encoding),
    GBM,
    load_object(s, String, :model_type)
  )
end

@register_serialization_type SmallGroupBasedModel

function type_params(sgbm::SmallGroupBasedModel) 
  I = vanishing_ideal(sgbm)
  if isnothing(I)
    OT = Nothing
  else
    OT = typeof(ordering(generating_system(I)))
  end
  return TypeParams(
    SmallGroupBasedModel,
    :model => group_based_phylogenetic_model(sgbm),
    :ideal_ordering_type => TypeParams(OT, nothing)
  )
end

function save_object(s::SerializerState, sgbm::SmallGroupBasedModel)
  save_data_dict(s) do
    save_object(s, sgbm._id, :model_encoding)
    save_object(s, sgbm.model_type, :model_type)
    save_object(s, n_leaves(sgbm), :n_leaves)
    save_object(s, level(sgbm), :level)
    save_object(s, n_cycles(sgbm), :n_cycles)
    save_object(s, dim(sgbm), :dimension)
    save_object(s, degree(sgbm), :degree)
    save_object(s, n_coordinates(sgbm), :n_coordinates)
    save_object(s, dimension_singular_locus(sgbm), :dimension_singular_locus)
    save_object(s, degree_singular_locus(sgbm), :degree_singular_locus)
    save_object(s, euclidean_distance_degree(sgbm), :euclidean_distance_degree)
    save_object(s, parametrization(sgbm), :parametrization)
    save_object(s, equivalent_classes(sgbm), :equivalent_classes)
    # here we store the vanishing ideal as it's generating systems as it will
    # store more information about the ideal
    I = vanishing_ideal(sgbm)
    if !isnothing(I)
      save_object(s, generating_system(I), :vanishing_ideal)
    else
      save_object(s, nothing, :vanishing_ideal)
    end
  end
end

function load_object(s::DeserializerState, ::Type{SmallGroupBasedModel}, params::Dict)
  GBM = params[:model]
  dom = base_ring(model_ring(GBM)[1])
  codom = parameter_ring(GBM)[1]
  fmr = base_ring(full_model_ring(GBM)[1])
  n_leaves = load_object(s, Int, :n_leaves)

  if params[:ideal_ordering_type] == Nothing
    I = nothing
  else
    I = ideal(load_object(s, Oscar.IdealGens, Dict(:base_ring => dom, :ordering_type => params[:ideal_ordering_type]), :vanishing_ideal))
  end
  return SmallGroupBasedModel(
    load_object(s, String, :model_encoding),
    GBM,
    load_object(s, String, :model_type),
    n_leaves,
    load_object(s, Int, :level),
    load_object(s, Int, :n_cycles),
    load_object(s, Union{Int, Nothing}, :dimension),
    load_object(s, Union{Int, Nothing}, :degree),
    load_object(s, Int, :n_coordinates),
    load_object(s, Union{Int, Nothing}, :dimension_singular_locus),
    load_object(s, Union{Int, Nothing}, :degree_singular_locus),
    load_object(s, Union{Int, Nothing}, :euclidean_distance_degree),
    load_object(s, Oscar.MPolyAnyMap, Dict(:domain => dom, :codomain => codom), :parametrization),
    load_object(s, Dict{Tuple{[Int for _ in 1:n_leaves]...}, Vector{MPolyRingElem}}, Dict(:key_params => nothing, :value_params => fmr), :equivalent_classes),
    I
  )
end

@register_serialization_type SmallPhylogeneticModel

function type_params(spm::SmallPhylogeneticModel)
  I = vanishing_ideal(spm)
  if isnothing(I)
    OT = Nothing
  else
    # fix at some point when we have groebner basis here
    OT = Nothing
    #OT = typeof(ordering(generating_system(I)))
  end
  return TypeParams(
    SmallPhylogeneticModel,
    :model => phylogenetic_model(spm),
    :ideal_ordering_type => TypeParams(OT, nothing)
  )
end

function save_object(s::SerializerState, spm::SmallPhylogeneticModel)
  save_data_dict(s) do
    save_object(s, spm._id, :model_encoding)
    save_object(s, spm.model_type, :model_type)
    save_object(s, n_leaves(spm), :n_leaves)
    save_object(s, level(spm), :level)
    save_object(s, n_cycles(spm), :n_cycles)
    save_object(s, dim(spm), :dimension)
    save_object(s, degree(spm), :degree)
    save_object(s, n_coordinates(spm), :n_coordinates)
    save_object(s, dimension_singular_locus(spm), :dimension_singular_locus)
    save_object(s, degree_singular_locus(spm), :degree_singular_locus)
    save_object(s, maximum_likelihood_degree(spm), :maximum_likelihood_degree)
    save_object(s, euclidean_distance_degree(spm), :euclidean_distance_degree)
    save_object(s, parametrization(spm), :parametrization)
    save_object(s, equivalent_classes(spm), :equivalent_classes)
    # here we store the vanishing ideal as it's generating systems as it will
    # store more information about the ideal
    I = vanishing_ideal(spm)
    # update when we can have GB
    save_object(s, I, :vanishing_ideal)
    #if !isnothing(I)
    #  save_object(s, generating_system(I), :vanishing_ideal)
    #else
    #  save_object(s, nothing, :vanishing_ideal)
    #end
  end
end

function load_object(s::DeserializerState, ::Type{SmallPhylogeneticModel}, params::Dict)
  SPM = params[:model]
  dom = base_ring(model_ring(SPM)[1])
  codom = parameter_ring(SPM)[1]
  fmr = base_ring(full_model_ring(SPM)[1])
  n_leaves = load_object(s, Int, :n_leaves)
  return SmallPhylogeneticModel(
    load_object(s, String, :model_encoding),
    SPM,
    load_object(s, String, :model_type),
    n_leaves,
    load_object(s, Int, :level),
    load_object(s, Int, :n_cycles),
    load_object(s, Union{Int, Nothing}, :dimension),
    load_object(s, Union{Int, Nothing}, :degree),
    load_object(s, Int, :n_coordinates),
    load_object(s, Union{Int, Nothing}, :dimension_singular_locus),
    load_object(s, Union{Int, Nothing}, :degree_singular_locus),
    load_object(s, Union{Int, Nothing}, :maximum_likelihood_degree),
    load_object(s, Union{Int, Nothing}, :euclidean_distance_degree),
    load_object(s, Oscar.MPolyAnyMap, Dict(:domain => dom, :codomain => codom), :parametrization),
    load_object(s, Dict{Tuple{[Int for _ in 1:n_leaves]...}, Vector{MPolyRingElem}}, Dict(:key_params => nothing, :value_params => fmr), :equivalent_classes),
    
    load_object(s, Union{Nothing, Ideal}, dom, :vanishing_ideal)
  )
end
