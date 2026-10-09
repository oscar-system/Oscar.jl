# This file implements explicit models of complex reflection groups.
#
# References:
#
# * Lehrer, G. I., & Taylor, D. E. (2009). Unitary reflection groups (Vol. 20, p. viii). Cambridge University Press, Cambridge.
#
# * Thiel, U. (2014). On restricted rational Cherednik algebras. TU Kaiserslautern.
#
# * Marin, I., & Michel, J. (2010). Automorphisms of complex reflection groups. Represent. Theory, 14, 747–788.
#
# Ulrich Thiel, 2023

function _natural_field_embedding(K, L)
  if K === L
    return id_hom(K)
  elseif K isa QQField
    return hom(QQ, L)
  elseif K isa AbsSimpleNumField
    return hom(K, L, L(gen(K)))
  elseif K isa Hecke.RelSimpleNumField
    base_embedding = _natural_field_embedding(base_field(K), L)
    base_embedding === nothing && return nothing
    return hom(K, L, base_embedding, L(gen(K)))
  end
  return nothing
end

function _existing_common_field(fields)
  candidates = sort(
    collect(fields);
    by=K -> K isa QQField ? 1 : absolute_degree(K),
    rev=true,
  )
  for common_field in candidates
    embeddings = Any[]
    for K in fields
      embedding = try
        _natural_field_embedding(K, common_field)
      catch
        nothing
      end
      embedding === nothing && break
      push!(embeddings, embedding)
    end
    length(embeddings) == length(fields) && return common_field, embeddings
  end
  return nothing
end

function _cyclotomic_common_field(fields)
  cyclotomic_orders = Int[]
  for K in fields
    if K isa QQField
      push!(cyclotomic_orders, 1)
      continue
    end
    if !(K isa AbsSimpleNumField)
      return nothing
    end
    is_cyclotomic, conductor = Hecke.is_cyclotomic_type(K)
    if !is_cyclotomic
      return nothing
    end
    push!(cyclotomic_orders, Int(conductor))
  end

  conductor = foldl(lcm, cyclotomic_orders; init=1)
  if conductor <= 2
    common_field = QQ
    embeddings = [
      K isa QQField ? id_hom(QQ) :
      hom(K, QQ, QQ(is_one(gen(K)) ? 1 : -1)) for K in fields
    ]
    return common_field, embeddings
  end

  common_field, root = cyclotomic_field(conductor)
  embeddings = Any[]
  for (K, order) in zip(fields, cyclotomic_orders)
    if K isa QQField
      push!(embeddings, hom(QQ, common_field))
    else
      push!(
        embeddings,
        hom(K, common_field, root^div(conductor, order)),
      )
    end
  end
  return common_field, embeddings
end

function _absolute_common_field(fields)
  absolute_fields = Any[]
  to_presentation_fields = Any[]

  for K in fields
    if K isa QQField
      continue
    end
    absolute_field, to_presentation = absolute_simple_field(K)
    push!(absolute_fields, absolute_field)
    push!(to_presentation_fields, to_presentation)
  end

  if isempty(absolute_fields)
    return QQ, Any[id_hom(QQ) for _ in fields]
  end

  common_field = absolute_fields[1]
  absolute_embeddings = Any[id_hom(common_field)]
  for absolute_field in absolute_fields[2:end]
    new_common_field, old_embedding, new_embedding =
      compositum(common_field, absolute_field)
    absolute_embeddings = [
      compose(embedding, old_embedding) for embedding in absolute_embeddings
    ]
    push!(absolute_embeddings, new_embedding)
    common_field = new_common_field
  end

  embeddings = Any[]
  nonrational_index = 0
  for K in fields
    if K isa QQField
      push!(embeddings, hom(QQ, common_field))
    else
      nonrational_index += 1
      push!(
        embeddings,
        compose(
          inv(to_presentation_fields[nonrational_index]),
          absolute_embeddings[nonrational_index],
        ),
      )
    end
  end

  return common_field, embeddings
end

function _common_field(fields)
  existing_data = _existing_common_field(fields)
  if existing_data !== nothing
    return existing_data
  end
  cyclotomic_data = _cyclotomic_common_field(fields)
  if cyclotomic_data !== nothing
    return cyclotomic_data
  end
  return _absolute_common_field(fields)
end

function _block_diagonal_product(groups::Vector{MatGroup})
  fields = base_ring.(groups)
  common_field, embeddings = _common_field(fields)
  mapped_groups = [
    map_entries(embedding, group) for (embedding, group) in zip(embeddings, groups)
  ]
  dimensions = degree.(mapped_groups)
  product_generators = elem_type(
    matrix_space(common_field, sum(dimensions), sum(dimensions))
  )[]

  for (i, group) in enumerate(mapped_groups)
    for generator in gens(group)
      blocks = [identity_matrix(common_field, dimension) for dimension in dimensions]
      blocks[i] = matrix(generator)
      push!(product_generators, block_diagonal_matrix(blocks))
    end
  end

  return matrix_group(product_generators), embeddings, mapped_groups
end

"""
    complex_reflection_group(G::ComplexReflectionGroupType, model::Symbol=:default)
    complex_reflection_group(i::Int, model::Symbol=:default)
    complex_reflection_group(m::Int, p::Int, n::Int, model::Symbol=:default)
    complex_reflection_group(types::Vector, model::Symbol=:default)

Construct an exact matrix representative of a known essential complex
reflection group type. The available named models are `:CHEVIE`, `:LT`, and `:Magma`;
`:default` chooses a supported model componentwise. For a reducible type the
components are put into a common exact coefficient field and their original
presentation fields and embeddings are retained.

This is a type-aware constructor, not a recognition routine for arbitrary
matrix groups. The empty essential type does not specify a positive ambient
dimension and has no matrix realization through this interface.
"""
function complex_reflection_group(G::ComplexReflectionGroupType, model::Symbol=:default)

  # this will be the list of matrix groups corresponding to the components of G
  component_groups = MatGroup[]

  # list of models of the components
  modellist = []

  for C in components(G)

    t = C.type[1]

    # Get default model for C
    if model == :default
      if isa(t,Int)
        Cmodel = :Magma
      else
        (m,p,n) = t
        if m == 1 && p == 1
          Cmodel = :CHEVIE
        else
          Cmodel = :Magma
        end
      end
    else
      Cmodel = model
    end

    if Cmodel == :LT
      matgrp = complex_reflection_group_LT(t)
    elseif Cmodel == :Magma
      matgrp = complex_reflection_group_Magma(t)
    elseif Cmodel == :CHEVIE
      matgrp = complex_reflection_group_CHEVIE(t)
    else
      error("Specified model not found")
    end

    # set attributes that are already known from type
    set_attribute!(matgrp, :order, order(C))
    set_attribute!(matgrp, :is_complex_reflection_group, true)
    set_attribute!(matgrp, :complex_reflection_group_type, C)
    set_attribute!(matgrp, :complex_reflection_group_model, [Cmodel])
    set_attribute!(matgrp, :reflection_hyperplane_orbit_marking_is_known, true)
    set_attribute!(matgrp, :is_irreducible, true)
    if Cmodel == :LT
      set_attribute!(
        matgrp,
        :invariant_hermitian_form,
        identity_matrix(base_ring(matgrp), degree(matgrp)),
      )
    end

    # add to list
    push!(component_groups, matgrp)
    push!(modellist, Cmodel)
  end

  if length(component_groups) == 1
    # treat this as a special case because direct_product([G]) will have type direct
    # product which looks weird for a single group.
    return component_groups[1]
  else
    @req !isempty(component_groups) "A concrete ambient space is required for the trivial reflection group type"
    matgrp, component_embeddings, mapped_component_groups =
      _block_diagonal_product(component_groups)

    set_attribute!(matgrp, :order, order(G))
    set_attribute!(matgrp, :is_complex_reflection_group, true)
    set_attribute!(matgrp, :complex_reflection_group_type, G)
    set_attribute!(matgrp, :complex_reflection_group_model, modellist)
    set_attribute!(matgrp, :reflection_hyperplane_orbit_marking_is_known, true)
    set_attribute!(matgrp, :is_irreducible, false)
    set_attribute!(matgrp, :complex_reflection_group_components, component_groups)
    set_attribute!(matgrp, :complex_reflection_group_mapped_components, mapped_component_groups)
    set_attribute!(matgrp, :complex_reflection_group_component_embeddings, component_embeddings)

    component_forms = invariant_hermitian_form.(component_groups)
    if all(form -> form !== nothing && is_one(form), component_forms)
      form = identity_matrix(base_ring(matgrp), degree(matgrp))
      if is_unitary(matgrp, form)
        set_attribute!(matgrp, :invariant_hermitian_form, form)
      end
    end

    return matgrp
  end

end

"""
    components(G::MatGroup)

Return the original component realizations retained by a reducible complex
reflection group model. For an irreducible marked realization, return `[G]`.
Return `nothing` when no component decomposition is supplied.

The component groups keep their own presentation fields. Use
[`complex_reflection_group_component_embeddings`](@ref) for the exact maps
from those fields to the common coefficient field of `G`.
"""
function components(G::MatGroup)
  if has_attribute(G, :complex_reflection_group_components)
    return get_attribute(G, :complex_reflection_group_components)
  end
  group_type = complex_reflection_group_type(G)
  if group_type !== nothing && is_irreducible(group_type)
    return [G]
  end
  return nothing
end

"""
    complex_reflection_group_component_embeddings(G::MatGroup)

Return the exact coefficient-field embeddings from the original component
realizations of `G` into `base_ring(G)`, or `nothing` if no component
decomposition is supplied.
"""
function complex_reflection_group_component_embeddings(G::MatGroup)
  if has_attribute(G, :complex_reflection_group_component_embeddings)
    return get_attribute(G, :complex_reflection_group_component_embeddings)
  end
  group_type = complex_reflection_group_type(G)
  if group_type !== nothing && is_irreducible(group_type)
    return [id_hom(base_ring(G))]
  end
  return nothing
end

# Convenience constructors
complex_reflection_group(i::Int, model::Symbol=:default) = complex_reflection_group(ComplexReflectionGroupType(i), model)

complex_reflection_group(m::Int, p::Int, n::Int, model::Symbol=:default) = complex_reflection_group(ComplexReflectionGroupType(m,p,n), model)

complex_reflection_group(X::Vector, model::Symbol=:default) = complex_reflection_group(ComplexReflectionGroupType(X), model)


############################################################################################
# Getter functions
############################################################################################
"""
    complex_reflection_group_type(G::MatGroup)

Return the assigned essential complex reflection group type of `G`, or
`nothing` if none is assigned. For an irreducible group, this means its
`GL(n, C)`-conjugacy class, represented by its normalized Shephard--Todd label.
This getter does not try to recognize an arbitrary matrix group.
"""
function complex_reflection_group_type(G::MatGroup)
  if has_attribute(G, :complex_reflection_group_type)
    return get_attribute(G, :complex_reflection_group_type)
  end
  return nothing
  # this should be upgraded later to work with a general matrix group (identifying the
  # type from scratch is not so easy though)
end

"""
    complex_reflection_group_model(G::MatGroup)

Return the realization-model labels assigned to the canonical components of
`G`, or `nothing` if none are assigned. Constructor models use `:CHEVIE`,
`:LT`, or `:Magma`; a group produced by
[`complex_reflection_group_dual`](@ref) uses `:dual` and retains the source
model separately.
"""
function complex_reflection_group_model(G::MatGroup)
  if has_attribute(G, :complex_reflection_group_model)
    return get_attribute(G, :complex_reflection_group_model)
  end
  return nothing
end

"""
    invariant_hermitian_form(G::MatGroup)

Return the Gram matrix of the invariant Hermitian form supplied with the
chosen model of the complex reflection group `G`, or `nothing` if that model
does not supply one.

This is realization data: every finite complex matrix group admits an
invariant Hermitian form, but it need not be the standard form and it is not a
property of the Shephard--Todd type alone.
"""
function invariant_hermitian_form(G::MatGroup)
  if has_attribute(G, :invariant_hermitian_form)
    return get_attribute(G, :invariant_hermitian_form)
  end
  return nothing
end

"""
    complex_reflection_group_dual(W::MatGroup)

Return the contragredient matrix realization of the complex reflection group
`W`, with generators `transpose(inv(g))` for `g` in `gens(W)`.

The complex reflection group type is retained and the realization model is
labelled `:dual`; use [`complex_reflection_group_dual_source`](@ref) and
[`complex_reflection_group_dual_source_model`](@ref) for its provenance. A
known reference reflection-class marking is transported explicitly. If `W`
has only an assigned type and no such marking, the result remains unmarked.
"""
function complex_reflection_group_dual(W::MatGroup)

  group_type = complex_reflection_group_type(W)
  if group_type === nothing && !is_complex_reflection_group(W)
    error("Group is not a complex reflection group")
  end
  WD = matrix_group([transpose(matrix(w^-1)) for w in gens(W)])
  group_model = complex_reflection_group_model(W)

  set_attribute!(WD, :order, group_type === nothing ? order(W) : order(group_type))
  set_attribute!(WD, :is_complex_reflection_group, true)
  set_attribute!(WD, :complex_reflection_group_dual_source, W)
  if group_model !== nothing
    set_attribute!(WD, :complex_reflection_group_dual_source_model, group_model)
  end
  if group_type !== nothing
    set_attribute!(WD, :complex_reflection_group_type, group_type)
    set_attribute!(
      WD,
      :complex_reflection_group_model,
      fill(:dual, number_of_components(group_type)),
    )
    set_attribute!(WD, :is_irreducible, is_irreducible(group_type))
  elseif has_attribute(W, :is_irreducible)
    set_attribute!(WD, :is_irreducible, get_attribute(W, :is_irreducible))
  end
  if invariant_hermitian_form(W) !== nothing && is_one(invariant_hermitian_form(W))
    set_attribute!(WD, :invariant_hermitian_form, invariant_hermitian_form(W))
  end

  # Transport the full reference marking, including the orientation of each
  # cyclic pointwise stabilizer. If s is the marked source reflection, then
  # transpose(s) = (s^-1)^(-T) has the same nontrivial eigenvalue in the dual
  # group, whereas the corresponding raw dual generator has the inverse one.
  if group_type !== nothing && _has_known_reflection_hyperplane_orbit_marking(W)
    source_orbits = reflection_hyperplane_orbits(W)
    dual_marking = [WD(transpose(matrix(representative(O)))) for O in source_orbits]
    set_reflection_hyperplane_orbit_marking!(WD, dual_marking)
  end

  component_groups = components(W)
  if component_groups !== nothing && length(component_groups) > 1
    set_attribute!(
      WD,
      :complex_reflection_group_components,
      complex_reflection_group_dual.(component_groups),
    )
    set_attribute!(
      WD,
      :complex_reflection_group_component_embeddings,
      complex_reflection_group_component_embeddings(W),
    )
  end

  return WD

end

"""
    complex_reflection_group_dual_source(G::MatGroup)

Return the source group of a realization constructed by
[`complex_reflection_group_dual`](@ref), or `nothing` if `G` has no such
provenance.
"""
function complex_reflection_group_dual_source(G::MatGroup)
  if has_attribute(G, :complex_reflection_group_dual_source)
    return get_attribute(G, :complex_reflection_group_dual_source)
  end
  return nothing
end

"""
    complex_reflection_group_dual_source_model(G::MatGroup)

Return the model labels of the source of a realization constructed by
[`complex_reflection_group_dual`](@ref), or `nothing` if no source model is
available.
"""
function complex_reflection_group_dual_source_model(G::MatGroup)
  if has_attribute(G, :complex_reflection_group_dual_source_model)
    return get_attribute(G, :complex_reflection_group_dual_source_model)
  end
  return nothing
end


###########################################################################################
# Checking if a matrix group is a complex reflection group.
###########################################################################################
function is_complex_reflection_group(G::MatGroup{T}) where T <: QQAlgFieldElem

  # First, check if we already know that G is a complex reflection group
  if has_attribute(G, :is_complex_reflection_group)
    return get_attribute(G, :is_complex_reflection_group)
  end

  # Now, check if the (fixed) generating set happens to consist of reflections.
  if all(is_complex_reflection, gens(G))
    set_attribute!(G, :is_complex_reflection_group, true)
    return true
  end

  # Last attempt: determine all reflections and see if they generate the group (this is
  # then equivalent to being a complex reflection group).
  refls = collect(complex_reflections(G))
  H,f = sub(G, [G(matrix(w)) for w in refls])
  if H == G
    set_attribute!(G, :is_complex_reflection_group, true)
    return true
  end

  return false
end

###########################################################################################
# Cartan matrix
###########################################################################################
function complex_reflection_group_cartan_matrix(W::MatGroup)

  if !is_complex_reflection_group(W)
    throw(ArgumentError("Group is not a complex reflection group"))
  end

  # We collect roots and coroots of the generators of W
  roots = []
  coroots = []

  for g in gens(W)
    b,g_data = is_complex_reflection_with_data(g)
    push!(roots, root(g_data))
    push!(coroots, coroot(g_data))
  end

  K = base_ring(W)
  n = length(roots)
  C = matrix(K,n,n,[ canonical_pairing(coroots[j], roots[i]) for i=1:n for j=1:n ])

  return C
end
