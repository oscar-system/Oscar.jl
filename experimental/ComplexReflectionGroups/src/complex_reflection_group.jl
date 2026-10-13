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

function _common_field(fields; field_presentation::Symbol=:native)
  existing_data = _existing_common_field(fields)
  data = if existing_data !== nothing
    existing_data
  else
    cyclotomic_data = _cyclotomic_common_field(fields)
    cyclotomic_data === nothing ? _absolute_common_field(fields) : cyclotomic_data
  end
  common_field, embeddings = data
  if field_presentation == :absolute &&
      !(common_field isa Union{QQField, AbsSimpleNumField})
    absolute_field, to_presentation = absolute_simple_field(common_field)
    to_absolute = inv(to_presentation)
    return absolute_field, [compose(f, to_absolute) for f in embeddings]
  end
  return common_field, embeddings
end

function _block_diagonal_product(
    groups::Vector{MatGroup}; fixed_space_dimension::Int=0,
    field_presentation::Symbol=:native
  )
  @assert fixed_space_dimension >= 0

  common_field, embeddings = if isempty(groups)
    QQ, Any[]
  else
    _common_field(base_ring.(groups); field_presentation)
  end
  mapped_groups = [
    map_entries(embedding, group) for (embedding, group) in zip(embeddings, groups)
  ]
  dimensions = Int[degree(group) for group in mapped_groups]
  total_dimension = fixed_space_dimension + sum(dimensions; init=0)
  product_generators = elem_type(
    matrix_space(common_field, total_dimension, total_dimension)
  )[]

  for (i, group) in enumerate(mapped_groups)
    for generator in gens(group)
      # Use the decomposition V^G \oplus V_1 \oplus ... \oplus V_r.  The fixed
      # block is omitted when it has dimension zero so that the old essential
      # models retain exactly their previous matrices.
      blocks = fixed_space_dimension == 0 ?
        [identity_matrix(common_field, dimension) for dimension in dimensions] :
        [
          identity_matrix(common_field, fixed_space_dimension),
          [identity_matrix(common_field, dimension) for dimension in dimensions]...,
        ]
      blocks[i + (fixed_space_dimension == 0 ? 0 : 1)] = matrix(generator)
      push!(product_generators, block_diagonal_matrix(blocks))
    end
  end

  return matrix_group(
    common_field, total_dimension, product_generators
  ), embeddings, mapped_groups
end

"""
    complex_reflection_group(G::ComplexReflectionGroupType, model::Symbol=:default;
                             field_presentation=:native)
    complex_reflection_group(; fixed_space_dimension=0, model=:default,
                             field_presentation=:native)
    complex_reflection_group(i::Int, model::Symbol=:default;
                             fixed_space_dimension=0, field_presentation=:native)
    complex_reflection_group(m::Int, p::Int, n::Int,
                             model::Symbol=:default;
                             fixed_space_dimension=0, field_presentation=:native)
    complex_reflection_group(types::Vector, model::Symbol=:default;
                             fixed_space_dimension=0, field_presentation=:native)

Construct an exact matrix representative of a known complex reflection group
type. The available named models are `:CHEVIE`, `:LT`, and `:Magma`;
`:default` chooses a supported model componentwise. The nontrivial components
are put into a common exact coefficient field and their original presentation
fields and embeddings are retained. A fixed summand is represented by an
identity block preceding the nontrivial component blocks.

The `:LT` models preserve the standard Hermitian form and follow the
orthonormal realizations in [LT09](@cite), Chapters 6--8; see the
[Lehrer--Taylor models](@ref crg-lt-models) for precise constructions and
coefficient fields. The arithmetic `:Magma` models need not preserve the
standard form.

With `field_presentation=:native`, retain a tower presentation when possible.
With `field_presentation=:absolute`, transport the same matrices to an absolute
simple number field (or `QQ`); this does not change their basis or minimize
their field of definition. The aliases `:LT_flat` and `Symbol("LT-flat")`
select `:LT` with an absolute presentation. Model labels remain `:LT`, and
[`components`](@ref) and
[`complex_reflection_group_component_embeddings`](@ref) retain the original
component realizations and their exact maps into the resulting field.

This is a type-aware constructor, not a recognition routine for arbitrary
matrix groups. Pure trivial types, including the zero-dimensional type, are
represented by a generator-free matrix group over `QQ`.
"""
function complex_reflection_group(
    G::ComplexReflectionGroupType, model::Symbol=:default;
    field_presentation::Symbol=:native
  )

  @req field_presentation in (:native, :absolute) "Unknown coefficient-field presentation"
  if model in (:LT_flat, Symbol("LT-flat"))
    model = :LT
    field_presentation = :absolute
  end
  @req model in (:default, :CHEVIE, :LT, :Magma) "Specified model not found"

  # this will be the list of matrix groups corresponding to the components of G
  component_groups = MatGroup[]

  # list of models of the components
  modellist = Symbol[]

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

  fixed_dimension = fixed_space_dimension(G)
  if length(component_groups) == 1 && fixed_dimension == 0 &&
      (field_presentation == :native ||
       base_ring(component_groups[1]) isa Union{QQField, AbsSimpleNumField})
    # treat this as a special case because direct_product([G]) will have type direct
    # product which looks weird for a single group.
    return component_groups[1]
  else
    matgrp, component_embeddings, mapped_component_groups =
      _block_diagonal_product(
        component_groups; fixed_space_dimension=fixed_dimension, field_presentation
      )

    set_attribute!(matgrp, :order, order(G))
    set_attribute!(matgrp, :is_complex_reflection_group, true)
    set_attribute!(matgrp, :complex_reflection_group_type, G)
    set_attribute!(matgrp, :complex_reflection_group_model, modellist)
    set_attribute!(matgrp, :reflection_hyperplane_orbit_marking_is_known, true)
    set_attribute!(matgrp, :is_irreducible, is_irreducible(G))
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

    @assert degree(matgrp) == degree(G)
    return matgrp
  end

end

"""
    components(G::MatGroup)

Return the original nontrivial component realizations retained by a complex
reflection group model. For an irreducible nontrivial marked realization,
return `[G]`. A constructor realization of a pure trivial type returns an empty
vector. Return `nothing` when no component decomposition is supplied.

The component groups keep their own presentation fields. Use
[`complex_reflection_group_component_embeddings`](@ref) for the exact maps
from those fields to the common coefficient field of `G`.
"""
function components(G::MatGroup)
  if has_attribute(G, :complex_reflection_group_components)
    return get_attribute(G, :complex_reflection_group_components)
  end
  group_type = complex_reflection_group_type(G)
  if group_type !== nothing &&
      is_irreducible(group_type) &&
      !is_trivial(group_type)
    return [G]
  end
  return nothing
end

"""
    complex_reflection_group_component_embeddings(G::MatGroup)

Return the exact coefficient-field embeddings from the original nontrivial
component realizations of `G` into `base_ring(G)`, or `nothing` if no component
decomposition is supplied. A constructor realization of a pure trivial type
returns an empty vector.
"""
function complex_reflection_group_component_embeddings(G::MatGroup)
  if has_attribute(G, :complex_reflection_group_component_embeddings)
    return get_attribute(G, :complex_reflection_group_component_embeddings)
  end
  group_type = complex_reflection_group_type(G)
  if group_type !== nothing &&
      is_irreducible(group_type) &&
      !is_trivial(group_type)
    return [id_hom(base_ring(G))]
  end
  return nothing
end

# Convenience constructors
complex_reflection_group(
    ; fixed_space_dimension::IntegerUnion=0, model::Symbol=:default,
    field_presentation::Symbol=:native
  ) = complex_reflection_group(
    ComplexReflectionGroupType(; fixed_space_dimension), model; field_presentation
  )

complex_reflection_group(
    i::Int, model::Symbol=:default; fixed_space_dimension::IntegerUnion=0,
    field_presentation::Symbol=:native
  ) = complex_reflection_group(
    ComplexReflectionGroupType(i; fixed_space_dimension), model; field_presentation
  )

complex_reflection_group(
    m::Int, p::Int, n::Int, model::Symbol=:default;
    fixed_space_dimension::IntegerUnion=0, field_presentation::Symbol=:native
  ) = complex_reflection_group(
    ComplexReflectionGroupType(m, p, n; fixed_space_dimension), model; field_presentation
  )

complex_reflection_group(
    X::Vector, model::Symbol=:default; fixed_space_dimension::IntegerUnion=0,
    field_presentation::Symbol=:native
  ) = complex_reflection_group(
    ComplexReflectionGroupType(X; fixed_space_dimension), model; field_presentation
  )


############################################################################################
# Getter functions
############################################################################################
"""
    complex_reflection_group_type(G::MatGroup)

Return the assigned complex reflection group type of `G`, or
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
    fixed_space_dimension(G::MatGroup)

Return the fixed-space dimension recorded by the assigned complex reflection
group type of `G`, or `nothing` if no type is assigned.

This getter does not compute a fixed space or recognize the type from the
matrices of `G`.
"""
function fixed_space_dimension(G::MatGroup)
  group_type = complex_reflection_group_type(G)
  return group_type === nothing ? nothing : fixed_space_dimension(group_type)
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
  K = base_ring(W)
  n = degree(W)
  matspace = matrix_space(K, n, n)
  dual_generators = elem_type(matspace)[]
  for w in gens(W)
    push!(dual_generators, matspace(transpose(matrix(w^-1))))
  end
  WD = matrix_group(K, n, dual_generators)
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

  # Only propagate an explicitly retained component decomposition.  The
  # fallback `components(W) == [W]` for an irreducible group must not recurse.
  if has_attribute(W, :complex_reflection_group_components)
    component_groups = get_attribute(W, :complex_reflection_group_components)
    dual_component_groups = complex_reflection_group_dual.(component_groups)
    component_embeddings = complex_reflection_group_component_embeddings(W)
    set_attribute!(
      WD,
      :complex_reflection_group_components,
      dual_component_groups,
    )
    set_attribute!(
      WD,
      :complex_reflection_group_component_embeddings,
      component_embeddings,
    )
    set_attribute!(
      WD,
      :complex_reflection_group_mapped_components,
      [
        map_entries(embedding, group) for
        (embedding, group) in zip(component_embeddings, dual_component_groups)
      ],
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
