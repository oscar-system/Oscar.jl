# This file implements explicit models of symplectic reflection groups.

@doc raw"""
    symplectic_doubling(G::MatGroup{T}) where T <: QQAlgFieldElem

Return the cotangent, or symplectic, doubling of the matrix group `G`.

If `G` acts on ``\mathfrak h`` by matrices ``g``, the returned group acts on
``\mathfrak h \oplus \mathfrak h^*`` by

```math
g^\star = \begin{pmatrix}g&0\\0&g^{-\mathsf T}\end{pmatrix}.
```

The result retains `G` as its source, together with the two block embeddings
and the standard alternating form

```math
\Omega = \begin{pmatrix}0&I\\-I&0\end{pmatrix}.
```

If `G` already has a complex reflection group type, that type and the chosen
model are retained as *source* metadata.  They are deliberately not assigned
as the complex reflection group type and model of the doubled representation.
"""
function symplectic_doubling(G::MatGroup{T}) where {T <: QQAlgFieldElem}
  K = base_ring(G)
  n = degree(G)
  matspace = matrix_space(K, 2*n, 2*n)
  generators = elem_type(matspace)[]

  for g in gens(G)
    doubled_generator =
      matspace(block_diagonal_matrix([matrix(g), transpose(matrix(g^-1))]))
    push!(generators, doubled_generator)
  end

  doubled_group = matrix_group(generators)

  identity_block = identity_matrix(K, n)
  primal_embedding = zero_matrix(K, n, 2*n)
  primal_embedding[1:n, 1:n] = identity_block
  dual_embedding = zero_matrix(K, n, 2*n)
  dual_embedding[1:n, n + 1:2*n] = identity_block

  omega = zero_matrix(K, 2*n, 2*n)
  omega[1:n, n + 1:2*n] = identity_block
  omega[n + 1:2*n, 1:n] = -identity_block

  source_type = complex_reflection_group_type(G)
  source_model = complex_reflection_group_model(G)
  set_attribute!(
    doubled_group,
    :order,
    source_type === nothing ? order(G) : order(source_type),
  )
  set_attribute!(doubled_group, :symplectic_doubling_source, G)
  set_attribute!(
    doubled_group,
    :symplectic_doubling_block_embeddings,
    (primal_embedding, dual_embedding),
  )
  set_attribute!(doubled_group, :symplectic_form, alternating_form(omega))

  if source_type !== nothing
    set_attribute!(doubled_group, :symplectic_doubling_source_type, source_type)
  end
  if source_model !== nothing
    set_attribute!(doubled_group, :symplectic_doubling_source_model, source_model)
  end

  source_is_known_reflection_group =
    source_type !== nothing ||
    (has_attribute(G, :is_complex_reflection_group) &&
     get_attribute(G, :is_complex_reflection_group))
  if source_is_known_reflection_group
    set_attribute!(doubled_group, :is_symplectic_reflection_group, true)
  end

  return doubled_group
end

"""
    symplectic_doubling_source(G::MatGroup)

Return the matrix group whose cotangent representation produced `G`, or
`nothing` if `G` is not a marked symplectic doubling.
"""
function symplectic_doubling_source(G::MatGroup)
  if has_attribute(G, :symplectic_doubling_source)
    return get_attribute(G, :symplectic_doubling_source)
  end
  return nothing
end

@doc raw"""
    symplectic_doubling_block_embeddings(G::MatGroup)

Return the row-space inclusions of ``\mathfrak h`` and ``\mathfrak h^*`` into
the marked decomposition ``\mathfrak h \oplus \mathfrak h^*`` of `G`, or
`nothing` if no such decomposition is supplied.
"""
function symplectic_doubling_block_embeddings(G::MatGroup)
  if has_attribute(G, :symplectic_doubling_block_embeddings)
    return get_attribute(G, :symplectic_doubling_block_embeddings)
  end
  return nothing
end

"""
    symplectic_form(G::MatGroup)

Return the alternating form supplied with the symplectic matrix group `G`, or
`nothing` if no form is supplied.
"""
function symplectic_form(G::MatGroup)
  if has_attribute(G, :symplectic_form)
    return get_attribute(G, :symplectic_form)
  end
  return nothing
end

"""
    symplectic_doubling_source_type(G::MatGroup)

Return the complex reflection group type assigned to the source of the marked
symplectic doubling `G`, or `nothing` if no source type is supplied.
"""
function symplectic_doubling_source_type(G::MatGroup)
  if has_attribute(G, :symplectic_doubling_source_type)
    return get_attribute(G, :symplectic_doubling_source_type)
  end
  source = symplectic_doubling_source(G)
  return source === nothing ? nothing : complex_reflection_group_type(source)
end

"""
    symplectic_doubling_source_model(G::MatGroup)

Return the model assigned to the source of the marked symplectic doubling `G`,
or `nothing` if no source model is supplied.
"""
function symplectic_doubling_source_model(G::MatGroup)
  if has_attribute(G, :symplectic_doubling_source_model)
    return get_attribute(G, :symplectic_doubling_source_model)
  end
  source = symplectic_doubling_source(G)
  return source === nothing ? nothing : complex_reflection_group_model(source)
end

"""
    symplectic_reflection_group(W::MatGroup)

Return the symplectic doubling of the complex reflection group `W`.

An assigned complex reflection group type supplies the conjugacy-type metadata
and avoids recognition of `W` from its matrices. Any reflection-class marking
is separate realization data.
"""
function symplectic_reflection_group(W::MatGroup)
  if complex_reflection_group_type(W) === nothing &&
     !is_complex_reflection_group(W)
    throw(ArgumentError("group must be a complex reflection group"))
  end

  doubled_group = symplectic_doubling(W)
  set_attribute!(doubled_group, :is_symplectic_reflection_group, true)
  return doubled_group
end

"""
    is_symplectic_reflection_group(G::MatGroup)

Return whether `G` is known to be a symplectic reflection group.

Recognition of an unmarked matrix group is not implemented.
"""
function is_symplectic_reflection_group(G::MatGroup)
  if has_attribute(G, :is_symplectic_reflection_group)
    return get_attribute(G, :is_symplectic_reflection_group)
  end
  error("recognition of an unmarked symplectic reflection group is not implemented")
end
