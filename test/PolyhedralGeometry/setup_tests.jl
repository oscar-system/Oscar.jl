if !isdefined(Main, :_prepare_scalar_types)
  function _prepare_scalar_types()
    NF, sr2 = quadratic_field(2)
    Qx, x = QQ[:x]
    K, a = Hecke.embedded_field(NF, real_embeddings(NF)[2])
    KK, (a1, a2) = embedded_number_field([x^2 - 2, x^3 - 5], [(0, 2), (0, 2)])
    return [(f, elem_type(f)) for f in (QQ, K, KK)]
  end

  function _check_im_perm_rows(inc::IncidenceMatrix, o)
    oinc = incidence_matrix(o)
    nr, nc = size(inc)
    (nr, nc) == size(oinc) &&
      issetequal(Polymake.row.(Ref(inc), 1:nr),
                 Polymake.row.(Ref(oinc), 1:nr))
  end

  _matrix_from_property(b::SubObjectIterator{<:Union{LinearHalfspace, LinearHyperplane}}) = permutedims(hcat([normal_vector(be) for be in b]...))

  _matrix_from_property(b::SubObjectIterator{<:Union{AffineHalfspace, AffineHyperplane}}) = permutedims(hcat([vcat(-negbias(be), normal_vector(be)) for be in b]...))

  # only used for cones that are linear halfspaces
  _matrix_from_property(b::SubObjectIterator{Cone{T}}) where T = _matrix_from_property(SubObjectIterator{LinearHalfspace{T}}(b.Obj, b.Acc, b.n))

  _matrix_from_property(b::SubObjectIterator) = permutedims(hcat(b...))

  _oscar_matrix_from_property(a, b::SubObjectIterator) = matrix(a, _matrix_from_property(b))

  function _polymake_matrix_from_property(b::SubObjectIterator)
    m = _matrix_from_property(b)
    return Polymake.Matrix{Oscar._scalar_type_to_polymake(eltype(m))}(m)
  end
end
