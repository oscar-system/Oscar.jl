if !isdefined(Main, :_prepare_scalar_types)
  function _prepare_scalar_types()
    NF, sr2 = quadratic_field(2)
    Qx, x = QQ[:x]
    K, a = Hecke.embedded_field(NF, real_embeddings(NF)[2])
    KK, (a1, a2) = embedded_number_field([x^2 - 2, x^3 - 5], [(0, 2), (0, 2)])
    return [(f, elem_type(f)) for f in (QQ, K, KK)]
  end
end
