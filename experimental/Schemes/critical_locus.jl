function ideal_sheaf_of_critical_locus(phi::AbsCoveredSchemeMorphism)
  phi_cov = covering_morphism(phi)
  dom_cov = domain(phi_cov)
  cod_cov = codomain(phi_cov)
  ideal_dict = IdDict{AbsAffineScheme, Ideal}()
  for U in patches(dom_cov)
    phi_loc = phi_cov[U]
    V = codomain(phi_loc)
    y = gens(OO(V))
    f = pullback(phi_loc).(y)
    df = jacobi_matrix(f)
    r = length(f)
    I = _degeneracy_locus(df, dim(codomain(phi)))
    ideal_dict[U] = radical(I)
  end
  IdealSheaf(domain(phi), ideal_dict, check=false)
end

function _degeneracy_locus(df::MatrixElem{T}, r::Int) where {T<:MPolyRingElem}
  R = base_ring(df)
  return ideal(R, minors(df, r))
end

function _degeneracy_locus(df::MatrixElem{T}, r::Int) where {T<:MPolyLocRingElem}
  R = base_ring(df)
  return ideal(R, minors(df, r))
end

function _degeneracy_locus(df::MatrixElem{T}, r::Int) where {T<:MPolyQuoRingElem}
  A = base_ring(df)
  R = base_ring(A)::MPolyRing
  M = map_entries(lift, dF)
  g = gens(modulus(A))
  dg = jacobi_matrix(g)
  M = hcat(M, dg)
  s = r + ngens(R) - dim(A)
  return ideal(A, minors(M, s))
end

function _degeneracy_locus(df::MatrixElem{T}, r::Int) where {T<:MPolyQuoLocRingElem}
  L = base_ring(df)
  R = localized_ring(L)::MPolyLocRing
  M = map_entries(lift, df)
  g = gens(modulus(L))
  dg = jacobi_matrix(g)
  M = hcat(M, dg)
  s = r + ngens(R) - dim(L)
  return ideal(L, minors(M, s))
end

