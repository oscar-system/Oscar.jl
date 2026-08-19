@doc raw"""
    embedding_in_unimodular_manyroots(S::ZZLat, pos::Int, neg::Int; primitive=true, even=true) -> ZZLat, ZZLat, AbstractSpaceMor, ZZLat
    

Compute an embedding $i_S \colon S \to L$ into some even unimodular lattice of the given signature such that the orthogonal complement $R:=i_S(S)^\perp$ has a root sublattice of the largest possible rank. The return value is the quadruple `(L, S', iS, R)`.

See also  [`embed_in_unimodular`](@ref).
  
"""
function embedding_in_unimodular_manyroots(S::ZZLat, pos::Int, neg::Int; primitive=true, even=true)
  @req iszero(mod(pos - neg,8)) "an even unimodular lattice of signature ($pos, $neg) does not exist"
  @vprintln :Lattice 1 "computing embedding in L_$(n)"
  pS, kS, nS = signature_tuple(S)
  @req kS == 0 "S must be non-degenerate"
  even || throw(NotImplementedError("for now we need the unimodular lattice to be even."))
  pR = pos - pS
  nR = neg - nS
  DS = discriminant_group(S)
  DR = rescale(DS, -1)  # discriminant group of R = S^\perp in L as predicted by Nikulin
  GR = genus(DR, (pR, nR)) # genus of R
  roots = biggest_root_sublattice(GR)
  R = _overlattice_orbits(roots, GR)[1] #don't trust the name of this function! 
  R = lll(R;  same_ambient=false)  # make R a bit nicer
  SR, inj = direct_sum(S, R)
  iS, iR = inj
  V = ambient_space(SR)
  S = lattice(V, basis_matrix(S) * iS.matrix)
  R = lattice(V, basis_matrix(R) * iR.matrix)
  fl, glue = is_anti_isometric_with_anti_isometry(discriminant_group(S), discriminant_group(R))
  @assert fl
  L = overlattice(glue)
  @assert V === ambient_space(SR)
  @hassert :Lattice 1 abs(det(L)) ==1
  @hassert :Lattice 1 is_integral(L)==1
  @hassert :Lattice 1 !even || iseven(L)
  return L, S, iS, R
end

function biggest_root_sublattice(g::ZZGenus) 
  # find a (negative definite) root sublattice of maximal rank represented by the genus g
  n = signature_pair(g)[2]
  for i = n:-1:1
    allroots = root_lattices(i)
    l = length(allroots)
    for ll = l:-1:1 # starts with the one with smallest discriminant
      RR = rescale(allroots[ll], -1)
      gr = genus(RR)
      if represents(g,gr)
        return RR
      end
    end
  end
end



