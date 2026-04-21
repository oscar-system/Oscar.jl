@doc raw"""
  embedding_in_unimodular_manyroots(S::ZZLat, pos::Int, neg::Int; primitive=true, even=true) -> ZZLat, ZZLat, AbstractSpaceMor{Hecke.QuadSpace{QQField, QQMatrix}, QQMatrix}, ZZLat
  This function is a modified version of embed_in_unimodular.
  It produces an embedding of a lattice S in an even unimodular lattice L of signature (pos, neg); the genus g of the orthogonal complement to S is unique, its discriminant form being opposite to that of S.
  The lattice R is chosen in the genus g as to maximize its root sublattice. R is computed using the function _overlattice_orbits(:ZZLat, :ZZGenus) (see the comments there).
  The output consists of: a lattice L' isometric to L, obtained as primitive extension of S+R; S; a primitive embedding of S in L'; R.
  
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
  R = lll(R)  # make R a bit nicer
  R = integer_lattice(; gram=gram_matrix(R), cached=false) # clear the history of R
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
  @hassert :Lattice 1 denominator(gram_matrix(L))==1
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



