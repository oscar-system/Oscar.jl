

function embedding_in_unimodular_manyroots(S::ZZLat, pos::Int, neg::Int; primitive=true, even=true,compute_overlattices=false)
  #this version computes a root lattice of the biggest rank possible represented by the genus of R
@vprintln :Lattice 1 "computing embedding in L_$(n)"
  pS, kS, nS = signature_tuple(S)
  @req kS == 0 "S must be non-degenerate"
  even || throw(NotImplementedError("for now we need the unimodular lattice to be even."))
  pR = pos - pS
  nR = neg - nS
  rkR=pR+nR
  DS = discriminant_group(S)
  DR = rescale(DS, -1)  # discriminant group of R = S^\perp in L as predicted by Nikulin
  GR = genus(DR, (pR, nR)) # genus of R
  roots=biggest_root_sublattice(GR);
  #now compute the actual embedding
  #first thing: try if it embeds primitively
  #em=primitive_embeddings(GR, roots; classification=:none) #maybe not necessary to do this as first step
  R=[]
  #if em[1]==true
    emb=primitive_embeddings(GR, roots; classification=:first)
    R = emb[2][1][1]
if compute_overlattices #compute all overlattices of "roots" at the same time
    M = ZZLat[]
    #println("overlattice timing")
    #@time 
    overs=Oscar._overlattice_orbits(roots)
    #lov=length(overs)
    for ov in overs
      em2=primitive_embeddings(GR, ov; classification=:first)
      if em2[1]==true
        M=em2[2][1][1]
        break
      end
    end
    R=M
  else #compute overlattices of increasing order of "roots", until one embeds primitively in GR
    d = ZZ(det(roots))
    D = discriminant_group(roots)
    idD = hom(D,D,gens(D))
    G,iG = image_in_Oq(roots)
    orders = [i for i in divisors(d) if divides(d,i^2)[1]]
    M = ZZLat[]
    for ord in orders 
    #@show ord, D
      b, l, p = is_prime_power_with_data(ord)
      if b && is_elementary(D, p)
        sg = first.(first.(_isotropic_subspaces_representatives_and_stabilizers_elementary(D, iG, valuation(ord,p);do_stab=false)))
      else 
      # slooow
        sg = domain.(first.(Oscar._subgroups_orbit_representatives_and_stabilizers(idD, G, ord)))
      end
      for S in sg 
        M = cover(S)
        if !is_integral(M) || (even && !is_even(M))
          continue
        end 
        em2=primitive_embeddings(GR, M; classification=:first)
        if em2[1]==true
          M=em2[2][1][1]
          break
        end
      end
      if M!=[]
        R=M
        break
      end
    end
  end
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
  #find a (negative definite) root sublattice of maximal rank represented by the genus g
  n=signature_pair(g)[2]
  for i=n:-1:1
    allroots=root_lattices(i)
    l=length(allroots)
    for ll=l:-1:1 #starts with the one with smallest discriminant
      RR=rescale(allroots[ll],-1)
      gr=genus(RR)
      if represents(g,gr)
        return RR
      end
    end
  end
end