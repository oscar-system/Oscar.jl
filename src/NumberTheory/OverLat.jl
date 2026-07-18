# compute representatives of the O(L) orbits 
# of overlattices L < M such that M and L have the same roots
function overlattices_no_new_roots(L::ZZLat; even::Bool=true)
  d = ZZ(det(L))
  primes = [p for (p,v) in factor(d) if v>1]
  result = ZZLat[L]
  for p in primes
    result2 = ZZLat[]
    for M in result
      tmp = _overlattices_no_new_roots_index_prime_power(M, p)
      append!(result2, tmp)
    end
    result = result2
  end
  return result
end

function _overlattices_no_new_roots_index_prime_power(L::ZZLat, p::ZZRingElem; even::Bool=true)
# computing overlattices of index p^n recursively, one risks to find the same overlattice multiple times by following different paths: 
# for each new overlattice T, we have to check if it is isomorphic to the ones already found (of the same index as T).
# We set up a dictionary of invariant functions: thus, we have to compare T only to overlattices having its same invariant function.
  unique_iso_classes = Dict()
  unique_iso_classes[oscar_invariant_function(L)] = [L]
  todo = [L]
  while !isempty(todo)
    M = pop!(todo)
    S = _overlattices_no_new_roots_index_p(M, p)
    for T in S
      kT = oscar_invariant_function(T)
      unique_isometry_classes_with_invariant_k = get!(unique_iso_classes, kT, ZZLat[])
      if !any(is_isometric(T,x) for x in unique_isometry_classes_with_invariant_k)
        push!(todo, T)
        push!(unique_isometry_classes_with_invariant_k, T)
      end
    end
  end
  return reduce(append!, values(unique_iso_classes); init=ZZLat[]) # put the values of the dict into one big list
end

function _overlattices_no_new_roots_index_p(L::ZZLat, p::ZZRingElem; even::Bool=true, lines::Bool=true)
  D = discriminant_group(L)
  timesp = hom(D, D, [p*x for x in gens(D)])
  Kp,iKp = kernel(timesp)
  OKp = orthogonal_group(Kp)
  qq = gram_matrix_quadratic(Kp)

  G, iG = image_in_Oq(L) #kinda slow, probably inevitable
  Gp,ip = restrict_automorphism_group(G, iKp; check=false)
  b,j = is_subgroup(Gp, OKp; check=false); @assert b  
  result = ZZLat[]
  nrootsL = length(short_vectors(L, 2))
  @assert OKp === codomain(j)
  if lines==false
    @vprintln :ZZLatWithIsom 1 "computing isotropic subspace orbits by permutation representation (GAP)"
    sg = first.(first.(Oscar._isotropic_subspaces_representatives_and_stabilizers_elementary(Kp, j, 1; do_stab=false)))
    @vprintln :ZZLatWithIsom 1 "found $(length(sg))"
  else @vprintln :ZZLatWithIsom 1 "computing isotropic subspace orbits by line orbits (Hecke)"
    sg = _isotropic_lines_representatives(p, Kp, Gp)
    @vprintln :ZZLatWithIsom 1 "found $(length(sg))"
  end
  for S in sg
    M = cover(S)
    if (even && p==2 && !is_even(M))
      continue
    end
    nrootsM = length(short_vectors(M,2))
    if nrootsL == nrootsM
      push!(result, M)
    end
  end
  return result
end

function _isotropic_lines_representatives(p::ZZRingElem, Kp::TorQuadModule, Gp::AutomorphismGroup{TorQuadModule})
  #input: a prime p, a p-elementary torsion quadratic module Kp, a subgroup Gp of Aut(K)
  @assert is_elementary(Kp,p)
  F = GF(p)
  Gp_gens = gens(Gp)
  Gp_gens_Hecke = FqMatrix[]
  for g in Gp_gens
    M = matrix(g)
    append!(Gp_gens_Hecke, [matrix(F, collect(M))])
  end
  lineorbs = Hecke.line_orbits(Gp_gens_Hecke)
  sg = []
  for l in lineorbs
    lambda = map_entries(x->lift(ZZ,x), transpose(matrix(l[1]))) #the orbit representative, as a ZZmatrix
    lambda_in_Kp, ii = sub(Kp, [sum(gen(Kp,j)*lambda[1,j] for j in 1:ncols(lambda); init=TorQuadModuleElem(Kp,zero(abelian_group(Kp))))]) #as submodule of Kp
    if  is_zero(gram_matrix_quadratic(lambda_in_Kp))
      append!(sg, [lambda_in_Kp])
    end
  end
  return sg
  end
