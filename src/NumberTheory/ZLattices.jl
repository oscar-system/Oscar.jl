function _overlattice_orbits(L::ZZLat; even=true)
  d = ZZ(det(L))
  D = discriminant_group(L)
  idD = hom(D,D,gens(D))
  G,iG = image_in_Oq(L)
  orders = [i for i in divisors(d) if divides(d,i^2)[1]]
  result = ZZLat[]
  for ord in orders 
    #@show ord, D
    b, l, p = is_prime_power_with_data(ord)
    if b && is_elementary(D, p)
      sg = first.(first.(_isotropic_subspaces_representatives_and_stabilizers_elementary(D, iG, valuation(ord,p);do_stab=false)))
    else 
      # slooow
      sg = domain.(first.(_subgroups_orbit_representatives_and_stabilizers(idD, G, ord)))
    end
    for S in sg 
      M = cover(S)
      if !is_integral(M) || (even && !is_even(M))
        continue
      end 
      push!(result,M)        
    end
  end
  return result
end

function root_overlattices(n::Int)
  result = ZZLat[]
  for R in root_lattices(n)
    for S in _overlattice_orbits(R)
      # only add the ones not adding new roots 
      RS = root_sublattice(S)
      if RS == R # this is terribly inefficient
        push!(result,S)
      end 
      if S != RS
        #@show "new"
      end
    end 
  end 
  return result
end



graph_hash(g::Graph{T}) where T<:Union{Directed,Undirected} = Polymake.graph.canonical_hash(Oscar.pm_object(g))

function invariant_function_graph_hash(L::ZZLat; max_size = 6000)
  @assert is_integral(L)
  if rank(L)==0 
    return UInt(0)
  end
  lb = minimum(L)
  ub = lb+2*scale(L)
  G = ZZ.(gram_matrix(L))
  # we just care about stuff modulo 2
  n = 0
  kk = GF(2)
  Gk = kk.(G)
  r = rank(L)
  sv = FqMatrix[]
  success = true
  gamma = graph(Undirected, 0)
  _v = zero_matrix(kk,1,r)
  tmp = zero_matrix(kk,r,1)
  tmp2 = zero_matrix(kk,1,1)
  cv = characteristic_vectors(L)
  if length(cv)<max_size
    for v in cv
      n = n+1
      for i in 1:r
        _v[1,i] = v[i]
        tmp[i,1] = v[i]
      end
      Gv = mul!(tmp,Gk,tmp)
      add_vertex!(gamma)
      for (i,w) in enumerate(sv)
        if !iszero(mul!(tmp2,w,Gv))
          add_edge!(gamma,i,n)
        end
      end
      push!(sv, _v)
    end
    return BigInt(graph_hash(gamma))
  end
  return BigInt(automorphism_group_order(L))
  #return Hecke.default_invariant_function(L)
end

function invariant_function_2_4(L::ZZLat; max_size = 10000)  
  v2 = short_vectors(L,2, Int)
  v4 = short_vectors(L,4,4,Int)
  G = Hecke._int_matrix_with_overflow(ZZ.(gram_matrix(L)),ZZ())
  m = MSet{MSet{Int}}()
  #Gi = zero_matrix(ZZ,rank(L),1)
  #v = zero_matrix(ZZ,rank(L),1)
  Gv = zeros(Int,rank(L))
  for (v,_) in v2
    Gv = G*v
    ii = MSet{Int}()
    for (j,_) in v4
      # put the absolute value
      # beause only v or -v is returned by the short vector functions 
      d = abs(dot(j,Gv))
      push!(ii,d)
    end
    push!(m, ii)
  end
  return m
end

send_lat(x::ZZLat) = collect(ZZ.(gram_matrix(x))) 
load_lat(x::Matrix{ZZRingElem}) = integer_lattice(gram=matrix(x),check=false)

aut_order(x) = automorphism_group_order(load_lat(x))

function _default_invariant_function(L::ZZLat)
  kn = kissing_number(L)::Int
  rlr, _ = root_lattice_recognition(L)
  R = BigInt
  m = R(ZZ(minimum(L)))
  d = R(ZZ(det(L)))
  ago = R(automorphism_group_order(L))
  if rank(L)>sum(i[2] for i in rlr;init=0)
    t = invariant_function_2_4(L; max_size = 10000)
  else
    t = multiset(multiset([0]))
  end 
  return (m, rlr, kn, ago, d, t)
end


function oscar_invariant_function(L::ZZLat)
  _invariants = Any[]
  _L = L
  while rank(_L) > 12
    M, P, _ = Hecke._shortest_vectors_sublattice(_L; check=false)
    i = index(P,M)
    push!(_invariants, (_default_invariant_function(rescale(P, 1//scale(P); cached=false)),i))
    _L = orthogonal_submodule(_L, P)
  end
  push!(_invariants, invariant_function_graph_hash(_L))
  push!(_invariants, invariant_function_2_4(L; max_size = 10000))
  push!(_invariants, BigInt(automorphism_group_order(L)))
  return Tuple(_invariants)
end


# return all characteristic vectors up to sign
# unfortunately still to many for a fast graph hash 
# at least in higher rank
# the idea follows https://arxiv.org/pdf/2004.14022
"""
    characteristic_vectors(L::ZZLat) -> Vector{ZZMatrix}
    
Return a set of characterisitc vectors of ``L`` up to sign.

The implementation follows ideas of  Sikirić Haensch, Voight and van Woerden.
"""
function characteristic_vectors(L::ZZLat)
  L = lattice(rational_span(L))
  S1,P1, v1  = Hecke._shortest_vectors_sublattice(L; check=false)
  cvL = v1
  #append!(cvL, [-i for i in cvL])
  B = coordinates(basis_matrix(S1), P1)
  A = abelian_group(ZZ.(B))
  BS1 = ZZ.(basis_matrix(S1))
  done = []
  for a in A
    -a in done && continue
    iszero(a) && continue
    push!(done,a)
    v = coordinates(a.coeff*basis_matrix(P1), S1)[1,:]
    tmp = [matrix(ZZ,1, degree(S1), (v -  j)*basis_matrix(S1)) for j in Hecke._closest_vectors(S1, v)[2]]
    append!(cvL, tmp)
  end
  if rank(S1) == rank(L)
    #@assert isone(hnf(reduce(vcat, cvL))[1:rank(L),:])
    return cvL
  end
  proj2 = orthogonal_projection(ambient_space(L), basis_matrix(P1))
  L2 = proj2(L)
  proj1 = orthogonal_projection(ambient_space(L), basis_matrix(L2))
  L1 = proj1(L)
  P_Z = ZZ.(solve(basis_matrix(L2),proj2.matrix;side=:left))    
  # recurse 
  for a in characteristic_vectors(L2)
    aL = a*basis_matrix(L2)
    if a*basis_matrix(L2) in L
      @assert rank(L) == ncols(aL)
      push!(cvL, ZZ.(aL))
      continue
    end
    # a vector in L projecting to a
    vL = solve(P_Z, a; side=:left)
    w_amb = vL*proj1.matrix
    w_amb == w_amb*proj1.matrix
    w = coordinates(w_amb[1,:], P1)
    if all(isone, denominator.(w))
      push!(cvL, w*basis_matrix(P1))
      continue 
    end
    _, cv = Hecke._closest_vectors(P1, w)
    tmp = [ZZ.(aL+matrix(QQ,1,length(w),w-j)*basis_matrix(P1)) for j in cv]
    append!(cvL, tmp)
  end
  @assert all(rank(L) == ncols(i) for i in cvL)
  @assert isone(hnf(reduce(vcat, cvL))[1:rank(L),:])
  return cvL
end 


