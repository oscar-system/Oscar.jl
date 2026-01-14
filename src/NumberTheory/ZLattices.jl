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
