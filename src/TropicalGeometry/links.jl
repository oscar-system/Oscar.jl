###
# Computing links on tropical varieties
# ======================================
#
# Used for the traversal of tropicalizations of algebraic varieties
# For details on the algorithm, see
#   T. Hofmann, Y. Ren: Computing tropical points and tropical links
###



#=======
homogeneity space
todo: proper documentation
Example:
Kx,(x,y,z) = PolynomialRing(QQ,3)
I = ideal([x+2*y,y+2*z])
homogeneity_space(I)
=======#
function homogeneity_space(I; compute_groebner_basis::Bool=false, return_linear_equations::Bool=false)

  if compute_groebner_basis
    GB = groebner_basis(I,complete_reduction=true)
  else
    GB = gens(I)
  end
  n = length(symbols(base_ring(I)))

  A = zeros(Int,0,n)
  b = zeros(Int,0)
  for g in GB
    leadexpv, tailexpvs = Iterators.peel(exponent_vectors(g))
    for tailexpv in tailexpvs
      A = vcat(A,transpose(tailexpv-leadexpv)) # todo: is there a better way of doing this line?
      push!(b,0)
    end
  end

  if return_linear_equations
    return Polyhedron((zeros(Int,0,n),zeros(Int,0)),(A,b)),A
  end

  return Polyhedron((zeros(Int,0,n),zeros(Int,0)),(A,b))
end
export homogeneity_space



function tropical_links(inI,val::ValuationMap{K,p} where{K,p}; coeff_bound::Int=1023, val_bound::Int=15)

  H = homogeneity_space(inI,compute_groebner_basis=true) # Polyhedron
  n = ambient_dim(H)
  H = affine_hull(H)

  A = zeros(Int,0,n)
  for h in H
    println(h.a)
  end

end
