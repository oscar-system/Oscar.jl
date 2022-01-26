###
# Computing multiplicities of maximal-dimensional polyhedra on tropical varieties
# ===============================================================================
#
# Uses the forumla in Lemma 3.4.7 of
#   D. Maclagan, B. Sturmfels: Introduction to tropical geometry
###


function multiplicity(I,val::ValuationMap{K,p} where {K,p},w::Vector{Int})
  return multiplicity(initial(I, val, w))
end


#=======
multiplicity
todo: proper documentation
Example:
Kx,(x1,x2,x3,x4,x5) = PolynomialRing(QQ,5)
# inI = ideal([-10*x2^2*x3+9*x1*x3^2+31*x2*x3^2-25*x1*x3*x4-14*x2^2*x5+6*x1*x5^2+46*x2*x5^2-32*x4*x5^2,
#              10*x1^3-2*x1^2*x2+38*x1^2*x3+6*x2^2*x3+34*x1*x2*x4-x3^2*x4+8*x1*x4^2+18*x2*x4^2+47*x3^2*x5-22*x4^2*x5+43*x2*x5^2+16*x3*x5^2,
#              49*x1^2*x4-50*x1*x2*x4-20*x1*x4^2+40*x1*x2*x5-27*x1*x3*x5+4*x1*x4*x5+22*x3*x5^2])
inI = ideal([23*x1*x3+33*x1^11*x5,-22*x3*x4+44*x2^8*x5,-5*x2^3*x4+33*x5^2])
multiplicity(inI)
=======#
function multiplicity(inI)

  ###
  # Step 1: Compute the homogeneity space and identify the pivots (and non-pivots) of its equation matrix in rref
  ###
  H = homogeneity_space(inI,compute_groebner_basis=true)
  _,Eqs = rref(affine_equation_matrix(affine_hull(H)))            # rref returns rank first, which is not required here
  pivotIndices = pivots(Eqs)
  nonpivotIndices = setdiff(collect(1:ncols(Eqs)-1),pivotIndices) # the final column of Eqs represents the RHS of the linear equations,
                                                                  # not a variable in the polynomial ring

  ###
  # Step 2: Construct linear equations and cut down the homogeneity space
  ###
  Kx = base_ring(inI)
  x = gens(Kx)
  inI0 = inI + ideal(Kx,[x[i]-1 for i in nonpivotIndices])        # dim(inI1) = 1, Trop(inI) = Trop(inI1)+H

  ###
  # Optional: Partially saturate inI0 using satstd
  ###
  singular_assure(inI0)    # defines necessary objects on the Singular side
  singularIdeal = inI0.gens.S
  singularRing = base_ring(singularIdeal)
  singularIdeal = Singular.satstd(singularIdeal,Singular.MaximalIdeal(singularRing,1))

  ###
  # Step 3: Saturate in0 with respect to the maximal ideal
  #   (using saturate as it is faster than repeated satstd,
  #    as ideal is expected to be simple and number of variables is expected to be high)
  ###
  prodx = ideal(Kx,prod(x))
  singular_assure(prodx)
  singularIdeal = Singular.std(Singular.saturation(singularIdeal,prodx.gens.S)[1])
  return Singular.vdim(singularIdeal)
end
