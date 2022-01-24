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
function homogeneity_space(I; compute_groebner_basis::Bool=false)

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

  return Polyhedron((zeros(Int,0,n),zeros(Int,0)),(A,b))
end
export homogeneity_space


function pivots(M)
  n = nrows(M)
  m = ncols(M)
  c = 1
  pivotsM = []
  for j = 1:m-1
    for i = c:min(j,n)
      if !iszero(M[i,j])
        push!(pivotsM,j)
        c += 1
        break
      end
    end
  end
  return pivotsM
end
export pivots


function rational_matrix_clear_denom(M)
end


#=======
tropical_link
todo: proper documentation
Example:
Kx,(x1,x2,x3,x4,x5) = PolynomialRing(QQ,5)
# inI = ideal([-10*x2^2*x3+9*x1*x3^2+31*x2*x3^2-25*x1*x3*x4-14*x2^2*x5+6*x1*x5^2+46*x2*x5^2-32*x4*x5^2,
#              10*x1^3-2*x1^2*x2+38*x1^2*x3+6*x2^2*x3+34*x1*x2*x4-x3^2*x4+8*x1*x4^2+18*x2*x4^2+47*x3^2*x5-22*x4^2*x5+43*x2*x5^2+16*x3*x5^2,
#              49*x1^2*x4-50*x1*x2*x4-20*x1*x4^2+40*x1*x2*x5-27*x1*x3*x5+4*x1*x4*x5+22*x3*x5^2])
inI = ideal([14*x1*x2-50*x1*x3+x2*x3+13*x3^2+40*x1*x4+16*x1*x5-27*x3*x5,-37*x1*x2+36*x3^2-x2*x4-12*x1*x5+37*x2*x5+12*x3*x5-20*x4*x5-26*x5^2,-2*x2*x3+39*x1*x4-5*x2*x5])
tropical_link(inI)
=======#
function tropical_link(inI; p_adic_prime=1000003)

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
  inI1 = inI + ideal(Kx,[x[i]-1 for i in nonpivotIndices])        # dim(inI1) = 1, Trop(inI) = Trop(inI1)+H

  ###
  # Optional: Compute a partially saturated GB using satstd
  ###
  singular_assure(inI1)    # defines necessary objects on the Singular side
  singularIdeal = inI1.gens.S
  singularRing = base_ring(singularIdeal)
  singularIdeal = Singular.satstd(singularIdeal,Singular.MaximalIdeal(singularRing,1))
  inI1 = ideal(Kx,singularIdeal) # cast the Singular ideal back to an Oscar ideal

  ###
  # Step 3.0: Create a p-adic field over a sufficiently large prime
  ###
  val_p = ValuationMap(QQ,p_adic_prime)   # todo: increase p when necessary

  ###
  # Step 3.1: Intersect the resulting one-dimensional ideal with hyperplanes p*x1-1, ..., p*xn-1, x1+...+xn-p
  ###
  hyperplanes = [val_p.uniformizer*x[i]-1 for i in pivotIndices]
  push!(hyperplanes,sum(x)-val_p.uniformizer)
  rayGenerators = [];
  # rayMultiplicities = []; # ray multiplicities cannot be generally computed using this method,
                            # however this method gives lower bounds on the multiplicities which may be used for sanity checking later
  # println("==================================")
  # println("inI1",inI1)
  # println("==================================")
  for hyperplane in hyperplanes
    inI0 = inI1+ideal(Kx,hyperplane)

    ###
    # Optional: Compute a partially saturated GB using satstd
    ###
    singular_assure(inI0)    # defines necessary objects on the Singular side
    singularIdeal = inI0.gens.S
    singularRing = base_ring(singularIdeal)
    singularIdeal = Singular.satstd(singularIdeal,Singular.MaximalIdeal(singularRing,1))
    inI0 = ideal(Kx,singularIdeal) # cast the Singular ideal back to an Oscar ideal

    ###
    # Step 3.2: bookkeeping points on slice and their multiplicities
    ###
    pointsOfSliceMatrix = tropical_points(inI0,val_p) # = rational matrix
    pointsOfSlice = []
    multsOfSlice = []
    for i in 1:size(pointsOfSliceMatrix,1)
      pointOfSlice = pointsOfSliceMatrix[i,:]        # = rational vector
      commonDenominator = lcm([denominator(pj) for pj in pointOfSlice])
      pointOfSlice = [numerator(commonDenominator*pj) for pj in pointOfSlice] # = integer vector
      j = findfirst(isequal(pointOfSlice),pointsOfSlice)
      if j == nothing
        push!(pointsOfSlice,pointOfSlice)
        push!(multsOfSlice,1)
      else
        multsOfSlice[j] += 1
      end
    end

    # println("==================================")
    # println("slice ",hyperplane)
    # println("points: ",pointsOfSlice)
    # println("mults: ",multsOfSlice)
    # println("==================================")

    ###
    # Step 3.3: merge points and multiplicities on slice and check for consistency
    ###
    for (pointOfSlice,m) in zip(pointsOfSlice,multsOfSlice)
      j = findfirst(isequal(pointOfSlice),rayGenerators)
      if j == nothing
        push!(rayGenerators,pointOfSlice)
        # push!(rayMultiplicities,m)
      end
    end
  end

  return rayGenerators,rayMultiplicities
end
export tropical_link
