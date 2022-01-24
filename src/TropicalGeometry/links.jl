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
inI = ideal([21*x1^2+31*x1*x2+46*x3^2-46*x1*x4+17*x2*x4+50*x3*x4-6*x4^2+30*x1*x5-28*x2*x5+46*x3*x5+x4*x5-7*x5^2,39*x1^2-25*x1*x2+37*x2^2-3*x1*x3+15*x2*x3+35*x3^2+44*x1*x5-8*x2*x5+17*x3*x5+37*x4*x5-9*x5^2,20*x1^2+42*x1*x2+8*x1*x3-42*x2*x3-33*x3^2-10*x2*x4+6*x3*x4-3*x1*x5+15*x3*x5-49*x4*x5-26*x5^2])
tropical_link(inI)
=======#
function tropical_link(inI)

  # Compute the homogeneity space and identify the pivots (and non-pivots) of its equation matrix in rref
  H = homogeneity_space(inI,compute_groebner_basis=true)
  _,Eqs = rref(affine_equation_matrix(affine_hull(H)))            # rref returns rank first, which is not required here
  pivotIndices = pivots(Eqs)
  nonpivotIndices = setdiff(collect(1:ncols(Eqs)-1),pivotIndices) # the final column of Eqs represents the RHS of the linear equations,
                                                                  # not a variable in the polynomial ring

  # Construct linear equations to cut down the homogeneity space
  Kx = base_ring(inI)
  x = gens(Kx)
  inI1 = inI + ideal(Kx,[x[i]-1 for i in nonpivotIndices])        # dim(inI1) = 1, Trop(inI) = Trop(inI1)+H

  # Quickly compute a partially saturated GB using satstd (copied from GITFans.jl)
  Oscar.singular_assure(inI1)    # defines necessary objects on the Singular side
  singularIdeal = inI1.gens.S
  singularRing = base_ring(singularIdealSinginI1)
  singularIdeal = Singular.satstd(singularIdeal,Singular.MaximalIdeal(singularRing,1))
  inI1 = ideal(Kx,singularIdeal) # cast the singular ideal back to a oscar ideal

  # Create a p-adic field over a sufficiently large prime
  val_p = ValuationMap(QQ,1000003)   # todo: increase p when necessary

  # Intersect the resulting one-dimensional ideal with hyperplanes p*x1-1, ..., p*xn-1, x1+...+xn-p
  hyperplanes = [val_p.uniformizer*x[i]-1 for i in pivotIndices]
  push!(hyperplanes,sum(x)-val_p.uniformizer)
  rayGenerators = [];
  for hyperplane in hyperplanes
    inI0 = inI1+ideal(Kx,hyperplane)
    inI0 = saturation(inI0,ideal(Kx,prod(x)))
    inI0 = ideal(groebner_basis(inI0))
    print(inI0)
    pointsOnRays = tropical_points(inI0,val_p,p_adic_precision=3) # = matrix of rational numbers
    pointsOnRays = [pointsOnRays[i,:] for i in 1:size(pointsOnRays,1)]                 # = array of vectors of rational numbers
    for pointOnRay in pointsOnRays
      pointOnRay = pointOnRay*lcm([denominator(pi) for pi in pointOnRay])
      pointOnRay = [numerator(pi) for pi in pointOnRay]                                # = vector of integers
      if findfirst(isequal(pointOnRay),rayGenerators) == nothing
        push!(rayGenerators,pointOnRay)
      end
    end
  end

  return rayGenerators
end
export tropical_link
