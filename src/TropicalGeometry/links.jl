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
val = ValuationMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
I = ideal([x+2*y,y+2*z])
tropical_link(tropical_link,val)
=======#
function tropical_link(inI,val::ValuationMap{K,p} where{K,p}; coeff_bound::Int=1023, val_bound::Int=15)

  # Compute the homogeneity space and identify the pivots (and non-pivots) of its equation matrix in rref
  H = homogeneity_space(inI,compute_groebner_basis=true)
  r,Eqs = rref(affine_equation_matrix(affine_hull(H)))
  pivotIndices = pivots(Eqs)
  nonpivotIndices = setdiff(collect(1:ncols(Eqs)-1),pivotIndices)
  print(Eqs)

  # Construct linear equations to cut down the homogeneity space
  Kx = base_ring(inI)
  vars = gens(Kx)
  inI0 = inI + ideal(Kx,[vars[i]-1 for i in pivotIndices])

  # Intersect the resulting one-dimensional ideal with hyperplanes
  hyperplanes = [val.uniformizer*vars[i]-1 for i in nonpivotIndices]
  push!(hyperplanes,sum(vars)-val.uniformizer)
  rayGenerators = [];
  for hyperplane in hyperplanes
    pointsOnRays = tropical_points(inI0+ideal(Kx,hyperplane),val)      # = matrix of rational numbers
    pointsOnRays = [pointsOnRays[i,:] for i in 1:size(pointsOnRays,1)] # = array of vectors of rational numbers
    for pointOnRay in pointsOnRays
      pointOnRay = pointOnRay*lcm([denominator(pi) for pi in pointOnRay])
      pointOnRay = [numerator(pi) for pi in pointOnRay]                # = vector of integers
      if findfirst(isequal(pointOnRay),rayGenerators) == nothing
        push!(rayGenerators,pointOnRay)
      end
    end
  end

  return rayGenerators
end
export tropical_link
