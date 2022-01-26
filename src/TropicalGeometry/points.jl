###
# Computing points on tropical varieties
# ======================================
#
# Used for computing starting points for the traversal of tropical varieties
# For details on the algorithm, see
#   T. Hofmann, Y. Ren: Computing tropical points and tropical links
###



#=======
random linear polynomials where coefficients have uniform valuation
todo: proper documentation

val_2 = ValuationMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
random_affine_linear_polynomials(3,Kx,val_2)

Kt,t = RationalFunctionField(QQ,"t")
val_t = ValuationMap(Kt,t)
Ktx,(x,y,z) = PolynomialRing(Kt,3)
random_affine_linear_polynomials(3,Ktx,val_t)
=======#
function random_affine_linear_polynomials(k::Int,Kx,val::ValuationMap{K,p} where{K,p}; coeff_bound::Int=1023, val_bound::Int=9)
  n = length(gens(Kx))
  p = val.uniformizer

  coeffs = rand(0:coeff_bound,k,n+1)
  vals = rand(0:val_bound,k,n+1)
  expvs = zeros(Int,n+1,n) # Question: is there a simpler way to construct identity matrix?
  for i in 1:n             #   identity_matrix(...), diagonal_matrix(ones(...)) do not seem to work
    expvs[i,i] = 1
  end

  lin_polys = []
  for i in 1:k
    lin_poly = MPolyBuildCtx(Kx)
    for (c,v,j) in zip(coeffs[i,:],vals[i,:],1:n+1)
      push_term!(lin_poly,c*p^v,expvs[j,:])
    end
    push!(lin_polys,finish(lin_poly))
  end

  return lin_polys
end
export random_affine_linear_polynomials


function contains_zero_entry(p)
  for pi in p
    if pi==0
      return true
    end
  end
  return false
end


#=======
tropical Groebner basis
todo: proper documentation
Example:

val_2 = ValuationMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
I = ideal([x+2*y,y+2*z,x+y+z+1])
tropical_points(I,val_2)

# Kt,t = RationalFunctionField(QQ,"t")
# val_t = ValuationMap(Kt,t)
# Ktx,(x,y,z) = PolynomialRing(Kt,3)
# I = ideal([x+t*y,y+t*z])
# tropical_points(I,val_t)
=======#
function tropical_points(I,val_p::ValuationMap{FlintRationalField, fmpz}; p_adic_precision::Int=9, remove_points_at_infinity::Bool=false) # currently only p-adic supported

  Kx = base_ring(I)
  K = base_ring(Kx)
  I0 = I

  while (dim(I0)!=0)
    I0 = I + ideal(Kx,random_affine_linear_polynomials(dim(I),Kx,val_p))
  end

  while true
    try
      Qp = PadicField(val_p.uniformizer,p_adic_precision)
      I0p = [change_base_ring(Qp,f) for f in groebner_basis(I0)]
      global TI0p = solve_macaulay(I0p,eigenvector_method = "tropical") #todo: pass groebner basis flag+ordering
      break
    catch
      p_adic_precision *= 2 # double precision if solver unsuccessful
    end
  end

  if remove_points_at_infinity
    return [TI0p[i,:] for i in 1:size(TI0p,1) if !contains_approximate_zero_entry(TI0p[i,:],-2^63)] # pAdicSolver returns valuation -2^63 if it encounters 0
  else
    return [TI0p[i,:] for i in 1:size(TI0p,1)]
  end
end
export tropical_points


function contains_approximate_zero_entry(tropical_point,p_adic_precision)
  for tropical_point_entry in tropical_point
    if tropical_point_entry>=p_adic_precision
      return true
    end
  end
  return false
end
