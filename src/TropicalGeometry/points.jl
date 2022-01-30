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
  p = val.uniformizer_field

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
points on the tropical variety
todo: proper documentation
Example:

Kx,(x,y,z) = PolynomialRing(QQ,3)
p = 32003
val_p = ValuationMap(QQ,32003)
I = ideal([x+p*y,y+p*z,x+y+z+1])
p_adic_precision=29
tropical_points(I,val_p,p_adic_precision=29)

# Kt,t = RationalFunctionField(QQ,"t")
# val_t = ValuationMap(Kt,t)
# Ktx,(x,y,z) = PolynomialRing(Kt,3)
# I = ideal([x+t*y,y+t*z])
# tropical_points(I,val_t)
=======#
function tropical_points(I,val_p::ValuationMap{FlintRationalField, fmpq}; p_adic_precision::Int=29) # currently only p-adic supported

  ###
  # Step 0: Check whether I has solutions.
  #   If I isn't zero-dimensional, make it zero-dimensional by adding random affine linear equations.
  ###
  Kx = base_ring(I)
  d = dim(I)
  if (d<0)
    error("input ideal is has no solutions")
  end
  while (d>0)
    # todo: change random_affine_linear_polynomials to be of the form x_i-z_i for i in an independent set and z_i of valuation 0 random
    I = I + ideal(random_affine_linear_polynomials(d,Kx,val_p))
    d = dim(I)
  end

  ###
  # Step 1: Construct a Groebner basis (previously computed for dim)
  #   and pass it to pAdicSolver.solve_affine_groebner_system with a fixed precision.
  #   Increase the precision if necessary.
  ###
  G = groebner_basis(I,complete_reduction=true)

  # while true
  #   try
      Qp = PadicField(val_p.uniformizer_ring,p_adic_precision)
      Gp = [change_base_ring(Qp,g) for g in G]
      global Tp = pAdicSolver.solve_affine_groebner_system(Gp, eigenvector_method=:tropical, ordering=:degrevlex)
  #     break
  #   catch
  #     p_adic_precision *= 2 # double precision if solver unsuccessful
  #   end
  # end

  ###
  # Step 2: remove solutions outside the torus and return
  ###
  T = [];
  for w in eachrow(Tp)
    is_finite = true
    for wj in w
      if max(wj)===Inf
        is_inf = false
        break
      end
    end
    if is_finite
      push!(T,(wi->-min(wi)).(w)) # replace every entry of w with its -min (MAX CONVENTION!!!)
    end
  end

  return T
end
export tropical_points
