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

val_2 = TropicalSemiringMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
random_affine_linear_polynomials(3,Kx,val_2)

Kt,t = RationalFunctionField(QQ,"t")
val_t = TropicalSemiringMap(Kt,t)
Ktx,(x,y,z) = PolynomialRing(Kt,3)
random_affine_linear_polynomials(3,Ktx,val_t)
=======#
function random_affine_linear_polynomials(k::Int,Kx,val::TropicalSemiringMap{K,p} where{K,p}; coeff_bound::Int=1023, val_bound::Int=9)
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



#=======
tropical Groebner basis
todo: proper documentation
Example:

val_2 = TropicalSemiringMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
I = ideal([x+2*y,y+2*z])
tropical_points(I,val_2)

Kt,t = RationalFunctionField(QQ,"t")
val_t = TropicalSemiringMap(Kt,t)
Ktx,(x,y,z) = PolynomialRing(Kt,3)
I = ideal([x+t*y,y+t*z])
tropical_points(I,val_t)
=======#
function tropical_points(I,val::TropicalSemiringMap{K,p} where{K,p})

  Kx = base_ring(I)
  I0 = I + ideal(Kx,random_linear_polynomials(codim(I),Kx,val))

  # todo: compute the tropical variety of I0

  return zeros(Int,0,length(symbols(Kx)))

end
export tropical_points
