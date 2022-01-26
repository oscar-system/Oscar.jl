###
# Groebner flip
# =============
###


#=======
groebner flip
Example:
val_2 = ValuationMap(QQ,2)
Kx,(x1,x2,x3,x4) = PolynomialRing(QQ,4)
I = ideal([x1-2*x2+3*x3,3*x2-4*x3+5*x4])
w = [3//2,0,3//2,0]
I = groebner_basis(I,val_2,w)
w = [2,0,2,0]
u = [-1,0,0,0]
groebner_flip(I,val_w,u)
=======#
function groebner_flip(I,val,w,u)
  vvI = simulate_valuation(I,val)
  w = simulate_valuation(w,val)

end
