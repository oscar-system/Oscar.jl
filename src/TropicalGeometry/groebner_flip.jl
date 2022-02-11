###
# Groebner flip
# =============
###


#=======
groebner flip
Example:
val = ValuationMap(QQ,2)
Kx,(x1,x2,x3,x4) = PolynomialRing(QQ,4)
I = ideal([x1-2*x2+3*x3,3*x2-4*x3+5*x4])
w = [3//2,0,3//2,0]
I = ideal(groebner_basis(I,val,w))
w = [2,0,2,0]
u = [-1,0,0,0]
groebner_flip(I,val,w,u)
=======#
function groebner_flip(I::MPolyIdeal,val::ValuationMap,w::Vector,u::Vector)
  vvI = simulate_valuation(I,val)
  w,u = simulate_valuation(w,u,val)

  Rtx = base_ring(vvI)
  R = singular_ring(coefficient_ring(Rtx))
  tx = map(string, Nemo.symbols(Rtx))

  r,_ = Singular.PolynomialRing(R, tx, ordering = Singular.ordering_a(w)*Singular.ordering_dp())
  rI = Singular.Ideal(r, [r(g) for g in gens(vvI)])
  rI = Singular.std(rI) # todo: this GB should be fast, as rI should already be a GB, check that this is indeed so
  rinI = simulated_initial(rI,w)

  s,_ = Singular.PolynomialRing(R, tx, ordering = Singular.ordering_a(w)*Singular.ordering_a(u)*Singular.ordering_dp())
  sinI = Singular.Ideal(s, [evaluate(g,gens(s)) for g in gens(rinI)])
  sinJ = Singular.std(sinI) # todo: this GB should be fast, as sinI is homogeneous and contains the uniformizer, check that this is indeed so

  println(sinJ)
  rinJ = Singular.Ideal(r, [evaluate(g,gens(r)) for g in gens(sinJ)])
  println("asdf")
  println(rinJ)
  println(rI)
  rJ = Singular.reduce(rinJ,rI)
  println("qwer")

  sJ = Singular.Ideal(s, [evaluate(g,gens(s)) for g in gens(rJ)])
  vvJ = ideal(Rtx,Singular.gens(sJ))

  return desimulate_valuation(vvJ,val)
end
export groebner_flip

function simulated_valued_weighted_degree(f::Singular.spoly, w::Vector; return_vector::Bool=false)
  vwds = [dot(w,alpha) for alpha in exponent_vectors(f)]
  vwd = max(vwds...)
  if return_vector
    return vwd,vwds
  end
  return vwd
end
export simulated_valued_weighted_degree

function simulated_initial(f::Singular.spoly, w::Vector)
  vwd,vwds = simulated_valued_weighted_degree(f, w, return_vector=true)

  initialf = parent(f)(0)
  for (vwdi,term) in zip(vwds,terms(f))
    if vwdi == vwd
      initialf += term
    end
  end

  return initialf
end
function simulated_initial(I::Singular.sideal, w::Vector)
  return Singular.Ideal(base_ring(I),[simulated_initial(f,w) for f in gens(I)])
end
export simulated_initial
