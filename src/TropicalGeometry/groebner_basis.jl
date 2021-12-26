###
# Computing (tropical) Groebner bases in Oscar
# ============================================
#
# For a definition of tropical Groebner basis see Section 2.4 in:
#   D. Maclagan, B. Sturmfels: Introduction to tropical geometry
# To see how they can be computed using standard bases see:
#   T. Markwig, Y. Ren: Computing tropical varieties over fields with valuation
###



#=======
return true if f is homogeneous (w.r.t. total degree)
return false otherwise
=======#
function is_homogeneous(f::Union{AbstractAlgebra.Generic.MPoly{K},fmpq_mpoly,fmpz_mpoly} where {K})
  d = sum(exponent_vector(f,1))
  for i in 2:length(f)
    if d!=sum(exponent_vector(f,i))
      return false
    end
  end
  return true
end
export is_homogeneous

function is_homogeneous(I::MPolyIdeal{K} where {K})
  # todo: test whether generators are interreduced
  @warn "is_homogeneous: merely checking whether given generators are homogeneous, can result in false negative"

  for f in gens(I)
    if !is_homogeneous(f)
      return false
    end
  end
  return true
end



#=======
tropical Groebner basis
todo: proper documentation
Example:

val_2 = ValuationMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
w = [0,0,0]
I = ideal([x+2*y,y+2*z])
groebner_basis(I,val_2,w)

Kt,t = RationalFunctionField(QQ,"t")
val_t = ValuationMap(Kt,t)
Ktx,(x,y,z) = PolynomialRing(Kt,3)
w = [0,0,0]
I = ideal([x+t*y,y+t*z])
groebner_basis(I,val_t,w,return_lead=true)
=======#
function groebner_basis(I,val::ValuationMap{valuedField,uniformizer} where{valuedField,uniformizer},w::Vector{Int}; complete_reduction::Bool=false, return_lead::Bool=false)
  vvI = simulate_valuation(I,val)
  w = vcat([-1],w)

  Rtx = base_ring(vvI)
  # todo: replace with groebner_bases in OSCAR once more orderings are supported
  S,_ = Singular.PolynomialRing(singular_ring(base_ring(Rtx)), map(string, Nemo.symbols(Rtx)), ordering = Singular.ordering_a(w)*Singular.ordering_dp())
  SI = Singular.Ideal(S, [S(g) for g in gens(vvI)])

  vvGB = Singular.std(SI,complete_reduction=complete_reduction)
  GB = desimulate_valuation(ideal(Rtx,Singular.gens(vvGB)),val)

  if return_lead
    vvLI = Singular.lead(vvGB)
    LI = desimulate_valuation(ideal(Rtx,Singular.gens(vvLI)),val)
    return gens(GB),gens(LI)
  end

  return gens(GB)
end
