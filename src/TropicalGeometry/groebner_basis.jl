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
function is_homogeneous(f::MPolyElem)
  leadexpv,tailexpvs = Iterators.peel(exponent_vectors(f))
  d = sum(leadexpv)
  for tailexpv in tailexpvs
    if d!=sum(tailexpv)
      return false
    end
  end
  return true
end
export is_homogeneous

function is_homogeneous(I::MPolyIdeal{K} where {K})
  # todo: test whether generators are interreduced
  @warn "is_homogeneous: merely checking whether given generators are homogeneous, can result in false negative"

  return all(is_homogeneous, gens(I))
end



#=======
tropical Groebner basis
todo: proper documentation
Example:

val_2 = TropicalSemiringMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
w = [0,0,0]
I = ideal([x+2*y,y+2*z])
groebner_basis(I,val_2,w)

Kt,t = RationalFunctionField(QQ,"t")
val_t = TropicalSemiringMap(Kt,t)
Ktx,(x,y,z) = PolynomialRing(Kt,3)
w = [0,0,0]
I = ideal([x+t*y,y+t*z])
groebner_basis(I,val_t,w,return_lead=true)
=======#
@doc Markdown.doc"""
    groebner_basis(I::Ideal, val::TropicalSemiringMap, w::Vector; complete_reduction::Bool, return_lead)

Compute a Groebner basis of `I` over a field with valuation `val` with respect
to weight vector `w`, that is a finite generating set of `I` whose initial
forms generate the initial ideal with respect to `w`.

For the definitions of initial form, initial ideal and Groebner basis see
Section 2.4 of [MS15](@cite).

# Warning
`I` must be homogeneous if `val` is non-trivial or `w` contains non-positive
entries. If `val` is trivla and `w` contains only non-negative entries, then
what is computed is a regular Groebner basis with respect to a weighted
ordering with weight vector `w`.

"""
function groebner_basis(I::MPolyIdeal,val::TropicalSemiringMap,w::Vector{<: Union{Int,Rational{Int}} }; complete_reduction::Bool=false, return_lead::Bool=false)

  ###
  # Step 1: Compute a standard basis in the simulation ring
  ###
  vvI = simulate_valuation(I,val)
  w = simulate_valuation(w,val)
  Rtx = base_ring(vvI)
  # todo: replace with groebner_bases in OSCAR once more orderings are supported
  S,_ = Singular.PolynomialRing(singular_ring(base_ring(Rtx)), map(string, Nemo.symbols(Rtx)), ordering = Singular.ordering_a(w)*Singular.ordering_dp())
  SI = Singular.Ideal(S, [S(g) for g in gens(vvI)])
  vvGB = Singular.gens(Singular.std(SI,complete_reduction=complete_reduction))


  ###
  # Step 2: tighten simulation so that no two monomials of the standard basis
  # elements have the same x-monomial
  ###
  vvGB = [S(tighten_simulation(Rtx(g),val)) for g in vvGB]


  ###
  # Step 3: if complete_reduction = true and val is non-trivial, eliminate
  # tail-monomials contained in the leading ideal in the tropical sense
  #  In the simulation, these monomials corresponds to tail-monomials contained
  #  in the leading ideal up to saturation by t and elimination means
  #  eliminating them after multiplying by a sufficiently high power in t
  ###
  if complete_reduction==true && is_valuation_nontrivial(val)
    sort!(vvGB,lt=x_monomial_lt) # sort vvGB by their leading x monomial from small to large
    Singular.libSingular.set_option("OPT_INFREDTAIL", true)
    for i in 1:length(vvGB)-1
      for j in i+1:length(vvGB)
        vvGB[j] = Singular.reduce(vvGB[j],Singular.std(Singular.Ideal(S,vvGB[i])))
      end
    end
    Singular.libSingular.set_option("OPT_INFREDTAIL", false)
  end

  GB = desimulate_valuation(ideal(Rtx,vvGB),val)
  if return_lead
    vvLI = Singular.lead(vvGB)
    LI = desimulate_valuation(ideal(Rtx,Singular.gens(vvLI)),val)
    return gens(GB),gens(LI)
  end

  return gens(GB)
end


function x_monomial_lt(f::Singular.spoly,g::Singular.spoly)
  expv_f = copy(Singular.leading_exponent_vector(f))
  expv_g = copy(Singular.leading_exponent_vector(g))
  popfirst!(expv_f)
  popfirst!(expv_g)
  return expv_f<expv_g
end
export x_monomial_lt
