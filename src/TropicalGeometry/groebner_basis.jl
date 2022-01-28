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
@doc Markdown.doc"""
    groebner_basis(I::Ideal, val::ValuationMap, w::Vector; complete_reduction::Bool, return_lead)

Computes a Groebner basis of `I` over a field with valuation `val` with respect to weight vector `w`, that is a finite generating set of `I` whose initial forms generate the initial ideal with respect to `w`.

For the definitions of initial form, initial ideal and Groebner basis see [Maclagan-Sturmfels, Section 2.4].

# Warning
`I` must be homogeneous if `val` is non-trivial or `w` contains non-positive entries. If `val` is trivla and `w` contains only non-negative entries, then what is computed is a regular Groebner basis with respect to a weighted ordering with weight vector `w`.

# Examples
```jldoctest
julia> Kx,(x0,x1,x2,x3,x4,x5) = PolynomialRing(QQ,6);

julia> Cyclic5Homogenized = ideal([x1+x2+x3+x4+x5,
                                   x1*x2+x2*x3+x3*x4+x1*x5+x4*x5,
                                   x1*x2*x3+x2*x3*x4+x1*x2*x5+x1*x4*x5+x3*x4*x5,
                                   x1*x2*x3*x4+x1*x2*x3*x5+x1*x2*x4*x5+x1*x3*x4*x5+x2*x3*x4*x5,
                                   -x0^5+x1*x2*x3*x4*x5]);

julia> Katsura5Homogenized = ideal([-x0+x1+2*x2+2*x3+2*x4+2*x5,
                                    -x0*x1+x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2,
                                    -x0*x2+2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5,
                                    x2^2-x0*x3+2*x1*x3+2*x2*x4+2*x3*x5,
                                    2*x2*x3-x0*x4+2*x1*x4+2*x2*x5]);

julia> val_2 = ValuationMap(QQ,2); # 2-adic valuation

julia> val_3 = ValuationMap(QQ,3); # 3-adic valuation

julia> w = [0,0,0,0,0,0];

julia> groebner_basis(Cyclic5Homogenized, val_2, w, complete_reduction=true)

julia> groebner_basis(Cyclic5Homogenized, val_3, w, complete_reduction=true) # same as for val_2

julia> groebner_basis(Katsura5Homogenized, val_2, w, complete_reduction=true)

julia> groebner_basis(Katsura5Homogenized, val_3, w, complete_reduction=true) # different to val_2

julia> Kt,t = RationalFunctionField(QQ,"t");

julia> Ktx,(x0,x1,x2,x3,x4,x5) = PolynomialRing(Kt,6);

julia> Cyclic5Homogenized_Kt = ideal([change_coefficient_ring(Kt,f) for f in gens(Cyclic5Homogenized)]);

julia> Katsura5Homogenized_Kt = ideal([change_coefficient_ring(Kt,f) for f in gens(Katsura5Homogenized)]);

julia> val_t = ValuationMap(Kt,t); # t-adic valuation

julia> groebner_basis(Cyclic5Homogenized_Kt, val_t, w, complete_reduction=true) # same leading monomials as for val_2 and val_3

julia> groebner_basis(Katsura5Homogenized_Kt, val_t, w, complete_reduction=true) # different leading monomials as for val_2
                                                                                 # same leading monomials as for val_3
```
"""
function groebner_basis(I::MPolyIdeal,val::ValuationMap,w::Vector{<: Union{Int,Rational{Int}} }; complete_reduction::Bool=false, return_lead::Bool=false)

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
  # Step 2: tighten simulation so that no two monomials of the standard basis elements have the same x-monomial
  ###
  vvGB = [S(tighten_simulation(Rtx(g),val)) for g in vvGB]


  ###
  # Step 3: if complete_reduction = true and val is non-trivial,
  #   eliminate tail-monomials contained in the leading ideal in the tropical sense
  #   Inside the tightened simulation, monomials to be eliminated are tail-monomials contained in the leading ideal up to saturation by t
  #   and elimination means eliminating them after multiplying the GB element by a sufficiently high power in t
  ###
  if complete_reduction==true && is_valuation_nontrivial(val)
    sort!(vvGB,lt=x_monomial_lt) # sort vvGB by their leading x monomial from small to large
    Singular.libSingular.set_option("OPT_INFREDTAIL", true)
    for i in 1:length(vvGB)-1
      for j in i+1:length(vvGB)
        t_ecart = x_monomial_ecart(vvGB[j],vvGB[i])
        if t_ecart>=0
          vvGB[j] = Singular.reduce(val.uniformizer_ring^t_ecart*vvGB[j],Singular.std(Singular.Ideal(S,vvGB[i])))
          vvGB[j] = S(tighten_simulation(Rtx(vvGB[j]),val))
        end
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


###
# returns true if (leading x-monomial of f) <_lex (leading x-monomial of g)
# returns false otherwise
###
function x_monomial_lt(f::Singular.spoly, g::Singular.spoly)
  exp_x_f = copy(Singular.leading_exponent_vector(f))
  exp_x_g = copy(Singular.leading_exponent_vector(g))
  popfirst!(exp_x_f)
  popfirst!(exp_x_g)
  return exp_x_f<exp_x_g
end
export x_monomial_lt


###
# returns true if x^expv_g divides x^expv_f
# returns false otherwise
###
function x_monomial_divides(exp_x_g::Vector,exp_x_f::Vector)
  for (eg,ef) in zip(exp_x_g,exp_x_f)
    if eg>ef
      return false
    end
  end
  return true
end
export x_monomial_divides


###
# if the leading x-monomial of g divides x-monomials in f
#   returns l=max(0, (t_g-exponent of g) - (t-exponents of f))
#   so that f*t^l can be reduced by g to eliminate all the x-monomials
# otherwise, returns -1
###
function x_monomial_ecart(f::Singular.spoly, g::Singular.spoly)
  exp_x_g = copy(Singular.leading_exponent_vector(g))
  exp_t_g = popfirst!(exp_x_g)
  e = 0
  dividend_found = false
  for exp_f in exponent_vectors(f)
    exp_x_f = copy(exp_f)
    exp_t_f = popfirst!(exp_x_f)
    if x_monomial_divides(exp_x_g,exp_x_f)
      e = max(e,exp_t_g-exp_t_f)
      dividend_found = true
    end
  end
  if dividend_found
    return e
  end
  return -1
end
export x_monomial_ecart
