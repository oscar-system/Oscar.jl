###
# Computing (tropical) Groebner polyhedra in Oscar
# ================================================
#
# For a definition of tropical Groebner cones see Section 2.5 in:
#   D. Maclagan, B. Sturmfels: Introduction to tropical geometry
# To see how they can be computed using standard bases see:
#   T. Markwig, Y. Ren: Computing tropical varieties over fields with valuation
###


#=======
tropical Groebner polyhedron
todo: proper documentation
Example:

val_2 = ValuationMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
w = [0,0,0]
I = ideal([x+2*y,y+2*z])
groebner_polyhedron(I,val_2,w)

Kt,t = RationalFunctionField(QQ,"t")
val_t = ValuationMap(Kt,t)
Ktx,(x,y,z) = PolynomialRing(Kt,3)
w = [0,0,0]
I = ideal([x+t*y,y+t*z])
groebner_polyhedron(I,val_t,w)
=======#
function groebner_polyhedron(I,val::ValuationMap{K,p} where {K,p},w::Vector{Int})
  GB,LI = groebner_basis(I,val,w,complete_reduction=true,return_lead=true)

  A = zeros(Int,0,length(w))
  b = zeros(Int,0)
  for (f,lf) in zip(GB,LI)
    leadcoeff,tailcoeffs = Iterators.peel(coefficients(f))
    leadexpv,tailexpvs = Iterators.peel(exponent_vectors(f))
    leadval = val(leadcoeff)
    for (tailcoeff,tailexpv) in zip(tailcoeffs,tailexpvs)
      tailval = val(tailcoeff)
      A = vcat(A,transpose(tailexpv-leadexpv)) # todo: is there a better way of doing this line?
      push!(b,tailval-leadval)
    end
  end

  return Polyhedron(A,b)
end
export groebner_polyhedron



#=======
homogeneity space
todo: proper documentation
Example:
Kx,(x,y,z) = PolynomialRing(QQ,3)
I = ideal([x+2*y,y+2*z])
homogeneity_space(I)
=======#
function homogeneity_space(I)
  GB = groebner_basis(I,complete_reduction=true)
  n = length(symbols(base_ring(I)))

  A = zeros(Int,0,n)
  b = zeros(Int,0)
  for f in GB
    leadexpv, tailexpvs = Iterators.peel(exponent_vectors(f))
    for tailexpv in tailexpvs
      A = vcat(A,transpose(tailexpv-leadexpv)) # todo: is there a better way of doing this line?
      push!(b,0)
    end
  end

  return Polyhedron((zeros(Int,0,n),zeros(Int,0)),(A,b))
end
export homogeneity_space
