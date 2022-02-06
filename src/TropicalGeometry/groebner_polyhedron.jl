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
function groebner_polyhedron(I::MPolyIdeal, val::ValuationMap, w::Vector; pertubation::Vector=[], skip_groebner_basis_computation::Bool=false)
  if skip_groebner_basis_computation
    GB = interreduce_tropically(gens(I),val,w,pertubation=pertubation)
  else
    GB = interreduce_groebner_basis(groebner_basis(I,val,w,pertubation=pertubation),val,w,pertubation=pertubation)
  end

  return groebner_polyhedron(GB,initial(GB,val,w,pertubation=pertubation),val)
end

function groebner_polyhedron(GB::Vector{<:MPolyElem}, inGB::Vector{<:MPolyElem}, val::ValuationMap) # GB entries can be MPolyElem and fmpq_mpoly
  eq_lhs = zeros(Int,0,nvars(parent(GB[1])))
  eq_rhs = zeros(Int,0)
  ineq_lhs = zeros(Int,0,nvars(parent(GB[1])))
  ineq_rhs = zeros(Int,0)

  for (f,inf) in zip(GB,inGB)
    ###
    # Step 0: collect the coefficients and exponent vectors of f
    ###
    coefficients_f = collect(coefficients(f))
    exponent_vectors_f = collect(exponent_vectors(f))

    ###
    # Step 1: construct weight equations enforcing that valued weighted degrees of inf are the same
    ###
    inf_leadexpv,inf_tailexpvs = Iterators.peel(exponent_vectors(inf))
    i = findfirst(isequal(inf_leadexpv),exponent_vectors_f)
    if i===nothing
      println(GB)
      println(inGB)
      error("initial forms have monomials which original polynomials do not")
    end
    inf_leadval = val(coefficients_f[i])

    for inf_tailexpv in inf_tailexpvs
      i = findfirst(isequal(inf_tailexpv),exponent_vectors_f)
      if i===nothing
        println(GB)
        println(inGB)
        error("initial forms have monomials which original polynomials do not")
      end
      inf_tailval = val(coefficients_f[i])
      eq_lhs = vcat(eq_lhs,transpose(inf_tailexpv-inf_leadexpv)) # todo: is there a better way of doing this line?
      push!(eq_rhs,inf_tailval-inf_leadval)
    end

    ###
    # Step 2: construct weight inequalities enforcing that valued weighted degree of inf is greater equal f
    ###
    for (f_coeff,f_expv) in zip(coefficients(f),exponent_vectors(f))
      f_val = val(f_coeff)
      ineq_lhs = vcat(ineq_lhs,transpose(f_expv-inf_leadexpv)) # todo: is there a better way of doing this line?
      push!(ineq_rhs,f_val-inf_leadval)
    end
  end

  return Polyhedron((ineq_lhs,ineq_rhs),(eq_lhs,eq_rhs))
end
export groebner_polyhedron
