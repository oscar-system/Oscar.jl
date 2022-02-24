###
# Computing links on tropical varieties
# ======================================
#
# Used for the traversal of tropicalizations of algebraic varieties
# For details on the algorithm, see
#   T. Hofmann, Y. Ren: Computing tropical points and tropical links
###



@doc Markdown.doc"""
    homogeneity_space(I::MPolyIdeal; skip_groebner_basis_computation::Bool=false)

Computes the homogeneity space of `I`, i.e., the space of all weight vectors with respect to which `I` is weighted homogeneous. The homogeneity space is the lineality space of the tropicalization independent of the valuation.

Requires a Groebner basis computation, set `skip_groebner_basis_computation=true` at your own risk.

# Examples
```jldoctest
julia> Kp,(p01,p02,p03,p04,p12,p13,p14,p23,p24,p34) = PolynomialRing(QQ,10);

julia> Grass25 = ideal([p03*p12-p02*p13+p01*p23,
                        p04*p12-p02*p14+p01*p24,
                        p04*p13-p03*p14+p01*p34,
                        p04*p23-p03*p24+p02*p34,
                        p03*p12-p02*p13+p01*p23,
                        p04*p12-p02*p14+p01*p24,
                        p04*p13-p03*p14+p01*p34,
                        p14*p23-p13*p24+p12*p34,
                        p03*p12-p02*p13+p01*p23,
                        p04*p12-p02*p14+p01*p24,
                        p04*p23-p03*p24+p02*p34,
                        p14*p23-p13*p24+p12*p34,
                        p03*p12-p02*p13+p01*p23,
                        p04*p13-p03*p14+p01*p34,
                        p04*p23-p03*p24+p02*p34,
                        p14*p23-p13*p24+p12*p34,
                        p04*p12-p02*p14+p01*p24,
                        p04*p13-p03*p14+p01*p34,
                        p04*p23-p03*p24+p02*p34,
                        p14*p23-p13*p24+p12*p34]);

julia> H = homogeneity_space(Grass25) # = linear space in form of a polyhedron

julia> display(rref(affine_equation_matrix(affine_hull(H)))[2]) # = equations in row-echelon form

julia> Kx,(x1,x2) = PolynomialRing(QQ,2);

julia> I = ideal([x1+x2,x1-x2])

julia> dim(homogeneity_space(I))

julia> dim(homogeneity_space(I,skip_groebner_basis_computation=true)) # if GB computation is skipped,
                                                                      # output may not be entire homogeneity space
```
"""
function homogeneity_space(I::MPolyIdeal; skip_groebner_basis_computation::Bool=false)

  if !skip_groebner_basis_computation
    GB = groebner_basis(I,complete_reduction=true)
  else
    GB = gens(I)
  end
  n = length(symbols(base_ring(I)))

  A = zeros(Int,0,n)
  b = zeros(Int,0)
  for g in GB
    leadexpv, tailexpvs = Iterators.peel(exponent_vectors(g))
    for tailexpv in tailexpvs
      A = vcat(A,transpose(tailexpv-leadexpv)) # todo: is there a better way of doing this line?
      push!(b,0)
    end
  end

  return Polyhedron((zeros(Int,0,n),zeros(Int,0)),(A,b))
end
export homogeneity_space


function pivots(M)
  r = 1
  pivotsM = []
  for c = 1:ncols(M)
    if !iszero(M[r,c])
      push!(pivotsM,c)
      r += 1
      if r>nrows(M)
        break
      end
    end
  end
  return pivotsM
end
export pivots



@doc Markdown.doc"""
    tropical_link(I::MPolyIdeal, val::ValuationMap, w::Vector)

Computes the tropical link of `I` around `w`, i.e., returns a minimal array of vectors `L` such that for every maximal Groebner polyhedron $\sigma\in\Trop(I)$ containing `w` there is a $u$ in `L` such that $w+\varepsilon*u\in\sigma$ for $\varepsilon>0$ sufficiently small.

# Warning
Requires that `w` lies in the relative interior of a one-codimensional Groebner polyhedron of `\Trop(I)`.

# Examples
```jldoctest
julia> Kx,(x1,x2,x3,x4) = PolynomialRing(QQ,4);

julia> I = ideal([x1 - 2*x2 + 3*x3, 3*x2 - 4*x3 + 5*x4]);

julia> val_2 = ValuationMap(QQ,2);

julia> w2_1 = [2,0,2,0];

julia> w2_2 = [0,1,0,1];

julia> val_3 = ValuationMap(QQ,3);

julia> w3_1 = [0,0,1,1];

julia> w3_2 = [1,1,0,0];

julia> val_5 = ValuationMap(QQ,5);

julia> w5 = [1,1,1,0];

julia> val_7 = ValuationMap(QQ,7);

julia> w7 = [0,0,0,0];

julia> tropical_link(I,val_2,w2_1)

julia> tropical_link(I,val_2,w2_2)

julia> tropical_link(I,val_3,w3_1)

julia> tropical_link(I,val_3,w3_2)

julia> tropical_link(I,val_4,w4)

julia> tropical_link(I,val_5,w5)

```
"""
function tropical_link(I::MPolyIdeal, val::ValuationMap, w::Vector; local_precision::Integer=19)
  return tropical_link(initial(I,val,w),local_precision=local_precision)
end


#=======
tropical_link
Example:
Kx,(x1,x2,x3,x4,x5) = PolynomialRing(QQ,5)
# inI = ideal([-10*x2^2*x3+9*x1*x3^2+31*x2*x3^2-25*x1*x3*x4-14*x2^2*x5+6*x1*x5^2+46*x2*x5^2-32*x4*x5^2,
#              10*x1^3-2*x1^2*x2+38*x1^2*x3+6*x2^2*x3+34*x1*x2*x4-x3^2*x4+8*x1*x4^2+18*x2*x4^2+47*x3^2*x5-22*x4^2*x5+43*x2*x5^2+16*x3*x5^2,
#              49*x1^2*x4-50*x1*x2*x4-20*x1*x4^2+40*x1*x2*x5-27*x1*x3*x5+4*x1*x4*x5+22*x3*x5^2])
inI = ideal([14*x1*x2-50*x1*x3+x2*x3+13*x3^2+40*x1*x4+16*x1*x5-27*x3*x5,-37*x1*x2+36*x3^2-x2*x4-12*x1*x5+37*x2*x5+12*x3*x5-20*x4*x5-26*x5^2,-2*x2*x3+39*x1*x4-5*x2*x5])
tropical_link(inI,local_precision=29)
=======#
function tropical_link(inI::MPolyIdeal; local_precision::Integer=19)

  ###
  # Step 1: Compute the homogeneity space and identify the pivots (and non-pivots) of its equation matrix in rref
  ###
  H = homogeneity_space(inI)
  _,Eqs = rref(matrix(QQ,lineality_space(H))) # rref returns rank first, which is not required here
  pivotIndices = pivots(Eqs)
  nonpivotIndices = setdiff(collect(1:ncols(Eqs)),pivotIndices)
  if (dim(inI)!=dim(H)+1)
    error("Homogeneity space not one-codimensional.")
  end


  ###
  # Step 2: Construct linear equations and cut down the homogeneity space
  ###
  Kx = base_ring(inI)
  x = gens(Kx)
  inI1 = inI + ideal(Kx,[x[i]-1 for i in pivotIndices]) # dim(inI1) = 1, Trop(inI) = Trop(inI1)+H

  ###
  # Optional: Compute a partially saturated GB using satstd
  ###
  singular_assure(inI1)    # defines necessary objects on the Singular side
  singularIdeal = inI1.gens.S
  singularRing = base_ring(singularIdeal)
  singularIdeal = Singular.satstd(singularIdeal,Singular.MaximalIdeal(singularRing,1))
  inI1 = ideal(Kx,singularIdeal) # cast the Singular ideal back to an Oscar ideal

  ###
  # Step 3.0: Create a t-adic field extension of K
  ###
  K = coefficient_ring(Kx)
  L,t = RationalFunctionField(K,"t")
  val = ValuationMap(L,t)
  Lx,x = PolynomialRing(L,symbols(Kx))
  inI1 = ideal([change_base_ring(L,g) for g in gens(inI1)])

  ###
  # Step 3.1: Intersect the resulting one-dimensional ideal with hyperplanes p*x1-1, ..., p*xn-1, x1+...+xn-p
  ###
  hyperplanes = [x[i]-t for i in nonpivotIndices] # todo: only x[i] for i in indep set
  append!(hyperplanes,[t*x[i]-1 for i in nonpivotIndices])
  rayGenerators = [];
  for hyperplane in hyperplanes
    inI0 = inI1+ideal(Lx,hyperplane)

    ###
    # Optional: Compute a partially saturated GB using satstd
    ###
    singular_assure(inI0)    # defines necessary objects on the Singular side
    singularIdeal = inI0.gens.S
    singularRing = base_ring(singularIdeal)
    Singular.libSingular.set_option("OPT_REDSB", true)
    singularIdeal = Singular.satstd(singularIdeal,Singular.MaximalIdeal(singularRing,1))
    Singular.libSingular.set_option("OPT_REDSB", false)
    inI0 = ideal(Lx,singularIdeal) # cast the Singular ideal back to an Oscar ideal

    if dim(inI0)!=0
      continue
    end

    ###
    # Step 3.2: compute tropical points on slice and merge them to rayGenerators
    ###
    for pointOfSlice in tropical_points(inI0,val,local_precision=local_precision) # todo: check how long it takes to saturate inI0
      #   so that points at infinity are hard erros
      commonDenominator = lcm([denominator(pj) for pj in pointOfSlice])
      pointOfSlice = [commonDenominator*pj for pj in pointOfSlice] # = integer vector
      commonFactor = gcd([numerator(pj) for pj in pointOfSlice])
      rayGenerator = [numerator(pj//commonFactor) for pj in pointOfSlice]
      j = findfirst(isequal(pointOfSlice),rayGenerators)
      if j == nothing
        push!(rayGenerators,pointOfSlice)
      end
    end

  end

  return rayGenerators
end
export tropical_link
