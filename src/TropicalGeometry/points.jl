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

Ft,t = RationalFunctionField(FiniteField(32003)[1],"t")
val_t = ValuationMap(Ft,t)
Ftx,(x,y,z) = PolynomialRing(Ft,3)
random_affine_linear_polynomials(3,Ftx,val_t)
=======#
function random_affine_linear_polynomials(k::Int,Kx,val::ValuationMap{K,p} where{K,p}; coeff_bound::Int=1023, val_bound::Int=19)
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


#=======
points on the tropical variety
todo: proper documentation
Example:

Kx,(x,y,z) = PolynomialRing(QQ,3)
p = 32003
val_p = ValuationMap(QQ,32003)
I = ideal([x+p*y,y+p*z,x+y+z+1])
tropical_points(I,val_p,local_precision=29)

Ks,s = RationalFunctionField(QQ,"s")
val_s = ValuationMap(Ks,s)
Ksx,(x,y,z) = PolynomialRing(Ks,3)
I = ideal([x+s*y,y+s*z,x+y+z+1])
tropical_points(I,val_s)


Fs,s = RationalFunctionField(GF(32003),"s")
val_s = ValuationMap(Fs,s)
Fsx,(x,y,z) = PolynomialRing(Fs,3)
I = ideal(Fsx,[x+s*y,y+s*z,x+y+z+1])
tropical_points(I,val_s)
=======#
function tropical_points(I::MPolyIdeal,val::ValuationMap; convention=:max, local_precision::Int=32, primes=[32003,32009,32027,32029,32051,32057,32059,32063])

  ###
  # Step 0: Check whether I has solutions.
  #   If I isn't zero-dimensional, make it zero-dimensional by adding random affine linear equations.
  ###
  Kx = base_ring(I)
  d = dim(I)
  if (d<0)
    println(I)
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
  if is_valuation_p_adic(val)
    ###
    # Case 1: p-adic valuation (handled by pAdicSolver natively)
    ###
    while true
      try
        Qp = PadicField(val.uniformizer_ring,local_precision)
        Gp = [change_base_ring(Qp,g) for g in G]
        global T = pAdicSolver.solve_affine_groebner_system(Gp, eigenvector_method=:tropical, ordering=:degrevlex)
        break
      catch
        println("tropical_point: insufficient precision ",local_precision)
        local_precision *= 2
      end
    end
  elseif is_valuation_t_adic(val) && characteristic(val.residue_field)>0
    ###
    # Case 2: t-adic valuation over finite field (handled by pAdicSolver natively)
    ###
    while true
      try
        Ft,t = LaurentSeriesRing(val.residue_field,local_precision,"t")
        Ftx,x = PolynomialRing(Ft,symbols(base_ring(I)))
        Gt = Vector{elem_type(Ftx)}()
        for g in G
          g = g*lcm([denominator(c) for c in coefficients(g)]) # make all denominators 1
          gFtx = MPolyBuildCtx(Ftx)
          for (cKs,expv) = zip(coefficients(g),exponent_vectors(g))
            cFt = Ft(cKs)
            push_term!(gFtx,cFt,expv)
          end
          push!(Gt,finish(gFtx))
        end

        global T = pAdicSolver.solve_affine_groebner_system(Gt, eigenvector_method=:tropical, ordering=:degrevlex)
        break
      catch
        println("tropical_point: insufficient precision ",local_precision)
        local_precision *= 2
      end
    end
  elseif is_valuation_t_adic(val) && characteristic(val.residue_field)==0
    ###
    # Case 3: t-adic valuation over rationals
    #   (using pAdicSolver over p-adics for p>>0)
    ###
    for p in primes
      try
        K = coefficient_ring(val.valued_ring) # = QQ
        Kx,x = PolynomialRing(K,symbols(base_ring(I)))
        Gp = []
        for g in G
          g = g*lcm([denominator(c) for c in coefficients(g)]) # make all denominators 1
          gKx = MPolyBuildCtx(Kx)
          for (cKs,expvKsx) = zip(coefficients(g),exponent_vectors(g))
            cK = evaluate(numerator(cKs),p)
            push_term!(gKx,cK,expvKsx)
          end
          push!(Gp,finish(gKx))
        end
        Qp = PadicField(p,local_precision)
        Gp = [change_base_ring(Qp,g) for g in Gp]

        # todo: why doesn't the following construction of Gp work for the example above?
        # Qp = PadicField(p,local_precision)
        # Qpx,x = PolynomialRing(Qp,symbols(base_ring(I)))
        # Gp = []
        # for g in G
        #   g = g*lcm([denominator(c) for c in coefficients(g)]) # make all denominators 1
        #   gQpx = MPolyBuildCtx(Qpx)
        #   for (cKs,expv) = zip(coefficients(g),exponent_vectors(g))
        #     cK = evaluate(numerator(cKs),Qp(p))
        #     push_term!(gQpx,cK,expv)
        #   end
        #   push!(Gp,finish(gQpx))
        # end

        global T = pAdicSolver.solve_affine_groebner_system(Gp, eigenvector_method=:tropical, ordering=:degrevlex)
        break
      catch
        println("tropical_point: insufficient precision ",local_precision)
        local_precision *= 2
      end
    end
  else
    error("unsupported valuation or valued field")
  end


  ###
  # Step 2: remove solutions outside the torus and return
  ###
  T_finite = [];
  if convention==:max
    minOrMax = -1
  else
    minOrMax = 1
  end
  for w in eachrow(T)
    is_finite = true
    for wj in w
      if max(wj)===Inf
        is_finite = false
        break
      end
    end
    if is_finite
      push!(T_finite,(wi->minOrMax*min(wi)).(w))
    end
  end

  return T_finite
end
export tropical_points

#===
Example:

K,tK = RationalFunctionField(GF(32003),"t")
L,tL = LaurentSeriesRing(GF(32003),32,"t")
L(tK)
L(tK^32)
===#
function (L::AbstractAlgebra.Generic.LaurentSeriesField)(f::AbstractAlgebra.Generic.Rat)
  @assert isone(denominator(f))
  f = numerator(f)
  tL = gen(L)

  fL = 0
  for (d,c) in enumerate(coefficients(f))
    if iszero(fL)
      fL = c*tL^(d-1) # to make sure that relative precision is optimal
    else
      fL += c*tL^(d-1)
    end
  end

  return fL
end
