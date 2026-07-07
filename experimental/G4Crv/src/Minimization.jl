


function minimize_and_reduce_cubic_surface(cubic::MPolyRingElem{T}) where T <:Union{QQFieldElem, ZZRingElem}
  
  mu = lcm(map(denominator, collect(coefficients(cubic))))
  R = parent(cubic)
  R, _ = grade(R)
  cubic *= mu
  cubic = change_coefficient_ring(ZZ, cubic)
  cubic /= content(cubic)

  bad_primes = find_bad_primes(cubic)

  M = identity_matrix(ZZ,4)
  for p in bad_primes
    cubic, M1 = minimize_cubic_at_p(cubic, p)
    #TODO: Compose M and M1
  end
  return cubic
end

function find_bad_primes(cubic::MPolyRingElem{T}) where T
  R_aff, X = polynomial_ring(ZZ, [:X1, :X2, :X3])

  #Find the bad primes of the surface
  bad_integer = ZZ(1)
  for i in (1:4)
    X = gens(R_aff)
    insert!(X, i, one(R_aff))
    cubic_aff = cubic(X...)
    I = [derivative(cubic_aff, i) for i in (1:3)]
    push!(I, cubic_aff)
    push!(I, 6*det(hessian_matrix(cubic_aff)))

    I = ideal(R_aff, I)
    B = groebner_basis(I, ordering = monomial_ordering(R_aff, :degrevlex))
    D = B[1]
    if total_degree(D) == 0 
      bad_integer = lcm(coeff(D, 1), bad_integer)
    else
      return cubic
    end
  end

  bad_primes = prime_divisors(bad_integer)
  return bad_primes
end

function minimize_cubic_at_p(cubic, p)
  M = identity_matrix(ZZ, 4)
  while true
    test, cub, M1 = try_minimize_cubic_at_p(cubic, p)
    #TODO: Compose M and M1
    if !test 
      return cubic, M
    else
      cubic = cub
    end
  end
end

function try_minimize_cubic_at_p(cubic, p)
  F = GF(p)
  cubic_p = change_coefficient_ring(F, cubic)
  R_p = parent(cubic_p)
  R_p, _ = grade(R_p)
  cubic_p = R_p(cubic_p)

  success, M = weight0001(cubic_p)

  if success
    cubic_base_change = cubic_new_basis(cubic, transpose(M))
    c = content(cubic_base_change)
    @req divides(c, p)[1] "There is an error in minimization of cubic surfaces."
    cubic = cubic_base_change/c
    return true, cubic, transpose(M)
  end

  success, M = weight0011(cubic_p)
  if success
    cubic_base_change = cubic_new_basis(cubic, transpose(M))
    c = content(cubic_base_change)
    @req divides(c, p)[1] "There is an error in minimization of cubic surfaces."
    if divides(c, p^2)[1]
      cubic = cubic_base_change/c
      return true, cubic, transpose(M)
    end
  end

  L1 = compute_relevant_singular_points(cubic, p)
  println("Trying for singular points")
  for P in L1
    K = transpose(matrix(P))
    M = map(x-> lift(ZZ, x), K)
    M = [M ; p*identity_matrix(ZZ, 4)]
    hnf!(M)
    M = M[1:4,:]
    cubic_base_change = cubic_new_basis(cubic, transpose(M))
    c = content(cubic_base_change)
    if divides(c, p^3)[1]
      cubic = cubic_base_change/c
      return true, cubic, transpose(M)
    end

    if divides(c, p^2)[1]
      cubic_temp = cubic_base_change/c
      success, cub, M = weight0122_and0223(cubic_temp,p)
      if success
        return cub
      end
    end

  end
  return false, cubic, transpose(M)

end

function _factor(f)
  Fx = parent(f)
  F = base_ring(f)
  FF = Nemo.Native.GF(characteristic(base_ring(f)))
  FFx, = polynomial_ring(FF, nvars(parent(f)))
  ff = map_coefficients(c -> FF(lift(ZZ, c)), f; parent = FFx)
  fffac = factor(ff)
  return Fac(Fx(lift(ZZ, constant_coefficient(unit(fffac)))), Dict(map_coefficients(c -> F(lift(ZZ, c)), g; parent = Fx) => e for (g, e) in fffac))
end

  #F = Nemo.Native.GF(ZZ(p))
function weight0001(cubic_p)
  println("try weight0001")
  R_p = parent(cubic_p)
  F = base_ring(R_p)
  p = characteristic(F)
  success = false
  facs = _factor(cubic_p)
  M = zero_matrix(F, 4, 4)
  for (g, e) in facs
    if total_degree(g) == 1
      success = true
      K = kernel(matrix([coeff(g, gens(R_p)[i]) for i in (1:4)]))
      M = map(x-> lift(ZZ, x), K)
      M = [M ; p*identity_matrix(ZZ, 4)]
      hnf!(M)
      M = M[1:4,:]
      break
    end
  end
  return success, M
end

function weight0011(cubic_p)
  println("Try weight 0011")
  R_p = parent(cubic_p)
  F = base_ring(R_p)
  p = characteristic(F)
   #Case: weight (0,0,1,1)
  M = zero_matrix(ZZ, 4, 4)
  cubic_sch = projective_scheme(ideal(R_p(cubic_p)))
  S = singular_locus(cubic_sch)
  decomposition = irreducible_components(S)
  success = false
  for l in decomposition
    l_red = reduced_scheme(l) 
    if degree(l_red) == 1 && dim(l_red) == 1 #We have found a line
      success = true
      points = Vector{FqFieldElem}[]
      while length(points) < 2
        v = rand(F, 4)
        if v == zeros(F,4)
          continue
        end
        if is_on_scheme(v, l_red)
          if (length(points) == 0) || rank(matrix([points[1], v])) == 2
            push!(points, v)
          end
        end
      end
      K = matrix(points)
      M = map(x-> lift(ZZ, x), K)
      M = [M ; p*identity_matrix(ZZ, 4)]
      hnf!(M)
      M = M[1:4,:]
      break
    end
  end
  return success, M
end

function weight0122_and0223(cubic_temp, p)
  println("0122 and 0223")
  F = GF(p)
  cubic_p_temp = change_base_ring(F, cubic_temp)
  @req !is_irreducible(cubic_p_temp) "This can't happen. There is a bug in the minimization algorithm."
  
  M = zero_matrix(ZZ, 4, 4)
  R_p = parent(cubic_p_temp)
  R_p, _ = grade(R_p)
  cubic_p_temp = R_p(cubic_p_temp)
  x, y, z, w, = gens(R_p)
  if valuation(cubic_p_temp, x) != 1
    return false, cubic_temp, transpose(M)
  end

  q_factor = cubic_p_temp/x
  decomposition = [f for (f,e) in factor(q_factor)]
  if length(decomposition) == 2 || length(decomposition) == 1 && degree(decomposition[1]) == 2
  #Quadratic case
    F = base_ring(R_p)
    q_sch = projective_scheme(ideal(R_p(q_factor)))
    S = singular_locus(q_sch)
    S_dec = irreducible_components(S)
    #Exactly one line is found.
    if length(S_dec) == 1 && degree(S_dec[1]) == 1 
      l = S_dec[1]
      points = Vector{FqFieldElem}[]
      while length(points) < 2
        v = rand(F, 4)
        if is_on_scheme(v, l)
          if (length(points) == 0 || v!= points[1]) && v != zeros(F, 4)
            push!(points, v)
          end
        end
      end
      K = matrix(points)
      M = map(x-> lift(ZZ, x), K)
      M = [M ; p*identity_matrix(ZZ, 4)]
      hnf!(M)
      M = M[1:4,:]
      cubic_base_change = cubic_new_basis(cubic_temp, transpose(M))
      c = content(cubic_base_change)
      if divides(c, p^2)[1]
        cubic_temp = cubic_base_change/c
        return true, cubic_base_change, transpose(M)
      else
        return false, cubic_temp, transpose(M)
      end
    else
      #Now we are in the case where q_factor is a square
      g = sqrt(q_factor)
      #Copied from weight0001
      K = kernel(matrix([coeff(g, gens(R_p)[i]) for i in (1:4)]))
      M = map(x-> lift(ZZ, x), K)
      M = [M ; p*identity_matrix(ZZ, 4)]
      hnf!(M)
      M = M[1:4,:]
      cubic_base_change = cubic_new_basis(cubic_temp, transpose(M))
      c = content(cubic_base_change)
      if divides(c, p^2)[1]
        cubic_temp = cubic_base_change/c
        return true, cubic_base_change, transpose(M)
      else
        cubic_temp = cubic_base_change/p 
        cubic_temp_p = change_coefficient_ring(F, cubic_temp)
        success, M = weight0001(cubic_p)

        if success
          cubic_base_change = cubic_new_basis(cubic, transpose(M))
          c = content(cubic_base_change)
          @req divides(c, p)[1] "There is an error in minimization of cubic surfaces."
          cubic = cubic_base_change/c
          return true, cubic_base_change, transpose(M)
        end
      end
    end
    #0223
    cub_sch = projective_scheme(ideal(R_p(cubic_temp)))
    S = singular_locus(q_sch)
    S_dec = irreducible_components(S)
    #Exactly one line is found.
    for l in S_dec
      l_red = reduced_scheme(l)
      if degree(l_red) == 1 && dim(l_red) == 1 
        points = Vector{FqFieldElem}[]
        while length(points) < 2
          v = rand(F, 4)
          if is_on_scheme(v, l_red)
            if (length(points) == 0 || v!= points[1]) && v != zeros(F, 4)
              push!(points, v)
            end
          end
        end
        K = matrix(points)
        M = map(x-> lift(ZZ, x), K)
        M = [M ; p*identity_matrix(ZZ, 4)]
        hnf!(M)
        M = M[1:4,:]
        cubic_base_change = cubic_new_basis(cubic, transpose(M))
        c = content(cubic_base_change)
        @req divides(c, p)[1] "There is an error in minimization of cubic surfaces."
        g = cubic_base_change/c
        if divides(c, p^2)
          return true, g, M 
        end
        g_p = change_base_ring(F, g)
        facs = factor(g_p)
        for (g, e) in facs
          if total_degree(g) == 1
            K = kernel(matrix([coeff(g, gens(R_p)[i]) for i in (1:4)]))
            K = Hecke.complete_to_basis(K)
            M = map(x-> lift(ZZ, x), K)
            M = [M ; p*identity_matrix(ZZ, 4)]
            hnf!(M)
            M = M[1:4,:]

            cubic_base_change = cubic_new_basis(g, transpose(M))
            c = content(cubic_base_change)
            @req divides(c, p)[1] "There is an error in minimization of cubic surfaces."
            if divides(c, p^2)
              return true, cubic_base_change/c, M 
            end
          end
        end
      end
    end
    return false, cubic_temp, transpose(M)
  end
  return false, cubic_temp, transpose(M)
end

function compute_relevant_singular_points(cubic, p)
  F = GF(p)
  cubic_p = change_base_ring(F, cubic)
  R_p = parent(cubic_p)
  R_p, _ = grade(R_p)
  cubic_p = cubic_p
  p = characteristic(F)
   #Case: weight (0,0,1,1)
  M = zero_matrix(ZZ, 4, 4)
  cubic_sch = projective_scheme(ideal(R_p(cubic_p)))
  S = singular_locus(cubic_sch)
  decomposition = irreducible_components(S)
  L1 = Vector{FqFieldElem}[]
  for l in decomposition
    l_red = reduced_scheme(l) 
    if dim(l_red) == 0
      J = defining_ideal(l_red)
      P = kernel(transpose(matrix([[coeff(gens(J)[j], gens(R_p)[i]) for i in (1:4)] for j in (1:3)])))
      push!(L1, P[1,:])
    end
    if degree(l_red) == 1 && dim(l_red) == 1 #We have found a line
      param = _parametrization(l_red)
      if p <=3
        V = Iterators.product(0:p-1, 0:p-1)
        for (s,t) in V
          L1 = push!(L1, map(p -> p(s,t) , param))
        end
        L1 = L1[2:end]
      else
        param_ZZ = _ZZ_parametrization(l_red)
        f_p_st = change_base_ring(F, cubic(param_ZZ...)/p)
        Rt, t = polynomial_ring(F, :t)
        for u in F 
          eq = f_p_st(u, t)
          if eq == zero(Rt)
            push!(L1, map(p -> p(u, F(1)) , param_ZZ))
            if u != zero(GF(p))
              push!(L1, map(p -> p(u, F(0)) , param_ZZ))
            end
          else
            for v in roots(eq)
              push!(L1, map(p -> p(u,v) , param_ZZ))
            end
          end
        end
      end
    end

    if degree(l_red) == 2
      println("Test this case?")
    end

    #Check if this case even makes sense
    if degree(l_red) == 3
      println("Test this case?")
    end
  end
  return L1
end

function _ZZ_parametrization(line)
  I = defining_ideal(line)
  R_p = base_ring(I)
  S, (s, t) = polynomial_ring(ZZ,  [:s,:t])
  par = kernel(map(x-> lift(ZZ, x), transpose(matrix([[coeff(I[j], gens(R_p)[i]) for i in (1:4)] for j in [1,2]]))))
  return [s*par[1,i] + t*par[2,i] for i in (1:4)]
end

function _parametrization(line)
  I = defining_ideal(line)
  R_p = base_ring(I)
  F = base_ring(R_p)
  S, (s, t) = polynomial_ring(F,  [:s,:t])
  par = kernel(transpose(matrix([[coeff(I[j], gens(R_p)[i]) for i in (1:4)] for j in [1,2]])))
  return [s*par[1,i] + t*par[2,i] for i in (1:4)]
end



function is_on_scheme(v::Vector, X::ProjectiveScheme)
  return all([iszero(evaluate(f, v)) for f in gens(defining_ideal(X))])
end
