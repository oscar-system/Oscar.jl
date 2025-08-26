# The unitary models of complex reflection groups by Lehrer & Taylor (2009).
#
# Ulrich Thiel, 2024


###########################################################################################
# Important remark: In OSCAR, matrices act by default from the right on vectors, so x*A. 
# For example the kernel of a matrix A is the space of all vectors such that x*A=0 
# (confusingly, this is usually referred to a s the *left kernel*....). The same holds for
# Magma and also CHEVIE.
# 
# This means that when we take matrices for the models from the literature in which 
# matrices act from the left (the usual convention for non-computer stuff), like the book
# by Lehrer & Taylor, we need to transpose these matrices! This is why down there in the
# code there are several transpose operation for the final generators.
###########################################################################################

function complex_reflection_group_LT(n::Int)

  if n == 4
    # See Lehrer & Taylor (2009), page 85-86
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is primitive 3rd root of unity
    i = K(i)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r1 = omega//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
    r1prime = omega//2 * matspace(K[-1-i -1+i ; 1+i -1+i])

    push!(gens, transpose(r1))
    push!(gens, transpose(r1prime))

    W = matrix_group(gens)

  elseif n == 5
     # See Lehrer & Taylor (2009), page 85-86
     R,x = polynomial_ring(QQ)
     K,i = number_field(x^2+1, "i")
     R,x = polynomial_ring(K)
     K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is primitive 3rd root of unity
     i = K(i)
     matspace = matrix_space(K, 2, 2)
     gens = elem_type(matspace)[]

     r1 = omega//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
     r2prime = omega//2 * matspace(K[-1+i 1-i ; -1-i -1-i])

     push!(gens, transpose(r1))
     push!(gens, transpose(r2prime))

     W = matrix_group(gens)

  elseif n == 6
    # See Lehrer & Taylor (2009), page 85-86
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is primitive 3rd root of unity
    i = K(i)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r = matspace(K[1 0 ; 0 -1])
    r1 = omega//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
  
    push!(gens, transpose(r))
    push!(gens, transpose(r1))

    W = matrix_group(gens)

  elseif n == 7
    # See Lehrer & Taylor (2009), page 85-86
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is primitive 3rd root of unity
    i = K(i)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r = matspace(K[1 0 ; 0 -1])
    r1 = omega//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
    r2 = omega//2 * matspace(K[-1+i -1+i ; 1+i -1-i])
  
    push!(gens, transpose(r))
    push!(gens, transpose(r1))
    push!(gens, transpose(r2))

    W = matrix_group(gens)

  elseif n == 8
    # See Lehrer & Taylor (2009), page 87-88
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r4 = matspace(K[1 0 ; 0 i])
    r4prime = 1//2 * matspace(K[1+i -1+i ; -1+i 1+i])
  
    push!(gens, transpose(r4))
    push!(gens, transpose(r4prime))

    W = matrix_group(gens)

  elseif n == 9
    # See Lehrer & Taylor (2009), page 87-88
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,sqrt2 = number_field(x^2-2, is_unicode_allowed() ? "√2" : "sqrt2")
    i = K(i)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])
    r4 = matspace(K[1 0 ; 0 i])
  
    push!(gens, transpose(r3))
    push!(gens, transpose(r4))

    W = matrix_group(gens)

  elseif n == 10
    # See Lehrer & Taylor (2009), page 87-88
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is primitive 3rd root of unity
    i = K(i)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r1 = omega//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
    r4prime = 1//2 * matspace(K[1+i -1+i ; -1+i 1+i])
  
    push!(gens, transpose(r1))
    push!(gens, transpose(r4prime))

    W = matrix_group(gens)

  elseif n == 11
    # See Lehrer & Taylor (2009), page 87-88
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is primitive 3rd root of unity
    R,x = polynomial_ring(K)
    K,sqrt2 = number_field(x^2-2, is_unicode_allowed() ? "√2" : "sqrt2")
    i = K(i)
    omega = K(omega)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r1 = omega//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
    r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])
    r4 = matspace(K[1 0 ; 0 i])

    push!(gens, transpose(r1))
    push!(gens, transpose(r3))
    push!(gens, transpose(r4))

    W = matrix_group(gens)

  elseif n == 12
    # See Lehrer & Taylor (2009), page 87-88
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,sqrt2 = number_field(x^2-2, is_unicode_allowed() ? "√2" : "sqrt2")
    i = K(i)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])
    r3prime = 1//sqrt2 * matspace(K[1 1 ; 1 -1])
    r3primeprime = 1//sqrt2 * matspace(K[0 1+i ; 1-i 0])

    push!(gens, transpose(r3))
    push!(gens, transpose(r3prime))
    push!(gens, transpose(r3primeprime))

    W = matrix_group(gens)

  elseif n == 13
    # See Lehrer & Taylor (2009), page 87-88
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,sqrt2 = number_field(x^2-2, is_unicode_allowed() ? "√2" : "sqrt2")
    i = K(i)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r = matspace(K[1 0 ; 0 -1])
    r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])
    r3primeprime = 1//sqrt2 * matspace(K[0 1+i ; 1-i 0])

    push!(gens, transpose(r))
    push!(gens, transpose(r3))
    push!(gens, transpose(r3primeprime))

    W = matrix_group(gens)

  elseif n == 14
    # See Lehrer & Taylor (2009), page 87-88
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is primitive 3rd root of unity
    R,x = polynomial_ring(K)
    K,sqrt2 = number_field(x^2-2, is_unicode_allowed() ? "√2" : "sqrt2")
    i = K(i)
    omega = K(omega)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r1 = omega//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
    r3prime = 1//sqrt2 * matspace(K[1 1 ; 1 -1])

    push!(gens, transpose(r1))
    push!(gens, transpose(r3prime))

    W = matrix_group(gens)

  elseif n == 15
    # See Lehrer & Taylor (2009), page 87-88
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is primitive 3rd root of unity
    R,x = polynomial_ring(K)
    K,sqrt2 = number_field(x^2-2, is_unicode_allowed() ? "√2" : "sqrt2")
    i = K(i)
    omega = K(omega)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r = matspace(K[1 0 ; 0 -1])
    r1 = omega//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
    r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])

    push!(gens, transpose(r))
    push!(gens, transpose(r1))
    push!(gens, transpose(r3))

    W = matrix_group(gens)

  elseif n == 16
    # See Lehrer & Taylor (2009), page 89-90
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,tau = number_field(x^2-x-1, is_unicode_allowed() ? "τ" : "tau") #tau is golden ratio
    i = K(i)
    R,x = polynomial_ring(K)
    K,zeta = number_field(x^2 + (-tau + 1)*x + 1, is_unicode_allowed() ? "ζ" : "zeta") #zeta is prim 5th root of uni
    i = K(i)
    tau = K(tau)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r5 = zeta^2//2 * matspace(K[-tau+i -tau+1 ; tau-1 -tau-i])
    r5prime = -zeta^2//2 * matspace(K[tau-i -tau+1 ; tau-1 tau+i])

    push!(gens, transpose(r5))
    push!(gens, transpose(r5prime))

    W = matrix_group(gens)

  elseif n == 17
    # See Lehrer & Taylor (2009), page 89-90
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,tau = number_field(x^2-x-1, is_unicode_allowed() ? "τ" : "tau") #tau is golden ratio
    i = K(i)
    R,x = polynomial_ring(K)
    K,zeta = number_field(x^2 + (-tau + 1)*x + 1, is_unicode_allowed() ? "ζ" : "zeta") #zeta is prim 5th root of uni
    i = K(i)
    tau = K(tau)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r5 = zeta^2//2 * matspace(K[-tau+i -tau+1 ; tau-1 -tau-i])
    r = matspace(K[1 0 ; 0 -1])

    push!(gens, transpose(r5))
    push!(gens, transpose(r))

    W = matrix_group(gens)

  elseif n == 18
    # See Lehrer & Taylor (2009), page 89-90
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,tau = number_field(x^2-x-1, is_unicode_allowed() ? "τ" : "tau") #tau is golden ratio
    i = K(i)
    R,x = polynomial_ring(K)
    K,zeta = number_field(x^2 + (-tau + 1)*x + 1, is_unicode_allowed() ? "ζ" : "zeta") #zeta is prim 5th root of uni
    i = K(i)
    tau = K(tau)
    R,x = polynomial_ring(K)
    K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is prim 3rd root of unity
    i = K(i)
    tau = K(tau)
    zeta = K(zeta)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r1 = omega//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
    r5 = zeta^2//2 * matspace(K[-tau+i -tau+1 ; tau-1 -tau-i])

    push!(gens, transpose(r1^2))
    push!(gens, transpose(r5))

    W = matrix_group(gens)

  elseif n == 19
    # See Lehrer & Taylor (2009), page 89-90
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,tau = number_field(x^2-x-1, is_unicode_allowed() ? "τ" : "tau") #tau is golden ratio
    i = K(i)
    R,x = polynomial_ring(K)
    K,zeta = number_field(x^2 + (-tau + 1)*x + 1, is_unicode_allowed() ? "ζ" : "zeta") #zeta is prim 5th root of uni
    i = K(i)
    tau = K(tau)
    R,x = polynomial_ring(K)
    K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is prim 3rd root of unity
    i = K(i)
    tau = K(tau)
    zeta = K(zeta)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r = matspace(K[1 0 ; 0 -1])
    r1 = omega//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
    r5 = zeta^2//2 * matspace(K[-tau+i -tau+1 ; tau-1 -tau-i])

    push!(gens, transpose(r))
    push!(gens, transpose(r1))
    push!(gens, transpose(r5))

    W = matrix_group(gens)

  elseif n == 20
    # See Lehrer & Taylor (2009), page 89-90
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,tau = number_field(x^2-x-1, is_unicode_allowed() ? "τ" : "tau") #tau is golden ratio
    i = K(i)
    R,x = polynomial_ring(K)
    K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is prim 3rd root of unity
    i = K(i)
    tau = K(tau)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r1 = omega//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
    r1primeprime = omega//2 * matspace(K[-tau*i-1 (-tau+1)*i ; (-tau+1)*i tau*i-1])

    push!(gens, transpose(r1))
    push!(gens, transpose(r1primeprime))

    W = matrix_group(gens)

  elseif n == 21
    # See Lehrer & Taylor (2009), page 89-90
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,tau = number_field(x^2-x-1, is_unicode_allowed() ? "τ" : "tau") #tau is golden ratio
    i = K(i)
    R,x = polynomial_ring(K)
    K,omega = number_field(x^2+x+1, is_unicode_allowed() ? "ω" : "omega") #omega is prim 3rd root of unity
    i = K(i)
    tau = K(tau)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r = matspace(K[1 0 ; 0 -1])
    r1primeprime = omega//2 * matspace(K[-tau*i-1 (-tau+1)*i ; (-tau+1)*i tau*i-1])

    push!(gens, transpose(r))
    push!(gens, transpose(r1primeprime))

    W = matrix_group(gens)

  elseif n == 22
    # See Lehrer & Taylor (2009), page 89-90
    R,x = polynomial_ring(QQ)
    K,i = number_field(x^2+1, "i")
    R,x = polynomial_ring(K)
    K,tau = number_field(x^2-x-1, is_unicode_allowed() ? "τ" : "tau") #tau is golden ratio
    i = K(i)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    r = matspace(K[1 0 ; 0 -1])
    rprime = matspace(K[0 1 ; 1 0])
    rprimeprime = 1//2 * matspace(K[tau (tau-1)*i+1 ; (-tau+1)*i+1 -tau])

    push!(gens, transpose(r))
    push!(gens, transpose(rprime))
    push!(gens, transpose(rprimeprime))

    W = matrix_group(gens)

  elseif n == 23
    # See Lehrer & Taylor (2009), page 110
    # This group is realized by the line system H3
    R,x = polynomial_ring(QQ)
    K,tau = number_field(x^2-x-1, is_unicode_allowed() ? "τ" : "tau")
    V = vector_space(K,3)
    matspace = matrix_space(K, 3, 3)

    # The line system H3
    cycshift = elem_type(matspace)[]
    push!(cycshift, identity_matrix(K, 3))
    push!(cycshift, matrix(K,3,3,[0 1 0; 0 0 1; 1 0 0])) #(2 3 1)
    push!(cycshift, matrix(K,3,3,[0 0 1; 1 0 0; 0 1 0])) #(3 1 2) 

    lines = [V([1,0,0])*A for A in cycshift]
    lines = vcat(lines, [1//2*V([K(1),tau,tau^-1])*A for A in cycshift])
    lines = vcat(lines, [1//2*V([K(1),tau,-tau^-1])*A for A in cycshift])
    lines = vcat(lines, [1//2*V([K(1),-tau,tau^-1])*A for A in cycshift])
    lines = vcat(lines, [1//2*V([K(1),-tau,-tau^-1])*A for A in cycshift])

    refls = Set([unitary_reflection(l) for l in lines])

    # I guessed the roots of generating reflections from the proof of 
    # Theorem 8.10 (page 144). 
    # The line system H3 is the star-closure of the following lines.
    a = V([0,1,0])
    b = 1//2*V([K(1),tau,tau^-1])
    p = V([1,0,0])  

    # Now, create the corresponding generators
    gens = elem_type(matspace)[]

    push!(gens, matrix(unitary_reflection(a)))
    push!(gens, matrix(unitary_reflection(b)))
    push!(gens, matrix(unitary_reflection(p)))

    W = matrix_group(gens)

    set_attribute!(W, :complex_reflections, refls)

  else

    error("LT model not yet implemented")

  end

  return W

end

function complex_reflection_group_LT(t::Tuple)
  (m,p,n) = t

  # See Lehrer & Taylor (2009), p 35-36
  
  # For the symmetric group case we do not want permutation matrices
  # since this is not an irreducible reflection representation
  if m == 1 && p == 1
    error("LT model not Specified for symmetric group case")
  end

  # See Lehrer & Taylor (2009), p 35-36
  K, z = cyclotomic_field(m)
  matspace = matrix_space(K, n, n)

  #the transpositions
  transp = elem_type(matspace)[]
  for i = 1:n-1
    r = matspace(1)
    r[i,i] = 0
    r[i+1,i+1] = 0
    r[i,i+1] = 1
    r[i+1,i] = 1
    push!(transp, r)
  end

  #the matrix t
  t = matspace(1)
  t[1,1] = z
                          
  #the matrix s
  if n > 1 
    s = t^-1*transp[1]*t
  end

  #now, create the list of generators
  if m == p
    gens = vcat([s], transp)
  elseif p == 1
    gens = vcat([t], transp)
  else
    gens = vcat([s], [t^p], transp)
  end

  W = matrix_group(gens)

  return W

end
