# The unitary models of complex reflection groups by Lehrer & Taylor (2009).
#
# Ulrich Thiel, 2024


###########################################################################################
# Important remark: In OSCAR, matrices act by default from the right on coordinate vectors.
# For example, the kernel of a matrix A is the space of all vectors x such that x * A = 0.
# (confusingly, this is usually referred to a s the *left kernel*....). The same holds for
# Magma and also CHEVIE.
#
# This means that when we take matrices for the models from the literature in which
# matrices act from the left (the usual convention for non-computer stuff), like the book
# by Lehrer & Taylor, we need to transpose these matrices! This is why down there in the
# code there are several transpose operation for the final generators.
###########################################################################################

# Construct reflections directly in OSCAR's right-action convention. For a row root a,
# the matrix is I - (1-zeta)/(a*conj(a)^t) * conj(a)^t*a; scaling a changes nothing.
# See Lehrer--Taylor, Definition 1.17 and equation (1.18) (p. 11), and
# Definition 8.18 (p. 147).
function _lt_root_group(K, roots, eigenvalues=fill(K(-1), length(roots)))
  V = vector_space(K, length(first(roots)))
  return matrix_group([
    matrix(unitary_reflection(V(root), K(zeta)))
    for (root, zeta) in zip(roots, eigenvalues)
  ])
end

function _lt_coordinate_root(K, d::Int, i::Int)
  root = fill(zero(K), d)
  root[i] = one(K)
  return root
end

function _lt_adjacent_roots(K, n::Int, d::Int=n)
  return [
    _lt_coordinate_root(K, d, i) - _lt_coordinate_root(K, d, i + 1)
    for i in 1:n-1
  ]
end

function _lt_d_roots(K, n::Int, d::Int=n)
  roots = _lt_adjacent_roots(K, n, d)
  push!(roots, _lt_coordinate_root(K, d, n - 1) + _lt_coordinate_root(K, d, n))
  return roots
end

# All coordinates in the constructions below refer to an orthonormal basis.
# The identifications with ST numbers are Lehrer--Taylor, Table D.2 (p. 273).
function _complex_reflection_group_LT_higher_rank(n::Int)
  if n == 24
    # J_3^(4): B_3 and the orbit of (lambda, lambda, 2), §7.6.2, pp. 108--109.
    R, x = polynomial_ring(QQ)
    K, lambda = number_field(x^2 + x + 2, "lambda")
    roots = _lt_adjacent_roots(K, 3)
    push!(roots, K.([0, 0, 1]), K.([lambda, lambda, 2]))
    return _lt_root_group(K, roots)
  elseif n in (25, 26, 32)
    # L_3, M_3, L_4: §8.5.3--5.5, p. 149. The order-three reflections have
    # eigenvalue omega as in Definition 8.18. Lemmas 8.23--8.24 (p. 150) give
    # L_3 from the first three roots and L_4 by adjoining (0,1,-1,-1).
    # Adjoining a transposition to L_3 gives M_3.
    R, x = polynomial_ring(QQ)
    K, omega = number_field(x^2 + x + 1, "omega")
    d = n == 32 ? 4 : 3
    roots = [
      _lt_coordinate_root(K, d, 1),
      _lt_coordinate_root(K, d, 2),
      vcat(K.([1, 1, 1]), fill(zero(K), d - 3)),
    ]
    eigenvalues = fill(omega, 3)
    if n == 26
      # Keep the order-two hyperplane orbit first, followed by the ternary orbit.
      pushfirst!(roots, K.([1, -1, 0]))
      pushfirst!(eigenvalues, K(-1))
    elseif n == 32
      push!(roots, K.([0, 1, -1, -1]))
      push!(eigenvalues, omega)
    end
    return _lt_root_group(K, roots, eigenvalues)
  elseif n in (27, 30)
    # J_3^(5) and H_4: §7.6.3, pp. 110--111. The three H_3 roots are those
    # of Theorem 8.10 (p. 144); the H_4 extension is Theorem 8.15 (p. 146).
    R, x = polynomial_ring(QQ)
    K, tau = number_field(x^2 - x - 1, "tau")
    if n == 27
      R, x = polynomial_ring(K)
      K, omega = number_field(x^2 + x + 1, "omega")
      tau = K(tau)
      roots = [
        K.([0, 1, 0]), K.([1, tau, inv(tau)]), K.([1, 0, 0]),
        K.([tau + omega, inv(tau) * omega - 1, 0]),
      ]
    else
      roots = [
        K.([0, 1, 0, 0]), K.([1, tau, inv(tau), 0]), K.([1, 0, 0, 0]),
        K.([1, 1, 1, 1]),
      ]
    end
    return _lt_root_group(K, roots)
  elseif n == 28
    # F_4, §7.6.2, p. 109. These are a simple system in the displayed line
    # system: the first two roots are long and the last two are short.
    return _lt_root_group(QQ, [
      QQ.([0, 1, -1, 0]), QQ.([0, 0, 1, -1]), QQ.([0, 0, 0, 1]),
      (1//2) .* QQ.([1, -1, -1, -1]),
    ])
  elseif n in (29, 31)
    # N_4 and O_4, §7.6.2, pp. 109--110: respectively D_4^(4) and B_4^(4),
    # with the orbit of the line (1,1,1,1) adjoined. Its scalar (1+i)/2 is
    # immaterial for the reflection. The additional coordinate reflection in
    # O_4 enlarges G(4,4,4) to G(4,2,4), not to G(4,1,4).
    R, x = polynomial_ring(QQ)
    K, i = number_field(x^2 + 1, "i")
    roots = _lt_adjacent_roots(K, 4)
    push!(roots, K.([1, -i, 0, 0]), K.([1, 1, 1, 1]))
    if n == 31
      push!(roots, K.([1, 0, 0, 0]))
    end
    return _lt_root_group(K, roots)
  elseif n in (33, 34)
    # K_5 and K_6, §7.6.1, pp. 107--108. Use the alternate K_6 realization
    # consisting of D_6^(2) and the orbit of (1,1,1,1,1,theta). For K_5 use
    # the book's orthonormal five-dimensional realization, not an arbitrary
    # basis of its six-dimensional hyperplane.
    R, x = polynomial_ring(QQ)
    K, omega = number_field(x^2 + x + 1, "omega")
    if n == 33
      R, x = polynomial_ring(K)
      K, sqrt2 = number_field(x^2 - 2, "sqrt2")
      omega = K(omega)
      theta = omega - omega^2
      roots = _lt_d_roots(K, 4, 5)
      push!(roots, K.([0, 0, 0, 0, 1]), K.([1, 1, 1, theta, sqrt2]))
    else
      theta = omega - omega^2
      roots = _lt_d_roots(K, 6)
      push!(roots, K.([1, 1, 1, 1, 1, theta]))
    end
    return _lt_root_group(K, roots)
  elseif n in (35, 36, 37)
    # E_6, E_7, E_8, §7.6.1, pp. 106--107. E_8 is D_8 plus the orbit of
    # (1,...,1). For E_6 take the subsystem orthogonal to e6-e7,e7-e8, and
    # use (e6+e7+e8)/sqrt(3) as the sixth orthonormal coordinate. For E_7
    # take the subsystem orthogonal to e7-e8 and use (e7+e8)/sqrt(2).
    # These are coordinate permutations of the E_6 subsystem in the book and
    # one of its equivalent E_7 subsystems; the omitted directions are fixed.
    if n == 35
      R, x = polynomial_ring(QQ)
      K, sqrt3 = number_field(x^2 - 3, "sqrt3")
      roots = _lt_d_roots(K, 5, 6)
      push!(roots, K.([1, 1, 1, 1, 1, sqrt3]))
    elseif n == 36
      R, x = polynomial_ring(QQ)
      K, sqrt2 = number_field(x^2 - 2, "sqrt2")
      roots = _lt_d_roots(K, 6, 7)
      push!(roots, K.([0, 0, 0, 0, 0, 0, 1]), K.([1, 1, 1, 1, 1, 1, sqrt2]))
    else
      K = QQ
      roots = _lt_d_roots(K, 8)
      push!(roots, fill(one(K), 8))
    end
    return _lt_root_group(K, roots)
  end
  error("Invalid exceptional Shephard--Todd number")
end

# Restrict the A_(n-1) permutation realization (§7.5, p. 104) to the sum-zero
# subspace. A Householder reflection H in u=(1,...,1,1-sqrt(n)) sends the row
# (1,...,1) to sqrt(n)*e_n. Its first n-1 rows are therefore an orthonormal
# basis of the sum-zero subspace. This exact basis is derived using the reflection
# formula in Definition 1.17, equation (1.18) (p. 11); it needs only Q(sqrt(n)).
function _lt_symmetric_group(n::Int)
  if n == 2
    return matrix_group([matrix(QQ, 1, 1, [-1])])
  end
  square, sqrt_n = is_square_with_sqrt(ZZ(n))
  if square
    K = QQ
    sqrt_n = QQ(sqrt_n)
  else
    R, x = polynomial_ring(QQ)
    K, sqrt_n = number_field(x^2 - n, "sqrt$n")
  end
  V = vector_space(K, n)
  u = fill(one(K), n)
  u[n] -= sqrt_n
  H = matrix(unitary_reflection(V(u)))
  roots = [[H[i, j] - H[i + 1, j] for j in 1:n-1] for i in 1:n-1]
  return _lt_root_group(K, roots)
end

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

    lines = [V([K(1), K(0), K(0)] * A) for A in cycshift]
    lines = vcat(lines, [1//2 * V([K(1), tau, tau^-1] * A) for A in cycshift])
    lines = vcat(lines, [1//2 * V([K(1), tau, -tau^-1] * A) for A in cycshift])
    lines = vcat(lines, [1//2 * V([K(1), -tau, tau^-1] * A) for A in cycshift])
    lines = vcat(lines, [1//2 * V([K(1), -tau, -tau^-1] * A) for A in cycshift])

    refls = Set([unitary_reflection(l) for l in lines])

    # The proof of Theorem 8.10 (page 144) gives these generating roots:
    # the line system H3 is the star-closure of the following lines.
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

  elseif 24 <= n <= 37

    W = _complex_reflection_group_LT_higher_rank(n)

  else

    error("LT model not yet implemented")

  end

  return W

end

function complex_reflection_group_LT(t::Tuple)
  (m,p,n) = t

  # G(m,p,n) is defined in Lehrer--Taylor, §2.2, pp. 24--25.

  # For the symmetric group case we do not want permutation matrices
  # since this is not an irreducible reflection representation
  if m == 1 && p == 1
    return _lt_symmetric_group(n)
  end

  # Reflection generators: Lehrer--Taylor, Lemma 2.8 and Proposition 2.10,
  # pp. 26--27. Their row-action matrices are decorated permutation matrices.
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
