
### POLYNOMIALRING
# EXAMPLE hasse_deriv (IZ == 0)
R, x = polynomial_ring(ZZ, 4, "x")
IZ = ideal(R, [zero(R)])
IX = ideal(R, [3*x[2], 4*x[4]^3])
hasse_deriv(IZ, IX)

R, x = polynomial_ring(ZZ, 8, "x")
IZ = ideal(R, [zero(R)])
IX = ideal(R, [5*x[2]*x[4]^2*x[5]])
hasse_deriv(IZ, IX)

# EXAMPLE hasse_deriv (IZ != 0)
R, x = polynomial_ring(ZZ, 4, "x")
IX = ideal(R, [3*x[2],4*x[4]^3])
IZ = ideal(R, [zero(0)])
# y = x
systemOfParameters = [1]
M = matrix(R, 1, 1, [R(0)])
hasse_deriv(IZ, IX, y, M) # working, since i call hasse_deriv(IZ, IX) for IZ = 0

R, x = polynomial_ring(ZZ, 4, "x")
IZ = ideal(R, [x[1]])
IX = IZ + ideal(R, [3*x[2], 4*x[4]^3])
# y = x[2:4]
systemOfParameters = [2, 3, 4]
M = matrix(R, 1, 1, [R(1)])
hasse_deriv(IZ, IX, y, M) # working 

R, x = polynomial_ring(ZZ, 8, "x")
IZ = ideal(R, [x[1]])
# y = x[2:8]
systemOfParameters = [2, 3, 4, 5, 6, 7, 8]
M = matrix(R, 1, 1, [R(1)])
IX = IZ + ideal(R, [x[2]*x[3]]) # working
IX = IZ + ideal(R, [x[2]^2]) # working
IX = IZ + ideal(R, [x[2]^2 + 3]) # working
IX = IZ + ideal(R, [x[2]^2*x[3]]) # working
IX = IZ + ideal(R, [x[4]^5, x[2]*x[3]*x[7]]) # working
hasse_deriv(IZ, IX, y, M)

R, x = polynomial_ring(ZZ, 8, "x")
IZ = ideal(R, [x[1]*x[2] - 1])
# y = x[2:8]
systemOfParameters = [2, 3, 4, 5, 6, 7, 8]
M = matrix(R, 1, 1, [x[2]])
IX = IZ # working # error (IZ and IX cannot be equal.)
IX = IZ + ideal(R, [x[1]]) # working
IX = IZ + ideal(R, [x[3]]) # working
IX = IZ + ideal(R, [x[2]^2]) # working
IX = IZ + ideal(R, [x[3]^2]) # working
IX = IZ + ideal(R, [x[2]*x[3]]) # working
IX = IZ + ideal(R, [x[4]^2*x[3]]) # working
IX = IZ + ideal(R, [x[5]^3, x[6]^4]) # working
IX = IZ + ideal(R, [x[6]^4, 7]) # working
IX = IZ + ideal(R, [x[8]^6]) # working
IX = IZ + ideal(R, [x[3]^3 + x[2]^2]) # not working # expected result might be wrong
IX = IZ + ideal(R, [x[4]^3 + x[6]^4]) # working
hasse_deriv(IZ, IX, y, M)

R, x = polynomial_ring(ZZ, 8, "x")
IZ = ideal(R, [x[1]*x[2] - 1])
# y = vcat(x[1], x[3:8])
systemOfParameters = [1, 3, 4, 5, 6, 7, 8]
M = matrix(R, 1, 1, [x[1]])
IX = IZ + ideal(R, [x[2]]) # working
IX = IZ + ideal(R, [x[2], x[3]]) # working
IX = IZ + ideal(R, [x[3]]) # working
IX = IZ + ideal(R, [x[3]^2]) # working
IX = IZ + ideal(R, [x[3]^2 - 5]) # working
IX = IZ + ideal(R, [x[4]^2*x[3]]) # working
IX = IZ + ideal(R, [x[5]^3, x[6]^4]) # working
IX = IZ + ideal(R, [x[2]*x[3]]) # working
IX = IZ + ideal(R, [x[2]*x[3]^2]) # working
IX = IZ + ideal(R, [x[6]^4, 7]) # working
IX = IZ + ideal(R, [x[6]^4]) # working
IX = IZ + ideal(R, [x[3]^3 + x[4]^2]) # working
IX = IZ + ideal(R, [x[3]^3 + x[2]^2]) # working
IX = IZ + ideal(R, [x[4]^3 + x[6]^4]) # working
hasse_deriv(IZ, IX, y, M)

### FAKTORRING / QUOTIENTENRING
R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
I = ideal(R, [x - 1])
RQ, _ = quo(R, I)
IZ = ideal(RQ, [y - 1])
IX = IZ + ideal(RQ, [-y^3 + x^2])
systemOfParameters = [3]
M = matrix(RQ, 1, 1, [z])

hasse_deriv(IZ, IX, systemOfParameters, M) # debuggen (endlosschleife), ? x-1 und y-1 sorgen für IX = 0 in RQ ?
# sollte ich auf IX = 0 prüfen? wo prüfe ich bisher auf IX = 0 ? kommt hasse_deriv mit IX = 0 nicht klar?

### LOKALISIERTER RING
R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
P = ideal(R, [x])
U = complement_of_prime_ideal(P)
Rloc, iota = localization(R, U)
IZ = ideal(Rloc, [y - 1])
IX = IZ + ideal(Rloc, [-y^3 + x^2])
systemOfParameters = [1, 3]
M = matrix(Rloc, 1, 1, [z])

hasse_deriv(IZ, IX, systemOfParameters, M)

# Oscar.MPolyLocalizedIdeal
R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
P = ideal(R, [x])
U = complement_of_prime_ideal(P)
Rloc, iota = localization(R, U)
IZ = ideal(Rloc, [zero(Rloc)])
IX = ideal(Rloc, [-y^3 + x^2])


# MPolyQuoIdeal
# T, t = polynomial_ring(QQ, "t")
# K, a =  number_field(2*t^2-1, "a")
R, (x, y) = polynomial_ring(QQ, ["x", "y"])
I = ideal(R, [x - 1])  # I = ideal(R, [x-1, x-a])
RQ, _ = quo(R, I)
IZ = ideal(RQ, [zero(RQ)])
IX = ideal(RQ, [-y^3 + x^2])

R, (x, y) = polynomial_ring(QQ, ["x", "y"])
I = ideal(R, [x^3 - 1])
RQ, phi = quo(R, I) 
P = ideal(R, [y])
U = complement_of_prime_ideal(P)
RQL, iota = localization(RQ, U) 
IZ = ideal(RQL, [zero(RQL)])
IX = ideal(RQL, [-y^3 + x^2])  # RESOVLED: IZ and IX are equal (need to understand why ^^ ) # Need to find an example with IZ != IX 