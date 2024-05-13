

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
y = x
M = matrix(R, 1, 1, [R(0)])
hasse_deriv(IZ, IX, y, M) # working, since i call hasse_deriv(IZ, IX) for IZ = 0

R, x = polynomial_ring(ZZ, 4, "x")
IZ = ideal(R, [x[1]])
IX = IZ + ideal(R, [3*x[2], 4*x[4]^3])
y = x[2:4]
M = matrix(R, 1, 1, [R(1)])
hasse_deriv(IZ, IX, y, M) # working 

R, x = polynomial_ring(ZZ, 8, "x")
IZ = ideal(R, [x[1]])
y = x[2:8]
M = matrix(R, 1, 1, [R(1)])
IX = IZ + ideal(R, [x[2]*x[3]]) # working
IX = IZ + ideal(R, [x[2]^2]) # working
IX = IZ + ideal(R, [x[2]^2 + 3]) # working
IX = IZ + ideal(R, [x[2]^2*x[3]]) # working
IX = IZ + ideal(R, [x[4]^5, x[2]*x[3]*x[7]]) # working
hasse_deriv(IZ, IX, y, M)

R, x = polynomial_ring(ZZ, 8, "x")
IZ = ideal(R, [x[1]*x[2] - 1])
y = x[2:8]
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
y = vcat(x[1], x[3:8])
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