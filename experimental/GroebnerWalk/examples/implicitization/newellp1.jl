#=
    Newell's teapot, patch 1
    A 2-dimensional ideal from Tran. "Efficient Groebner walk conversion for implicitization of geometric objects" (2004)
=#
using Oscar

R, (x,y,z,u,v) = polynomial_ring(QQ, ["x","y","z","u","v"])

I = ideal([
    -x + 7//5 - 231//125 * v^2 + 39//80 * u^2 - 1//5 * u^3 + 99//400 * u * v^2 - 1287//2000 * u^2 * v^2 + 33//125 * u^3 * v^2 - 3//16 * u + 56//125 * v^3 - 3//50 * u * v^3 + 39//250 * u^2 * v^3 - 8//125 * u^3 * v^3,
    -y + 63//125 * v^2 - 294//125 * v + 56//125 * v^3 - 819//1000 * u^2 * v + 42//125 * u^3 * v - 3//50 * u * v^3 + 351//2000 * u^2 * v^2 + 39//250 * u^2 * v^3 - 9//125 * u^3 * v^2 - 8//125 * u^3 * v^3,
    -z + 12//5 - 63//160 * u^2 + 63//160 * u
])

Ginit = groebner_basis(I)
# "Optimal" choice of orderings, as stated in Tran 2004
o1 = matrix_ordering(R, [1 1 1 0 0; 0 0 0 1 1; 0 0 0 1 0; 1 1 0 0 0; 1 0 0 0 0])
o2 = matrix_ordering(R, [0 0 0 1 1; 1 1 1 0 0; 1 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0])

set_verbosity_level(:groebner_walk, 1)
G = groebner_walk(I, o2, o1; algorithm=:standard)

