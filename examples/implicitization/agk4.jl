#implicitization of Bezier surface


using Oscar
using GroebnerWalk
using BenchmarkTools

R, (x,y,z,u,v) = polynomial_ring(QQ, ["x","y","z","u","v"])

o1 = matrix_ordering(R, [1 1 1 0 0; 0 0 0 1 1; 0 0 0 1 0; 1 1 0 0 0; 1 0 0 0 0])
o2 = matrix_ordering(R, [0 0 0 1 1; 1 1 1 0 0; 1 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0])

I = ideal([
    u + u^2 - 2*v - 2*u^2*v + 2*u*v^2 - x,
    -6*u + 2*v + v^2 - 5*v^3 + 2*u*v^2 - 4*u^2*v^2 - y,
    -2 + 2*u^2 + 6*v - 3*u^2*v^2 - z
])
set_verbosity_level(:groebner_walk, 1)

G = groebner_basis(I, ordering = o2 )
t_standard = @elapsed G = groebner_walk(I, o1, o2; algorithm=:standard)
# t_generic = @elapsed G = groebner_walk(I, o2, o1; algorithm=:generic)
# t_perturbed = @elapsed G = groebner_walk(I, o2, o1; algorithm=:perturbed)