kQ = monoid_algebra_from_lattice([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
a, b, c, d = gens(kQ)

#some longer examples
I = ideal(kQ, [a^2*b, c^2])
@time inj_res = injective_res(M, 3)

I = ideal(kQ, [a^2*b, c^2, d*a^4])
@time inj_res = injective_res(I, 3)
