# Mimic the Laurent polynomials in one variable ℤ[t, t⁻¹] as a ℤ-module:
F = InfiniteDirectSum(QQ, ZZ)

v = 2*F[1] + 7*F[3] # 2 ⋅ t + 7 ⋅ t³
w = 9*F[1] + 122*F[5] # 9 ⋅ t + 122 ⋅ t⁵

# implement the general differentiation rule: tᵏ ↦ k ⋅ tᵏ⁻¹
diff_map(i::elem_type(ZZ)) = QQ(i)*F[i-1]

# construct the endomorphism of ℤ-modules
D = InfiniteDirectSumHom(F, F, diff_map)

# apply it to the above elements:
@show v
@show w
@show D(v+w)
