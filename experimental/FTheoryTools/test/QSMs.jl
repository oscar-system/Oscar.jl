qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 283))
as = ambient_space(qsm_model)
x1, x2, x3, x4, x5, x6, x7, v, e3, e2, u, e4, e1, w = gens(cohomology_ring(as))
X1 = cohomology_class(as, x1)
X2 = cohomology_class(as, x2)
X3 = cohomology_class(as, x3)
X4 = cohomology_class(as, x4)
X5 = cohomology_class(as, x5)
X6 = cohomology_class(as, x6)
X7 = cohomology_class(as, x7)
E3 = cohomology_class(as, e3)
E4 = cohomology_class(as, e4)
E1 = cohomology_class(as, e1) + E4
E2 = cohomology_class(as, e2) + E3
U = cohomology_class(as, u)
V = cohomology_class(as, v)
W = cohomology_class(as, w)
H = W+E1
cy = cohomology_class(toric_divisor_class(as, degree(hypersurface_equation(qsm_model))))
c2_B = cohomology_class(as, 4*x4*x6 + 8*x5*x6 + 14*x6*x7 + 5*x7^2)
Kbar = cohomology_class(ambient_space(qsm_model), x1+x2+x3+x4+x5+x6+x7)
Kbar3 = kbar3(qsm_model)
fg_not_breaking = special_flux_family(qsm_model, not_breaking = true, check = false)
g4_exp = flux_instance(fg_not_breaking, [3], [])

@testset "Test properties of the QSM" begin
  @test H == U + E1 + E2 + E4 == V + E2 + E3
  @test chern_classes(qsm_model)[3] == c2_B - 7 * E4^2 - E1 * Kbar - E2 * Kbar - E3 * Kbar - E4 * Kbar + 3 * H * Kbar
  @test euler_characteristic(qsm_model) == integrate((3 * Kbar * (4 * c2_B + 5 * Kbar^2)) * V * cy)
  @test d3_tadpole_constraint(g4_exp) == 12 + 5//8 * Kbar3 - 45//(2 * Kbar3)
end
