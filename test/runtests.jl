using FTheoryTools
using Test
using Oscar

P3 = projective_space(NormalToricVariety,3)
Kbar = anticanonical_bundle(P3)
f = sum([rand(Int)*b for b in basis_of_global_sections(Kbar^4)]);
g = sum([rand(Int)*b for b in basis_of_global_sections(Kbar^6)]);
w = GlobalWeierstrassModel(f,g)

@testset "Dummy test" begin
    @test poly_f(w) == f
    @test poly_g(w) == g
end
