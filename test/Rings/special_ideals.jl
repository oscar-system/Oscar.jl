@testset "Special Ideals" begin
    @testset "Katsura Ideals" begin
        I = katsura(ZZ, 5)
        R = parent(gens(I)[1])
        (x1, x2, x3, x4, x5, x6) = gens(R)
        @test gens(I) == [
            x1 + 2 * x2 + 2 * x3 + 2 * x4 + 2 * x5 + 2 * x6 - 1,
            x1^2 - x1 + 2 * x2^2 + 2 * x3^2 + 2 * x4^2 + 2 * x5^2 + 2 * x6^2,
            2 * x1 * x2 + 2 * x2 * x3 - x2 + 2 * x3 * x4 + 2 * x4 * x5 + 2 * x5 * x6,
            2 * x1 * x3 + x2^2 + 2 * x2 * x4 + 2 * x3 * x5 - x3 + 2 * x4 * x6,
            2 * x1 * x4 + 2 * x2 * x3 + 2 * x2 * x5 + 2 * x3 * x6 - x4,
            2 * x1 * x5 + 2 * x2 * x4 + 2 * x2 * x6 + x3^2 - x5,
        ]
    end
end
