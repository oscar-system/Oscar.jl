@testset "src/TropicalGeometry/matrix.jl" begin

    for minOrMax in (min,max)
        T = tropical_semiring(minOrMax)
        @testset "det(A::Generic.MatSpaceElem{<:TropicalSemiringElem})" begin
            A = matrix(T,[1 2; 3 5])
            @test det(A) == (minOrMax==min ? T(5) : T(6))
        end
    end

end
