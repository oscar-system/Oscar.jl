@testset "src/TropicalGeometry/matrix.jl" begin

    for minOrMax in (min,max)
        T = tropical_semiring(minOrMax)
        @testset "det(A::Generic.MatSpaceElem{<:TropicalSemiringElem})" begin
            A = matrix(T,[1 2; 3 5]) # matrix over tropical numbers
            @test det(A) == (minOrMax==min ? T(5) : T(6))
            R, (x,y) = T[:x,:y]
            A = identity_matrix(R,2) # matrix over tropical polynomials
            @test det(A) == one(R)
        end
    end

end
