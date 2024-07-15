@testset "src/TropicalGeometry/matrix.jl" begin

    for minOrMax in (min,max)
        T = tropical_semiring(minOrMax)
        @testset "det(A::Generic.MatSpaceElem{<:TropicalSemiringElem})" begin
            A = matrix(T,[1 2; 3 5])
            @test det(A) == (minOrMax==min ? T(5) : T(6))
        end
    end

end

@testset "tropical general position" begin
	A = matrix(tropical_semiring(),[1 0;0 1])
	@test tropically_generic(A,min) == true
	A = matrix(tropical_semiring(max),[1 0;0 1])
	@test_throws ArgumentError tropically_generic(A,min)
end
