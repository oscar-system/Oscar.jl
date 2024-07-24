@testset "src/TropicalGeometry/matrix.jl" begin

    for minOrMax in (min,max)
        T = tropical_semiring(minOrMax)
        @testset "det(A::Generic.MatSpaceElem{<:TropicalSemiringElem})" begin
            A = Matrix(T[1 2; 3 5]) # (julia) matrix over tropical numbers
            @test det(A) == (minOrMax==min ? T(5) : T(6))
        end
    end

end

@testset "tropical general position" begin
  A = matrix(tropical_semiring(),[1 0;0 1])
  @test is_tropically_generic(A) == true
  A = matrix(tropical_semiring(max),[1 5;0 0;0 0])
  @test is_tropically_generic(A) == false
  @test is_tropically_generic(transpose(A)) == false
end
