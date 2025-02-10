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

@testset "testing for polytropes" begin
  do_rows_form_polytrope = Oscar.do_rows_form_polytrope
  A = tropical_semiring()[0 0 1; 0 1 0; 0 3 3]
  @test do_rows_form_polytrope(A) == true
  B = tropical_semiring()[0 1 0; 0 0 1; 1 0 0]
  @test do_rows_form_polytrope(B) == false
  C = tropical_semiring()[0 0 1; 0 1 0; 0 3 3; 0 2 3; 0 3 2; 0 0 0]
  @test do_rows_form_polytrope(C) == true
  D = tropical_semiring()[0 1 0; 0 0 1; 2 1 0]
  @test do_rows_form_polytrope(D) == false
  E = tropical_semiring()[0 1 0; 0 0 0; 0 -1 0]
  @test do_rows_form_polytrope(E) == true
  F = tropical_semiring()[0 3 1 4 ;0 -1 -1 3; 0 1 -5 1 ; 0 2 -2 0]
  @test do_rows_form_polytrope(F) == true
  G = tropical_semiring()[0 3 1 4; 0 -1 -1 3; 0 1 -5 1; 0 2 -2 0; 0 0 -3 2]
  @test do_rows_form_polytrope(G) == true
  H = tropical_semiring()[0 0 0 0 0; 0 1 2 3 4; 0 2 4 6 8; 0 4 8 12 16; 0 5 10 15 20]
  @test do_rows_form_polytrope(H) == false
  end
