@testset "StandardWalk" begin
    R1, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
            
    I1 = ideal([x^2 + y*z, x*y + z^2]) 
    G1 = groebner_basis(I1)

    I2 = ideal([z^3 - z^2 - 4])
    G2 = groebner_basis(I2)

    I3 = ideal([x^6 - 1])
    G3 = groebner_basis(I3)

    @test is_groebner_basis(groebner_walk(I1); ordering=lex(R1))
    #@test bounding_vectors(G1) == Vector{ZZRingElem}.([[1,1,-2],[2,-1,-1], [-1,2,-1]])
    #@test next_weight(G1, ZZ.([1,1,1]), ZZ.([1,0,0])) == ZZ.([1,1,1])

    @test is_groebner_basis(groebner_walk(I2); ordering=lex(R1))
    #@test bounding_vectors(G2) == ZZ.([[0,0,1],[0,0,3]])

    @test is_groebner_basis(groebner_walk(I3); ordering=lex(R1))
    #@test bounding_vectors(G3) == ZZ.([[6,0,0]])
    #@test next_weight(G3, ZZ.([1,1,1]), ZZ.([1,0,0])) == ZZ.([1,0,0])
    
    G = [x + 2*y + 2*z - 1, x^2 - x + 2*y^2 + 2*z^2, 2*x*y + 2*y*z - y];
    I = ideal(G);
    @test leading_ideal(groebner_walk(I, lex(R1))) == leading_ideal(groebner_basis(I, ordering = lex(R1)))
end 

