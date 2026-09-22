# Tests for /src/rings/sagbi.jl
using Test
using Oscar
include("../../src/Rings/sagbi.jl")

@testset "Subduction modulo generators" begin
    
    R, (x,y) = polynomial_ring(QQ, ['x', 'y'])

    B = [x^2 - x, y+1]
    @test subduct(x^2*y + x*y - 1, B) == 2*x*y - 1 
    @test subduct(x^2*y-x*y-y, B) == 0

    B = [x+y^2, x*y+y^3]
    @test subduct(x*y, B) == x*y
    @test subduct(x*y+y^3, B) == 0
    @test subduct(x*y, B; ordering=lex(R)) == -y^3
    @test subduct(x*y+y^3, B; ordering=lex(R)) == 0
end

@testset "Checking SAGBI bases" begin
    Qx, x = QQ["x"];
    K, a = number_field(x^2-2, "a")
    R, (x,y,z) = polynomial_ring(K, ['x','y','z'])

    @test is_sagbi([x,y]) === true
    @test is_sagbi([x^2-x, y+1]) === true
    @test is_sagbi([x+y, x^2 + y^2, a*z ]) === false

    # Plucker coordinates for Gr(4,2)
    # Plucker coordinates w/ increasing index form a SAGBI basis
    R, vars = polynomial_ring(QQ, [
        :x11, :x12, :x13, :x14,
        :x21, :x22, :x23, :x24
    ])

    x11, x12, x13, x14, x21, x22, x23, x24 = vars
    p12 = x11*x22 - x12*x21
    p13 = x11*x23 - x13*x21
    p14 = x11*x24 - x14*x21
    p23 = x12*x23 - x13*x22
    p24 = x12*x24 - x14*x22
    p34 = x13*x24 - x14*x23

    B = [p12, p13, p14, p23, p24, p34]
    @test is_sagbi([b for b in B]) === true

    # An example which depends on the ordering:
    R, (x,y) = polynomial_ring(QQ, ['x','y'])
    B = [x+y^2, x*y+y^3]
    @test is_sagbi(B, ordering=degrevlex(R)) === false
    @test is_sagbi(B, ordering=lex(R)) === true
end

@testset "Computing SAGBI bases" begin
    R, (x,y) = polynomial_ring(QQ, ['x','y'])
    B = sagbi([x+y^2, x*y+y^3]; degree_bound=10)
    @test B.elements == [x+y^2, x*y+y^3, x^3 + 2*x^2*y^2 + x*y^4]
    @test B.sagbi_degree === -1
    @test is_sagbi(B) === true

    # QQ[x, xy-y^2, xy^2] has no finite SAGBI basis:
    # see section 4.11 from Robbiano and Sweedler "Subalgebra bases" 
    A = [x, x*y-y^2, x*y^2]
    B = sagbi(A; degree_bound=5)
    B_6 = [
        x, x*y - y^2, x*y^2, x*y^3 - 1//2*y^4,
        x*y^5 - 1//3*y^6,  x*y^4
    ]
    @test B.elements == B_6
    @test B.sagbi_degree === 6

    B = _initialize_sagbi_candidate(B_6)
    compute_sagbi_degree!(B)
    @test B.sagbi_degree === 6
end