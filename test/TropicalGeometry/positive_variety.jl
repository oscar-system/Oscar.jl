@testset "src/TropicalGeometry/positive_variety.jl" begin

    C = matrix(QQ,[[-3,1,-1,-2,2],[-1,1,-1,-1,1]])
    R,x = polynomial_ring(QQ,ncols(C))
    nu = tropical_semiring_map(QQ)
    I = ideal(C*gens(R))
    TropPlusI = positive_tropical_variety(I,nu)
    @test n_maximal_polyhedra(TropPlusI) == 5

    I = ideal([x[1]^2-x[2]^2,x[3]^3-x[4]^3])
    TropPlusI = positive_tropical_variety(I,nu)
    @test n_maximal_polyhedra(TropPlusI) == 1

end

@testset "non-positive binomial ideal" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  nu = tropical_semiring_map(QQ)
  TropPlusI = positive_tropical_variety(ideal(R, [x+y]), nu)
  @test ambient_dim(TropPlusI) == 2
  @test n_maximal_polyhedra(TropPlusI) == 0
end
