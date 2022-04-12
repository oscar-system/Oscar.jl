@testset "Chow rings for matroids" begin

    #test ring by extracting the characteristic_polynomial
    for M in (uniform_matroid(4,4), uniform_matroid(3,5), fano_matroid(), vamos_matroid())
        A = chow_ring(M)
        @test A isa MPolyQuo
        proper_flats = flats(M)[2:size(flats(M))[1]-1]
        e = matroid_groundset(M)[1]
        a = sum([A[i] for i in _select([e],[],proper_flats)])
        b = sum([A[i] for i in _select([],[e],proper_flats)])
        list = [a^i*b^(rank(M)-i-1) for i in (0:rank(M)-1)]
        L1 = [abs(coeff(M.f,1)) for M in list]
        f = reduced_characteristic_polynomial(M)
        d = degree(f)
        L2 = [(-1)^(i+d)*coeff(f,i) for i in 0:d]
        @test L1==L2
    end

    @test_throws ErrorException chow_ring(uniform_matroid(0,2))
    @test_throws ErrorException chow_ring(uniform_matroid(1,1))

end
