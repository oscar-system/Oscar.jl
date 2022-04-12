@testset "Matroids" begin
    @testset "standard examples" begin
        for (M, values) in ((uniform_matroid(3,5), [3, 5, 10]),
			    (uniform_matroid(0,2), [0, 2, 1]),
			    (fano_matroid(), [3,7,28]),
			    (non_fano_matroid(), [3,7,29]),
			    (pappus_matroid(), [3,9,75]),
			    (non_pappus_matroid(), [3,9,76]),
			    (vamos_matroid(), [4,8,65]))
	    @test M isa Matroid
	    @test rank(M) == values[1]
        @test length(M) == values[2]
	    @test length(bases(M)) == values[3]
        end      
    end

    @testset "constructors (bases)" begin
	    mb1 = matroid_from_revlex_basis_encoding("*****0",2,4)
        @test mb1 isa Matroid
        @test rank(mb1) == 2
        @test length(mb1) == 4
        @test length(bases(mb1)) == 5
        @test rank(mb1,[3,4]) == 1

	    @test_throws ErrorException matroid_from_revlex_basis_encoding("*****",2,4)
	    @test_throws ErrorException matroid_from_revlex_basis_encoding("123",2,4)
	
	    #bases
	    mb2 = matroid_from_bases([[1,2],[1,3],[1,4],[2,3],[2,4]], 4)
        @test mb2 isa Matroid
        @test sort(bases(mb1)) == sort(bases(mb2))

        B = Set{Set{Int}}([Set([1,2]), Set([1,3]), Set([1,4]), Set([2,3]), Set([2,4])])
        mb3 = matroid_from_bases(B, 4)
        @test mb3 isa Matroid
        @test sort(bases(mb1)) == sort(bases(mb3))

        mb4 = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j'])
        @test mb4 isa Matroid
        @test rank(mb4) == 2
        @test length(mb4) == 4
        @test length(bases(mb4)) == 5
        @test rank(mb4,['i','j']) == 1

        mb5 = matroid_from_bases(Set{Set{Int}}([Set()]),1)
        @test mb5 isa Matroid
        @test length(bases(mb5)) == 1

        mb6 = matroid_from_bases(Set{Set{Int}}([Set()]),0)
        @test mb6 isa Matroid
        @test length(bases(mb6)) == 1

	    @test_throws ErrorException matroid_from_bases(Set{Set{Int}}([]),1)
	    @test_throws ErrorException matroid_from_bases([[1,2],[3]], 3)
	    @test_throws ErrorException matroid_from_bases([[1,2],[1,3],[1,4],[3,4]], 4)
	    @test_throws ErrorException matroid_from_bases([[1,2],[1,3],[2,4],[3,4]], 3)
	    @test_throws ErrorException matroid_from_bases([[1,2]], [1,2,1])

	    #TODO matroids from nonbases
    end

    #TODO circuits
    #TODO rest

    @testset "properties" begin
    #isomorphic_matroid
    @test rank(isomorphic_matroid(fano_matroid(), ["001","010","011","100","101","110","111"])) == 3
    @test_throws ErrorException isomorphic_matroid(fano_matroid(), [1,2,3])
    @test_throws ErrorException isomorphic_matroid(fano_matroid(), [1,2,3,1,2,3,4])

    @test isisomorphic(fano_matroid(), isomorphic_matroid(fano_matroid(), ["001","010","011","100","101","110","111"]))
    @test  !isisomorphic(fano_matroid(), non_fano_matroid())


    N = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j'])
	O = matroid_from_bases(Set{Set{Int}}([Set()]),0)

    R, q = PolynomialRing(ZZ, 'q')
    for (M, values) in ((uniform_matroid(3,5), [17, 2, 3, 4, inf, R(q^3-5q^2+10q-6)]),
                        (uniform_matroid(0,2), [1, 1, 0, 1, 1, R(0)]),
                        (vamos_matroid(), [79, 7, 3, 4, 3, R(q^4-8q^3+28q^2-51q+30)]),
                        (N, [5, 3, 2, 2, 2, R(q^2-3q+2)]),
                        (direct_sum(N,N), [25, 9, 2, 2 ,2, R(q^4-6q^3+13q^2-12q+4)]),
                        (O, [1, 1, 0, inf, inf, R(1)]))

        @test length(flats(M)) == values[1]
        @test length(cyclic_flats(M)) == values[2]
        r = rank(M)
        if(r!=0)
            @test Set{Vector{Any}}(flats(M,rank(M)-1)) == Set{Vector{Any}}(hyperplanes(M))
        end
        @test vertical_connectivity(M) == values[3]
        @test girth(M) == values[4]
        @test tutte_connectivity(M) == values[5]
        @test characteristic_polynomial(M) == values[6]
    end      

    @test_throws ErrorException flats(uniform_matroid(2,3),4)
    @test_throws ErrorException cyclic_flats(uniform_matroid(2,3),-1)

    @test fundamental_circuit(N,Set([1, 2]),'i') == [1, 2, 'i']
    @test fundamental_circuit(N,Set([1,'i',]),'j') == ['i','j']
    @test_throws ErrorException fundamental_circuit(N,Set(['i','j']),1)
    @test_throws ErrorException fundamental_circuit(N,Set([1, 2]), 3)
    @test_throws ErrorException fundamental_circuit(N,Set([1, 2]), 1)


    end


end
