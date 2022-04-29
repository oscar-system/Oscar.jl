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

        mb4 = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],Set([1,2,'i','j']))
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

        mb7 = matroid_from_nonbases([[3,4]],4)
        @test mb7 isa Matroid
        @test sort(bases(mb2)) == sort(bases(mb7))

        mb8 = matroid_from_nonbases(Set([Set([3,4])]),Set([3,4,5,6]))
        @test is_isomorphic(mb7,mb8) == true

        @test_throws ErrorException matroid_from_nonbases(Set{Set{Int}}([]),1)
        @test_throws ErrorException matroid_from_nonbases([[1,2],[3]], 3)
        @test_throws ErrorException matroid_from_nonbases([[2,3],[2,4]], 4)
        @test_throws ErrorException matroid_from_nonbases([[1,2],[1,3],[2,4],[3,4]], 3)
        @test_throws ErrorException matroid_from_nonbases([[1,2]], [1,2,1])
    end

    @testset "constructors (circuts, hyperplanes)" begin
        mc1 = matroid_from_circuits([[3,4],[1,2,3],[1,2,4]],4)
        @test length(bases(mc1)) == 5

        mc2 = matroid_from_circuits([[3,4],[1,2,3],[1,2,4]],Set([2,3,1,4]))
        @test is_isomorphic(mc1,mc2)

        @test_throws ErrorException matroid_from_circuits([[1,2],[1,3],[2,4],[3,4]], 3)
        @test_throws ErrorException matroid_from_circuits([[1,2]], [1,2,1])

        mc3 = matroid_from_hyperplanes([[1],[2],[3,4]],4)
        @test length(bases(mc3)) == 5

        mc4 = matroid_from_hyperplanes(Set([Set([1]),Set([2]),Set(['a','b'])]),Set(['b',2,1,'a']))
        @test length(cyclic_flats(mc3)) == 3
        @test is_isomorphic(mc3,mc4)

        @test_throws ErrorException matroid_from_hyperplanes(Set{Set{Int}}([]),1)
        @test_throws ErrorException matroid_from_hyperplanes([[1,2],[1,3],[2,4],[3,4]], 3)
        @test_throws ErrorException matroid_from_hyperplanes([[1,2]], [1,2,1])
    end

    @testset "constructors (matroid from matrices and graphs)" begin
        mm1 = matroid_from_matrix_columns(matrix(QQ,[1 0 0 1 1; 0 1 0 0 1]))
        @test length(bases(mm1)) == 5
        @test loops(mm1) == [3]

        mm2 = matroid_from_matrix_columns(matrix(GF(5),[1 2]'))
        @test bases(mm2) == [[1]]
        @test circuits(mm2) == []

        mm3 = matroid_from_matrix_columns(matrix(GF(3),[1 0 0 1 1]))
        @test bases(mm3) == [[1],[4],[5]]

        mm4 = matroid_from_matrix_rows(matrix(GF(5),[1 2]))
        @test is_isomorphic(mm2,mm4)

        mm5 = matroid_from_matrix_columns(matrix(QQ,[[]]))
        @test mm5 isa Matroid
        @test rank(mm5) == 0


        mm6 = matroid_from_matrix_columns(matrix(GF(2),[0 0; 0 0]))
        @test mm6 isa Matroid
        @test rank(mm6) == 0

        g1 = Oscar.Graphs.complete_graph(5);
        mg1 = cycle_matroid(g1)
        @test mg1 isa Matroid
        @test length(bases(mg1)) == 125

        g2 = Oscar.Graphs.Graph{Oscar.Graphs.Undirected}(1)
        Oscar.Graphs.add_edge!(g2,1,1)
        mg2 = cycle_matroid(g2)
        @test length(mg2) == 1
        @test rank(mg2) == 0

        mg3 = bond_matroid(g1)
        @test mg3 isa Matroid
        @test  length(circuits(mg3)) == 15

        mg4 = cocycle_matroid(g1)
        @test mg4 isa Matroid
        @test sort(bases(mg3)) == sort(bases(mg4))

        mg5 = dual_matroid(mg1)
        @test mg5 isa Matroid
        # the next two lines delete the parent object of the dual matroid
        # by reassiging mg1 and calling the garbage collector
        mg1 = nothing
        GC.gc()
        @test bases(mg5) isa Vector
        @test length(hyperplanes(mg5)) == 37
        @test matroid_groundset(mg5) == 1:10
    end

    @testset "matroid modifications" begin
        g3 = Oscar.Graphs.complete_graph(4);
        m1 = cycle_matroid(g3)
        m2 = direct_sum(m1,m1)
        @test m2 isa Matroid
        @test length(hyperplanes(m2)) == 14

        m3 = direct_sum([m1,m1])
        @test m3 isa Matroid
        @test sort(bases(m2)) == sort(bases(m3))

        m4 = deletion(m2, Set([2,"1'",5]))
        @test length(bases(m4)) == 32

        m5 = deletion(m2, 5)
        @test length(m5) == 11

        m6 = deletion(m2,Set([]))
        @test sort(bases(m6)) == sort(bases(m2))
        @test_throws ErrorException deletion(m2,Set([1,9]))
        @test_throws ErrorException deletion(m2,42)

        m7 = restriction(m2,Set([1,2]))
        @test m7 isa Matroid
        @test rank(m7) == 2

        m8 = contraction(m2,Set([1,2]))
        @test m8 isa Matroid
        @test rank(m8) == 4

        m9 = contraction(m2,[1,1])
        m10 = contraction(m2,[1])
        @test sort(bases(m9)) == sort(bases(m10))

        m11 = minor(m2,Set([1,2]),[3,4,3])
        @test m11 isa Matroid
        @test length(bases(m11)) == 32

        @test_throws ErrorException minor(m2,Set([1,2]),[1,4,3])

        m12 = principal_extension(fano_matroid(),Set([1,2,3]),8)
        @test m12 isa Matroid
        @test length(nonbases(principal_extension(fano_matroid(),Set([1,2,3]),8))) == 10

        @test_throws ErrorException principal_extension(fano_matroid(),Set([1,2]),5)

        @test length(circuits(free_extension(fano_matroid(),8))) == 42
        @test length(circuits(series_extension(fano_matroid(), 1,8))) == 14
        @test length(circuits(parallel_extension(fano_matroid(), 1,8))) == 22
        @test_throws ErrorException parallel_extension(fano_matroid(), 8,9)

        @test length(circuits(all_subsets_matroid(3))) == 17

        @test length(bases(projective_plane(3))) == 234
        @test_throws ErrorException projective_plane(6)

        @test length(bases(projective_geometry(2,3))) == 234
        @test length(bases(projective_geometry(3,3))) == 63180
        @test_throws ErrorException projective_geometry(1,3)
        @test_throws ErrorException projective_geometry(3,6)

        @test length(bases(affine_geometry(2,3))) == 72
        @test length(bases(affine_geometry(3,3))) == 40743
        @test_throws ErrorException affine_geometry(1,3)
        @test_throws ErrorException affine_geometry(3,6)
    end

    @testset "properties" begin
    #isomorphic_matroid
        @test rank(isomorphic_matroid(fano_matroid(), ["001","010","011","100","101","110","111"])) == 3
        @test_throws ErrorException isomorphic_matroid(fano_matroid(), [1,2,3])
        @test_throws ErrorException isomorphic_matroid(fano_matroid(), [1,2,3,1,2,3,4])

        @test is_isomorphic(fano_matroid(), isomorphic_matroid(fano_matroid(), ["001","010","011","100","101","110","111"]))
        @test  !is_isomorphic(fano_matroid(), non_fano_matroid())


        N = matroid_from_bases([[1,2],[1,'i'],[1,'j'],[2,'i'],[2,'j']],[1,2,'i','j'])
        NN = direct_sum(N,N)
        O = matroid_from_bases(Set{Set{Int}}([Set()]),0)

        R, q = PolynomialRing(ZZ, 'q')
        for (M, values) in ((uniform_matroid(3,5), [2, 17, 2, 3, 4, inf, R(q^3-5q^2+10q-6), false, false, false, true]),
                            (uniform_matroid(0,2), [2, 1, 1, 0, 1, 1, R(0), true, true, true, false]),
                            (vamos_matroid(), [4, 79, 7, 3, 4, 3, R(q^4-8q^3+28q^2-51q+30), false, false, false, true]),
                            (N, [2, 5, 3, 2, 2, 2, R(q^2-3q+2), true, true, true, false]),
                            (NN, [4, 25, 9, 2, 2 ,2, R(q^4-6q^3+13q^2-12q+4), true, true, true, false]),
                            (uniform_matroid(5,5), [0, 32, 1, 1, inf, 1, R(q^5-5q^4+10q^3-10q^2+5q-1), true, true, true, true]),
                            (O, [0, 1, 1, 0, inf, inf, R(1), true, true, true, true]))

            @test nullity(M,matroid_groundset(M)) == values[1]
            @test length(flats(M)) == values[2]
            @test length(cyclic_flats(M)) == values[3]
            r = rank(M)
            if(r!=0)
                @test Set{Vector{Any}}(flats(M,rank(M)-1)) == Set{Vector{Any}}(hyperplanes(M))
                @test corank(M,[1,2]) == rank(dual_matroid(M),[1,2])
            end
            @test corank(M,matroid_groundset(M)) == rank(dual_matroid(M),matroid_groundset(M))
            @test vertical_connectivity(M) == values[4]
            @test girth(M) == values[5]
            @test tutte_connectivity(M) == values[6]
            @test characteristic_polynomial(M) == values[7]
            @test length(cobases(M)) == length(bases(M))
            @test cohyperplanes(M) == [setdiff(matroid_groundset(M),set) for set in circuits(M)]
            @test is_regular(M) == values[8]
            @test is_binary(M) == values[9]
            @test is_ternary(M) == values[10]
            @test n_connected_components(M) == length(connected_components(M))
            @test is_connected(M) == (n_connected_components(M)<2)
            @test is_simple(M) == values[11]
        end      
    
        @test_throws ErrorException flats(uniform_matroid(2,3),4)
        @test_throws ErrorException cyclic_flats(uniform_matroid(2,3),-1)
        @test nullity(N,['i','j']) == 1
        @test nullity(N,Set(['i','j'])) == 1

        @test fundamental_circuit(N,Set([1, 2]),'i') == [1, 2, 'i']
        @test fundamental_circuit(N,Set([1,'i',]),'j') == ['i','j']
        @test_throws ErrorException fundamental_circuit(N,Set(['i','j']),1)
        @test_throws ErrorException fundamental_circuit(N,Set([1, 2]), 3)
        @test_throws ErrorException fundamental_circuit(N,Set([1, 2]), 1)

        @test fundamental_cocircuit(N,Set([1, 2]),1) == [1,'i','j']
        @test fundamental_cocircuit(N,Set([1,'i',]),'i') == [2,'i','j']
        @test_throws ErrorException fundamental_cocircuit(N,Set(['i','j']),'i')
        @test_throws ErrorException fundamental_cocircuit(N,Set([1, 2]), 'i')
        @test_throws ErrorException fundamental_cocircuit(N,Set([1, 2]), 3)

        @test independent_sets(N) == [[],[1],[2],['i'],['j'],[1,'j'], [2,'j'], [1,'i'], [2,'i'], [1,2]]
        @test independent_sets(O) == [[]]
        @test independent_sets(uniform_matroid(0,2)) == [[]]

        @test spanning_sets(N) == [[1,'j'], [2,'j'], [1,2], [1,'i'], [2,'i'], [1,2,'i'], [1,2,'j'], [1,'i','j'], [2,'i','j'], [1,2,'i','j']]
        @test spanning_sets(O) == [[]]
        @test spanning_sets(uniform_matroid(0,2)) == [[], [1], [2], [1, 2]]

        @test cobases(N) == [['i','j'], [2,'j'], [2,'i'], [1,'j'], [1,'i']]

        @test charpoly(N) == R(q^2-3q+2)

        @test is_clutter(bases(N)) == true
        @test is_clutter(circuits(N)) == true
        @test is_clutter(hyperplanes(N)) == true
        @test is_clutter([[1],[1,2],[3,4]]) == false

        @test is_regular(uniform_matroid(2,4)) == false
        @test is_binary(uniform_matroid(2,4)) == false
        @test is_ternary(uniform_matroid(2,4)) == true

        @test [bases(M) for M in direct_sum_components(NN)] == [bases(N),bases(isomorphic_matroid(N,["1'","2'","i'","j'"]))]

        @test loops(direct_sum(uniform_matroid(0,2),uniform_matroid(2,2))) == [1,2]
        @test coloops(direct_sum(uniform_matroid(0,2),uniform_matroid(2,2))) == ["1'","2'"]
        @test is_loopless(NN)
        @test is_coloopless(NN)
        @test !is_loopless(direct_sum(uniform_matroid(0,2),N))
        @test !is_coloopless(direct_sum(uniform_matroid(2,2),N))

        @test connectivity_function(N,[1,'i']) == 2
        @test connectivity_function(N,Set([1,'i'])) == 2
        @test connectivity_function(NN,[1,"1'"]) == 2
        @test connectivity_function(N,matroid_groundset(N)) == 0
        @test connectivity_function(N,[]) == 0
        @test connectivity_function(O,[]) == 0
        @test_throws KeyError connectivity_function(N,[1,3])

        K23 = cycle_matroid(Graphs.complete_bipartite_graph(2,3));
        @test !is_vertical_k_separation(K23,2,[1,3,6])
        @test is_vertical_k_separation(K23,3,[1,3,6])
        @test !is_vertical_k_separation(K23,4,[1,3,6])
        @test is_vertical_k_separation(K23,2,[1,2])
        @test !is_vertical_k_separation(K23,3,[1,2])

        @test !is_k_separation(K23,2,[1,3,6])
        @test is_k_separation(K23,3,[1,3,6])
        @test !is_k_separation(K23,4,[1,3,6])
        @test is_k_separation(K23,2,[1,2])
        @test !is_k_separation(K23,3,[1,2])

        @test is_isomorphic(matroid_from_revlex_basis_encoding("**0***",2,4),N)
        @test !is_isomorphic(matroid_from_nonbases([[1,2,3],[3,4,5]],6),matroid_from_nonbases([[1,2,3],[4,5,6]],6))
        @test !is_isomorphic(N,principal_extension(N,[],3))

        @test is_minor(N,N)
        @test is_minor(uniform_matroid(1,3),N)
        @test is_minor(direct_sum(uniform_matroid(1,2),uniform_matroid(0,1)),N)
        @test is_minor(N,NN)
        @test !is_minor(NN,N)
        @test !is_minor(fano_matroid(), uniform_matroid(2,4))
   end
end
