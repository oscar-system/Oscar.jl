@testset "Matroid Strata in the Grassmannian" begin

    mb1 = uniform_matroid(2,5)
    mb2 = fano_matroid()
    mb3 = matroid_from_nonbases([[1,2,6],[1,4,5],[1,7,8],[2,3,5],[2,4,8],[3,4,7],[3,6,8],[5,6,7]], 8)
    mb4 = matroid_from_nonbases([[1,2,3,4,5]], 7)

    bmc1 = bases_matrix_coordinates(bases(mb1), [1,2])
    bmc2 = bases_matrix_coordinates(bases(mb2), [1,2,4])
    bmc3a = bases_matrix_coordinates(bases(mb3), [1,2,3])
    bmc3b = bases_matrix_coordinates(bases(mb3), [1,2,4])
    bmc4a = bases_matrix_coordinates(bases(mb4), [1,2,3,4,6])
    bmc4b = bases_matrix_coordinates(bases(mb4), [1,2,3,6,7]) 
    
    @testset "bases_matrix_coordinates" begin
	
        @test bmc1 == [[1,1],[2,1],[1,2],[2,2],[1,3],[2,3]]
        @test bmc2 == [[1,1],[2,1],[1,2],[3,2],[2,3],[3,3],[1,4],[2,4],[3,4]]
        @test bmc3a == [[1,1],[2,1],[3,1],[2,2],[3,2],[1,3],[2,3],[1,4],[2,4],[3,4],[1,5],[2,5],[3,5]]
        @test bmc3b == [[1,1],[2,1],[3,1],[1,2],[3,2],[1,3],[2,3],[1,4],[2,4],[3,4],[2,5],[3,5]]
        @test bmc4a == [[1,1],[2,1],[3,1],[4,1],[1,2],[2,2],[3,2],[4,2],[5,2]]
        @test bmc4b == [[1,1],[2,1],[3,1],[4,1],[5,1],[1,2],[2,2],[3,2],[4,2],[5,2]]
            
    end

    T1, t1 = PolynomialRing(QQ, "t")
    T2, (s2, t2) = PolynomialRing(GF(3), ["s", "t"]);
    QT1 = FractionField(T1)
    QT2 = FractionField(T2)
    K, a = NumberField(t1^2 + 1, "a")
    
    R1, x1, xd1 = make_polynomial_ring(bases(mb1), [1,2], ZZ)
    R2, x2, xd2 = make_polynomial_ring(bases(mb2), [1,2,4], GF(2))
    R3a, x3a, xd3a = make_polynomial_ring(bases(mb3), [1,2,3], QQ)
    R3b, x3b, xd3b = make_polynomial_ring(bases(mb3), [1,2,4], QT1)
    R4a, x4a, xd4a = make_polynomial_ring(bases(mb4), [1,2,3,4,6], QT2)
    R4b, x4b, xd4b = make_polynomial_ring(bases(mb4), [1,2,3,6,7], K)

    @testset "make_polynomial_ring" begin

        @test all(xd1[bmc1[a]] == x1[a]  for a in 1:6)
        @test all(xd2[bmc2[a]] == x2[a]  for a in 1:9)
        @test all(xd3a[bmc3a[a]] == x3a[a]  for a in 1:13)
        @test all(xd3b[bmc3b[a]] == x3b[a]  for a in 1:12)
        @test all(xd4a[bmc4a[a]] == x4a[a]  for a in 1:9)
        @test all(xd4b[bmc4b[a]] == x4b[a]  for a in 1:10)
        
    end

    S1 = MatrixSpace(R1, 2, 3); M1 = make_coordinate_matrix_no_identity(2, 5, bmc1, R1, x1, xd1)
    S2 = MatrixSpace(R2, 3, 4); M2 = make_coordinate_matrix_no_identity(3, 7, bmc2, R2, x2, xd2)
    S3a = MatrixSpace(R3a, 3, 5); M3a = make_coordinate_matrix_no_identity(3, 8, bmc3a, R3a, x3a, xd3a)
    S3b = MatrixSpace(R3b, 3, 5); M3b = make_coordinate_matrix_no_identity(3, 8, bmc3b, R3b, x3b, xd3b)
    S4a = MatrixSpace(R4a, 5, 2); M4a = make_coordinate_matrix_no_identity(5, 7, bmc4a, R4a, x4a, xd4a)
    S4b = MatrixSpace(R4b, 5, 2); M4b = make_coordinate_matrix_no_identity(5, 7, bmc4b, R4b, x4b, xd4b)    
    
    @testset "make_coordinate_matrix_no_identity" begin
        @test M1 == S1([x1[1] x1[3] x1[5];
                        x1[2] x1[4] x1[6]])
        @test M2 == S2([x2[1] x2[3] R2(0) x2[7];
                        x2[2] R2(0) x2[5] x2[8];
                        R2(0) x2[4] x2[6] x2[9]])
        @test M3a == S3a([x3a[1] R3a(0) x3a[6] x3a[8] x3a[11];
                          x3a[2] x3a[4] x3a[7] x3a[9] x3a[12];
                          x3a[3] x3a[5] R3a(0) x3a[10] x3a[13] ])
        @test M3b == S3b([x3b[1] x3b[4] x3b[6] x3b[8] R3b(0);
                          x3b[2] R3b(0) x3b[7] x3b[9] x3b[11];
                          x3b[3] x3b[5] R3b(0) x3b[10] x3b[12] ])
        @test M4a == S4a([x4a[1] x4a[5];
                          x4a[2] x4a[6];
                          x4a[3] x4a[7];
                          x4a[4] x4a[8];
                          R4a(0) x4a[9]])
        @test M4b == S4b([x4b[1] x4b[6];
                          x4b[2] x4b[7];
                          x4b[3] x4b[8];
                          x4b[4] x4b[9];
                          x4b[5] x4b[10]])
    end

    

#    @testset "interlace_columns" begin

#        @test
        
#    end

    S1 = MatrixSpace(R1, 2, 5); X1 = make_coordinate_matrix(2, 5, bmc1, [1,2], R1, x1, xd1)
    S2 = MatrixSpace(R2, 3, 7); X2 = make_coordinate_matrix(3, 7, bmc2, [1,2,4], R2, x2, xd2)
    S3a = MatrixSpace(R3a, 3, 8); X3a = make_coordinate_matrix(3, 8, bmc3a, [1,2,3], R3a, x3a, xd3a)
    S3b = MatrixSpace(R3b, 3, 8); X3b = make_coordinate_matrix(3, 8, bmc3b, [1,2,4], R3b, x3b, xd3b)
    S4a = MatrixSpace(R4a, 5, 7); X4a = make_coordinate_matrix(5, 7, bmc4a, [1,2,3,4,6], R4a, x4a, xd4a)
    S4b = MatrixSpace(R4b, 5, 7); X4b = make_coordinate_matrix(5, 7, bmc4b, [1,2,3,6,7], R4b, x4b, xd4b)
    
    
    @testset "make_coordinate_matrix" begin

        @test X1 == S1([R1(1) R1(0) x1[1] x1[3] x1[5];
                  R1(0) R1(1) x1[2] x1[4] x1[6]])

        @test X2 == S2([R2(1) R2(0) x2[1] R2(0) x2[3] R2(0) x2[7];
                  R2(0) R2(1) x2[2] R2(0) R2(0) x2[5] x2[8];
                  R2(0) R2(0) R2(0) R2(1) x2[4] x2[6] x2[9]])

        @test X3a == S3a([R3a(1) R3a(0) R3a(0) x3a[1] R3a(0) x3a[6] x3a[8] x3a[11];
                          R3a(0) R3a(1) R3a(0) x3a[2] x3a[4] x3a[7] x3a[9] x3a[12];
                          R3a(0) R3a(0) R3a(1) x3a[3] x3a[5] R3a(0) x3a[10] x3a[13]])

        @test X3b == S3b([R3b(1) R3b(0) x3b[1] R3b(0) x3b[4] x3b[6] x3b[8] R3b(0);
                          R3b(0) R3b(1) x3b[2] R3b(0) R3b(0) x3b[7] x3b[9] x3b[11];
                          R3b(0) R3b(0) x3b[3] R3b(1) x3b[5] R3b(0) x3b[10] x3b[12] ])

        @test X4a == S4a([R4a(1) R4a(0) R4a(0) R4a(0) x4a[1] R4a(0) x4a[5];
                          R4a(0) R4a(1) R4a(0) R4a(0) x4a[2] R4a(0) x4a[6];
                          R4a(0) R4a(0) R4a(1) R4a(0) x4a[3] R4a(0) x4a[7];
                          R4a(0) R4a(0) R4a(0) R4a(1) x4a[4] R4a(0) x4a[8];
                          R4a(0) R4a(0) R4a(0) R4a(0) R4a(0) R4a(1) x4a[9]])

        @test X4b == S4b([R4b(1) R4b(0) R4b(0) x4b[1] x4b[6] R4b(0) R4b(0);
                          R4b(0) R4b(1) R4b(0) x4b[2] x4b[7] R4b(0) R4b(0);
                          R4b(0) R4b(0) R4b(1) x4b[3] x4b[8] R4b(0) R4b(0);
                          R4b(0) R4b(0) R4b(0) x4b[4] x4b[9] R4b(1) R4b(0);
                          R4b(0) R4b(0) R4b(0) x4b[5] x4b[10] R4b(0) R4b(1)])        
    end

    D1 = bases_determinants(2, 5, bases(mb1), bmc1, [1,2], R1, x1, xd1)
    D2 = bases_determinants(3, 7, bases(mb2), bmc2, [1,2,4], R2, x2, xd2)
    D3a = bases_determinants(3, 8, bases(mb3), bmc3a, [1,2,3], R3a, x3a, xd3a)
    D3b = bases_determinants(3, 8, bases(mb3), bmc3b, [1,2,4], R3b, x3b, xd3b)
    D4a = bases_determinants(5, 7, bases(mb4), bmc4a, [1,2,3,4,6], R4a, x4a, xd4a)
    D4b = bases_determinants(5, 7, bases(mb4), bmc4b, [1,2,3,6,7], R4b, x4b, xd4b)


    @testset "bases_determinants" begin

        @test D1 == [R1(1), xd1[[2,1]], xd1[[2,2]], xd1[[2,3]], -xd1[[1,1]], -xd1[[1,2]], -xd1[[1,3]], xd1[[1,1]] * xd1[[2,2]] - xd1[[2,1]] * xd1[[1,2]], xd1[[1,1]] * xd1[[2,3]] - xd1[[2,1]] * xd1[[1,3]],  xd1[[1,2]] * xd1[[2,3]] - xd1[[2,2]] * xd1[[1,3]]]

        @test D2 == [xd2[[3, 3]]*xd2[[1, 4]], xd2[[1, 1]]*xd2[[2, 3]]*xd2[[3, 4]] + xd2[[1, 1]]*xd2[[3, 3]]*xd2[[2, 4]] + xd2[[2, 1]]*xd2[[3, 3]]*xd2[[1, 4]], xd2[[2, 3]]*xd2[[1, 4]], xd2[[1, 2]]*xd2[[2, 3]]*xd2[[3, 4]] + xd2[[1, 2]]*xd2[[3, 3]]*xd2[[2, 4]] + xd2[[3, 2]]*xd2[[2, 3]]*xd2[[1, 4]], xd2[[3, 2]]*xd2[[2, 4]], xd2[[1, 1]]*xd2[[3, 2]]*xd2[[2, 4]] + xd2[[2, 1]]*xd2[[1, 2]]*xd2[[3, 4]] + xd2[[2, 1]]*xd2[[3, 2]]*xd2[[1, 4]], xd2[[1, 2]]*xd2[[2, 4]], xd2[[2, 4]], xd2[[1, 4]], xd2[[2, 1]]*xd2[[3, 4]], xd2[[1, 1]]*xd2[[3, 4]], xd2[[3, 4]], xd2[[3, 2]]*xd2[[2, 3]], xd2[[1, 2]]*xd2[[3, 3]], xd2[[1, 2]]*xd2[[2, 3]], xd2[[2, 3]], xd2[[1, 1]]*xd2[[2, 3]], xd2[[2, 1]]*xd2[[3, 3]], xd2[[1, 1]]*xd2[[3, 3]], xd2[[3, 3]], xd2[[1, 2]], xd2[[2, 1]]*xd2[[1, 2]], xd2[[2, 1]]*xd2[[3, 2]], xd2[[1, 1]]*xd2[[3, 2]], xd2[[3, 2]], xd2[[2, 1]], xd2[[1, 1]], R2(1)]

        @test D3a ==  [R3a(1), xd3a[[3, 1]],  xd3a[[3, 2]],  xd3a[[3, 4]],  xd3a[[3, 5]],  -xd3a[[2, 1]],  -xd3a[[2, 2]],  -xd3a[[2, 3]],  -xd3a[[2, 4]],  -xd3a[[2, 5]],  -xd3a[[3, 1]]*xd3a[[2, 3]],  xd3a[[2, 1]]*xd3a[[3, 4]] - xd3a[[3, 1]]*xd3a[[2, 4]],  xd3a[[2, 1]]*xd3a[[3, 5]] - xd3a[[3, 1]]*xd3a[[2, 5]],  -xd3a[[3, 2]]*xd3a[[2, 3]],  xd3a[[2, 2]]*xd3a[[3, 4]] - xd3a[[3, 2]]*xd3a[[2, 4]],  xd3a[[2, 2]]*xd3a[[3, 5]] - xd3a[[3, 2]]*xd3a[[2, 5]],  xd3a[[2, 3]]*xd3a[[3, 4]],  xd3a[[2, 3]]*xd3a[[3, 5]],  xd3a[[1, 1]],  xd3a[[1, 3]],  xd3a[[1, 4]],  xd3a[[1, 5]],  -xd3a[[1, 1]]*xd3a[[3, 2]],  xd3a[[3, 1]]*xd3a[[1, 3]],  -xd3a[[1, 1]]*xd3a[[3, 4]] + xd3a[[3, 1]]*xd3a[[1, 4]],  xd3a[[3, 2]]*xd3a[[1, 3]],  xd3a[[3, 2]]*xd3a[[1, 4]],  xd3a[[3, 2]]*xd3a[[1, 5]],  -xd3a[[1, 3]]*xd3a[[3, 4]],  -xd3a[[1, 3]]*xd3a[[3, 5]],  -xd3a[[1, 4]]*xd3a[[3, 5]] + xd3a[[3, 4]]*xd3a[[1, 5]],  xd3a[[1, 1]]*xd3a[[2, 2]],  xd3a[[1, 1]]*xd3a[[2, 3]] - xd3a[[2, 1]]*xd3a[[1, 3]],  xd3a[[1, 1]]*xd3a[[2, 5]] - xd3a[[2, 1]]*xd3a[[1, 5]],  -xd3a[[2, 2]]*xd3a[[1, 3]],  -xd3a[[2, 2]]*xd3a[[1, 4]],  -xd3a[[2, 2]]*xd3a[[1, 5]],  xd3a[[1, 3]]*xd3a[[2, 4]] - xd3a[[2, 3]]*xd3a[[1, 4]],  xd3a[[1, 4]]*xd3a[[2, 5]] - xd3a[[2, 4]]*xd3a[[1, 5]],  -xd3a[[1, 1]]*xd3a[[3, 2]]*xd3a[[2, 3]] + xd3a[[2, 1]]*xd3a[[3, 2]]*xd3a[[1, 3]] - xd3a[[3, 1]]*xd3a[[2, 2]]*xd3a[[1, 3]],  xd3a[[1, 1]]*xd3a[[2, 2]]*xd3a[[3, 4]] - xd3a[[1, 1]]*xd3a[[3, 2]]*xd3a[[2, 4]] + xd3a[[2, 1]]*xd3a[[3, 2]]*xd3a[[1, 4]] - xd3a[[3, 1]]*xd3a[[2, 2]]*xd3a[[1, 4]],  xd3a[[1, 1]]*xd3a[[2, 2]]*xd3a[[3, 5]] - xd3a[[1, 1]]*xd3a[[3, 2]]*xd3a[[2, 5]] + xd3a[[2, 1]]*xd3a[[3, 2]]*xd3a[[1, 5]] - xd3a[[3, 1]]*xd3a[[2, 2]]*xd3a[[1, 5]],  xd3a[[1, 1]]*xd3a[[2, 3]]*xd3a[[3, 4]] - xd3a[[2, 1]]*xd3a[[1, 3]]*xd3a[[3, 4]] + xd3a[[3, 1]]*xd3a[[1, 3]]*xd3a[[2, 4]] - xd3a[[3, 1]]*xd3a[[2, 3]]*xd3a[[1, 4]],  xd3a[[1, 1]]*xd3a[[2, 3]]*xd3a[[3, 5]] - xd3a[[2, 1]]*xd3a[[1, 3]]*xd3a[[3, 5]] + xd3a[[3, 1]]*xd3a[[1, 3]]*xd3a[[2, 5]] - xd3a[[3, 1]]*xd3a[[2, 3]]*xd3a[[1, 5]],  xd3a[[1, 1]]*xd3a[[2, 4]]*xd3a[[3, 5]] - xd3a[[1, 1]]*xd3a[[3, 4]]*xd3a[[2, 5]] - xd3a[[2, 1]]*xd3a[[1, 4]]*xd3a[[3, 5]] + xd3a[[2, 1]]*xd3a[[3, 4]]*xd3a[[1, 5]] + xd3a[[3, 1]]*xd3a[[1, 4]]*xd3a[[2, 5]] - xd3a[[3, 1]]*xd3a[[2, 4]]*xd3a[[1, 5]],  -xd3a[[2, 2]]*xd3a[[1, 3]]*xd3a[[3, 5]] + xd3a[[3, 2]]*xd3a[[1, 3]]*xd3a[[2, 5]] - xd3a[[3, 2]]*xd3a[[2, 3]]*xd3a[[1, 5]],  -xd3a[[2, 2]]*xd3a[[1, 4]]*xd3a[[3, 5]] + xd3a[[2, 2]]*xd3a[[3, 4]]*xd3a[[1, 5]] + xd3a[[3, 2]]*xd3a[[1, 4]]*xd3a[[2, 5]] - xd3a[[3, 2]]*xd3a[[2, 4]]*xd3a[[1, 5]],  xd3a[[1, 3]]*xd3a[[2, 4]]*xd3a[[3, 5]] - xd3a[[1, 3]]*xd3a[[3, 4]]*xd3a[[2, 5]] - xd3a[[2, 3]]*xd3a[[1, 4]]*xd3a[[3, 5]] + xd3a[[2, 3]]*xd3a[[3, 4]]*xd3a[[1, 5]]]

        @test D3b == [xd3b[[3, 1]], R3b(1), xd3b[[3, 2]], xd3b[[3, 4]], xd3b[[3, 5]], xd3b[[2, 1]], xd3b[[2, 1]]*xd3b[[3, 2]], -xd3b[[3, 1]]*xd3b[[2, 3]], xd3b[[2, 1]]*xd3b[[3, 4]] - xd3b[[3, 1]]*xd3b[[2, 4]], xd3b[[2, 1]]*xd3b[[3, 5]] - xd3b[[3, 1]]*xd3b[[2, 5]], -xd3b[[2, 3]], -xd3b[[2, 4]], -xd3b[[2, 5]], -xd3b[[3, 2]]*xd3b[[2, 3]], -xd3b[[3, 2]]*xd3b[[2, 4]], -xd3b[[3, 2]]*xd3b[[2, 5]], xd3b[[2, 3]]*xd3b[[3, 4]], xd3b[[2, 3]]*xd3b[[3, 5]], -xd3b[[1, 1]], xd3b[[3, 1]]*xd3b[[1, 3]], -xd3b[[1, 1]]*xd3b[[3, 4]] + xd3b[[3, 1]]*xd3b[[1, 4]], -xd3b[[1, 1]]*xd3b[[3, 5]], xd3b[[1, 2]], xd3b[[1, 3]], xd3b[[1, 4]], xd3b[[3, 2]]*xd3b[[1, 3]], -xd3b[[1, 2]]*xd3b[[3, 4]] + xd3b[[3, 2]]*xd3b[[1, 4]], -xd3b[[1, 2]]*xd3b[[3, 5]], -xd3b[[1, 3]]*xd3b[[3, 4]], -xd3b[[1, 3]]*xd3b[[3, 5]], -xd3b[[1, 4]]*xd3b[[3, 5]], xd3b[[2, 1]]*xd3b[[1, 2]], -xd3b[[1, 1]]*xd3b[[2, 3]] + xd3b[[2, 1]]*xd3b[[1, 3]], -xd3b[[1, 1]]*xd3b[[2, 5]], -xd3b[[1, 1]]*xd3b[[3, 2]]*xd3b[[2, 3]] + xd3b[[2, 1]]*xd3b[[3, 2]]*xd3b[[1, 3]] + xd3b[[3, 1]]*xd3b[[1, 2]]*xd3b[[2, 3]], -xd3b[[1, 1]]*xd3b[[3, 2]]*xd3b[[2, 4]] - xd3b[[2, 1]]*xd3b[[1, 2]]*xd3b[[3, 4]] + xd3b[[2, 1]]*xd3b[[3, 2]]*xd3b[[1, 4]] + xd3b[[3, 1]]*xd3b[[1, 2]]*xd3b[[2, 4]], -xd3b[[1, 1]]*xd3b[[3, 2]]*xd3b[[2, 5]] - xd3b[[2, 1]]*xd3b[[1, 2]]*xd3b[[3, 5]] + xd3b[[3, 1]]*xd3b[[1, 2]]*xd3b[[2, 5]], xd3b[[1, 1]]*xd3b[[2, 3]]*xd3b[[3, 4]] - xd3b[[2, 1]]*xd3b[[1, 3]]*xd3b[[3, 4]] + xd3b[[3, 1]]*xd3b[[1, 3]]*xd3b[[2, 4]] - xd3b[[3, 1]]*xd3b[[2, 3]]*xd3b[[1, 4]], xd3b[[1, 1]]*xd3b[[2, 4]]*xd3b[[3, 5]] - xd3b[[1, 1]]*xd3b[[3, 4]]*xd3b[[2, 5]] - xd3b[[2, 1]]*xd3b[[1, 4]]*xd3b[[3, 5]] + xd3b[[3, 1]]*xd3b[[1, 4]]*xd3b[[2, 5]], xd3b[[1, 2]]*xd3b[[2, 3]], xd3b[[1, 2]]*xd3b[[2, 4]], xd3b[[1, 2]]*xd3b[[2, 5]], xd3b[[1, 3]]*xd3b[[2, 4]] - xd3b[[2, 3]]*xd3b[[1, 4]], xd3b[[1, 3]]*xd3b[[2, 5]], xd3b[[1, 4]]*xd3b[[2, 5]], xd3b[[1, 2]]*xd3b[[2, 3]]*xd3b[[3, 5]] + xd3b[[3, 2]]*xd3b[[1, 3]]*xd3b[[2, 5]], xd3b[[1, 2]]*xd3b[[2, 4]]*xd3b[[3, 5]] - xd3b[[1, 2]]*xd3b[[3, 4]]*xd3b[[2, 5]] + xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[2, 5]], xd3b[[1, 3]]*xd3b[[2, 4]]*xd3b[[3, 5]] - xd3b[[1, 3]]*xd3b[[3, 4]]*xd3b[[2, 5]] - xd3b[[2, 3]]*xd3b[[1, 4]]*xd3b[[3, 5]]]

        @test D4a == [R4a(1), xd4a[[5, 2]], xd4a[[4, 1]], xd4a[[4, 1]]*xd4a[[5, 2]], 2*xd4a[[4, 2]], 2*xd4a[[3, 1]], 2*xd4a[[3, 1]]*xd4a[[5, 2]], xd4a[[3, 2]], 2*xd4a[[3, 1]]*xd4a[[4, 2]] + xd4a[[4, 1]]*xd4a[[3, 2]], xd4a[[2, 1]], xd4a[[2, 1]]*xd4a[[5, 2]], 2*xd4a[[2, 2]], xd4a[[2, 1]]*xd4a[[4, 2]] + 2*xd4a[[4, 1]]*xd4a[[2, 2]], 2*xd4a[[2, 1]]*xd4a[[3, 2]] + xd4a[[3, 1]]*xd4a[[2, 2]], 2*xd4a[[1, 1]], 2*xd4a[[1, 1]]*xd4a[[5, 2]], xd4a[[1, 2]], 2*xd4a[[1, 1]]*xd4a[[4, 2]] + xd4a[[4, 1]]*xd4a[[1, 2]], xd4a[[1, 1]]*xd4a[[3, 2]] + 2*xd4a[[3, 1]]*xd4a[[1, 2]], 2*xd4a[[1, 1]]*xd4a[[2, 2]] + xd4a[[2, 1]]*xd4a[[1, 2]] ]

        @test D4b == [-xd4b[[5, 1]], xd4b[[4, 1]], -xd4b[[5, 2]], xd4b[[4, 2]], R4b(1), -xd4b[[3, 1]]*xd4b[[5, 2]] + xd4b[[5, 1]]*xd4b[[3, 2]], xd4b[[3, 1]]*xd4b[[4, 2]] - xd4b[[4, 1]]*xd4b[[3, 2]], xd4b[[3, 1]], xd4b[[3, 2]], xd4b[[2, 1]]*xd4b[[5, 2]] - xd4b[[5, 1]]*xd4b[[2, 2]], -xd4b[[2, 1]]*xd4b[[4, 2]] + xd4b[[4, 1]]*xd4b[[2, 2]], -xd4b[[2, 1]], -xd4b[[2, 2]], xd4b[[2, 1]]*xd4b[[3, 2]] - xd4b[[3, 1]]*xd4b[[2, 2]], -xd4b[[1, 1]]*xd4b[[5, 2]] + xd4b[[5, 1]]*xd4b[[1, 2]], xd4b[[1, 1]]*xd4b[[4, 2]] - xd4b[[4, 1]]*xd4b[[1, 2]], xd4b[[1, 1]], xd4b[[1, 2]], -xd4b[[1, 1]]*xd4b[[3, 2]] + xd4b[[3, 1]]*xd4b[[1, 2]], xd4b[[1, 1]]*xd4b[[2, 2]] - xd4b[[2, 1]]*xd4b[[1, 2]]]
 
    end


    SG1 = localizing_semigroup(2, 5, bases(mb1), bmc1, [1,2], R1, x1, xd1)
    SG2 = localizing_semigroup(3, 7, bases(mb2), bmc2, [1,2,4], R2, x2, xd2)
    SG3a = localizing_semigroup(3, 8, bases(mb3), bmc3a, [1,2,3], R3a, x3a, xd3a)
    SG3b = localizing_semigroup(3, 8, bases(mb3), bmc3b, [1,2,4], R3b, x3b, xd3b)
    SG4a = localizing_semigroup(5, 7, bases(mb4), bmc4a, [1,2,3,4,6], R4a, x4a, xd4a)
    SG4b = localizing_semigroup(5, 7, bases(mb4), bmc4b, [1,2,3,6,7], R4b, x4b, xd4b)
    
    @testset "localizing_semigroup" begin
        
        @test all(a in SG1 for a in D1)
        @test all(a in SG2 for a in D2)
        @test all(a in SG3a for a in D3a)
        @test all(a in SG3b for a in D3b)
        @test all(a in SG4a for a in D4a)
        @test all(a in SG4b for a in D4b)

        
    end

    MS1 = matroid_stratum_matrix_coordinates_given_ring(2,5,mb1,ZZ,[1,2],R1,x1,xd1)
    MS2 = matroid_stratum_matrix_coordinates_given_ring(3,7,mb2,GF(2),[1,2,4],R2,x2,xd2)
    MS3a = matroid_stratum_matrix_coordinates_given_ring(3,8,mb3,QQ,[1,2,3],R3a,x3a,xd3a)
    MS3b = matroid_stratum_matrix_coordinates_given_ring(3,8,mb3,QT1,[1,2,4],R3b,x3b,xd3b)
    MS4a = matroid_stratum_matrix_coordinates_given_ring(5,7,mb4,QT2,[1,2,3,4,6],R4a,x4a,xd4a)
    MS4b = matroid_stratum_matrix_coordinates_given_ring(5,7,mb4,K,[1,2,3,6,7],R4b,x4b,xd4b)

    
    @testset "matroid_stratum_matrix_coordinates_given_ring" begin
        @test MS1[2] == ideal(R1, [])

        @test MS2[2] == ideal(R2, [xd2[[2, 3]]*xd2[[3, 4]] + xd2[[3, 3]]*xd2[[2, 4]], xd2[[1, 2]]*xd2[[3, 4]] + xd2[[3, 2]]*xd2[[1, 4]], xd2[[1, 1]]*xd2[[2, 4]] + xd2[[2, 1]]*xd2[[1, 4]], xd2[[1, 1]]*xd2[[3, 2]]*xd2[[2, 4]] + xd2[[2, 1]]*xd2[[1, 2]]*xd2[[3, 4]], xd2[[1, 1]]*xd2[[3, 2]]*xd2[[2, 3]] + xd2[[2, 1]]*xd2[[1, 2]]*xd2[[3, 3]]])

        @test MS3a[2] == ideal(R3a, [-xd3a[[2, 4]]*xd3a[[3, 5]] + xd3a[[3, 4]]*xd3a[[2, 5]], -xd3a[[1, 3]]*xd3a[[2, 5]] + xd3a[[2, 3]]*xd3a[[1, 5]], -xd3a[[1, 1]]*xd3a[[3, 5]] + xd3a[[3, 1]]*xd3a[[1, 5]], -xd3a[[1, 1]]*xd3a[[2, 4]] + xd3a[[2, 1]]*xd3a[[1, 4]], -xd3a[[2, 1]]*xd3a[[3, 2]] + xd3a[[3, 1]]*xd3a[[2, 2]], -xd3a[[1, 1]]*xd3a[[2, 3]]*xd3a[[3, 5]] + xd3a[[3, 1]]*xd3a[[1, 3]]*xd3a[[2, 5]], -xd3a[[2, 2]]*xd3a[[1, 4]]*xd3a[[3, 5]] + xd3a[[2, 2]]*xd3a[[3, 4]]*xd3a[[1, 5]] + xd3a[[3, 2]]*xd3a[[1, 4]]*xd3a[[2, 5]], -xd3a[[1, 1]]*xd3a[[2, 4]]*xd3a[[3, 5]] + xd3a[[2, 1]]*xd3a[[3, 4]]*xd3a[[1, 5]] + xd3a[[3, 1]]*xd3a[[1, 4]]*xd3a[[2, 5]], -xd3a[[2, 2]]*xd3a[[1, 4]]*xd3a[[3, 5]] + xd3a[[3, 2]]*xd3a[[2, 4]]*xd3a[[1, 5]], -xd3a[[1, 1]]*xd3a[[2, 2]]*xd3a[[3, 5]] + xd3a[[2, 1]]*xd3a[[3, 2]]*xd3a[[1, 5]], -xd3a[[1, 1]]*xd3a[[2, 3]]*xd3a[[3, 4]] + xd3a[[3, 1]]*xd3a[[1, 3]]*xd3a[[2, 4]], xd3a[[2, 2]]*xd3a[[1, 3]]*xd3a[[3, 4]] - xd3a[[3, 2]]*xd3a[[1, 3]]*xd3a[[2, 4]] + xd3a[[3, 2]]*xd3a[[2, 3]]*xd3a[[1, 4]], -xd3a[[1, 1]]*xd3a[[2, 3]]*xd3a[[3, 4]] + xd3a[[2, 1]]*xd3a[[1, 3]]*xd3a[[3, 4]] + xd3a[[3, 1]]*xd3a[[2, 3]]*xd3a[[1, 4]], -xd3a[[2, 2]]*xd3a[[1, 4]]*xd3a[[2, 5]]*xd3a[[3, 5]] + xd3a[[2, 2]]*xd3a[[2, 4]]*xd3a[[1, 5]]*xd3a[[3, 5]] + xd3a[[3, 2]]*xd3a[[1, 4]]*xd3a[[2, 5]]^2, -xd3a[[1, 1]]*xd3a[[2, 4]]*xd3a[[2, 5]]*xd3a[[3, 5]] + xd3a[[2, 1]]*xd3a[[2, 4]]*xd3a[[1, 5]]*xd3a[[3, 5]] + xd3a[[3, 1]]*xd3a[[1, 4]]*xd3a[[2, 5]]^2, xd3a[[2, 2]]^2*xd3a[[3, 5]]^2 - xd3a[[2, 2]]*xd3a[[3, 2]]*xd3a[[2, 5]]*xd3a[[3, 5]] + xd3a[[3, 2]]^2*xd3a[[2, 5]]^2, xd3a[[2, 1]]*xd3a[[2, 2]]*xd3a[[3, 5]]^2 - xd3a[[2, 1]]*xd3a[[3, 2]]*xd3a[[2, 5]]*xd3a[[3, 5]] + xd3a[[3, 1]]*xd3a[[3, 2]]*xd3a[[2, 5]]^2, -xd3a[[1, 1]]*xd3a[[2, 2]]*xd3a[[2, 5]]*xd3a[[3, 5]] + xd3a[[1, 1]]*xd3a[[3, 2]]*xd3a[[2, 5]]^2 + xd3a[[2, 1]]*xd3a[[2, 2]]*xd3a[[1, 5]]*xd3a[[3, 5]], xd3a[[2, 1]]^2*xd3a[[3, 5]]^2 - xd3a[[2, 1]]*xd3a[[3, 1]]*xd3a[[2, 5]]*xd3a[[3, 5]] + xd3a[[3, 1]]^2*xd3a[[2, 5]]^2, -xd3a[[1, 1]]*xd3a[[2, 1]]*xd3a[[2, 5]]*xd3a[[3, 5]] + xd3a[[1, 1]]*xd3a[[3, 1]]*xd3a[[2, 5]]^2 + xd3a[[2, 1]]^2*xd3a[[1, 5]]*xd3a[[3, 5]], -xd3a[[2, 2]]*xd3a[[2, 3]]*xd3a[[1, 4]]*xd3a[[3, 5]] + xd3a[[3, 2]]*xd3a[[1, 3]]*xd3a[[2, 4]]*xd3a[[2, 5]], xd3a[[2, 2]]^2*xd3a[[3, 4]]*xd3a[[3, 5]] - xd3a[[2, 2]]*xd3a[[3, 2]]*xd3a[[2, 4]]*xd3a[[3, 5]] + xd3a[[3, 2]]^2*xd3a[[2, 4]]*xd3a[[2, 5]], xd3a[[2, 1]]*xd3a[[2, 2]]*xd3a[[3, 4]]*xd3a[[3, 5]] - xd3a[[2, 1]]*xd3a[[3, 2]]*xd3a[[2, 4]]*xd3a[[3, 5]] + xd3a[[3, 1]]*xd3a[[3, 2]]*xd3a[[2, 4]]*xd3a[[2, 5]], xd3a[[2, 1]]^2*xd3a[[3, 4]]*xd3a[[3, 5]] - xd3a[[2, 1]]*xd3a[[3, 1]]*xd3a[[2, 4]]*xd3a[[3, 5]] + xd3a[[3, 1]]^2*xd3a[[2, 4]]*xd3a[[2, 5]], xd3a[[1, 1]]*xd3a[[2, 2]]*xd3a[[3, 4]]*xd3a[[3, 5]] - xd3a[[1, 1]]*xd3a[[3, 2]]*xd3a[[2, 4]]*xd3a[[3, 5]] + xd3a[[3, 1]]*xd3a[[3, 2]]*xd3a[[1, 4]]*xd3a[[2, 5]], xd3a[[1, 1]]*xd3a[[2, 1]]*xd3a[[3, 4]]*xd3a[[3, 5]] - xd3a[[1, 1]]*xd3a[[3, 1]]*xd3a[[2, 4]]*xd3a[[3, 5]] + xd3a[[3, 1]]^2*xd3a[[1, 4]]*xd3a[[2, 5]], -xd3a[[1, 1]]*xd3a[[2, 2]]*xd3a[[2, 3]]*xd3a[[3, 5]] + xd3a[[1, 1]]*xd3a[[3, 2]]*xd3a[[2, 3]]*xd3a[[2, 5]] + xd3a[[2, 1]]*xd3a[[2, 2]]*xd3a[[1, 3]]*xd3a[[3, 5]], -xd3a[[1, 1]]*xd3a[[2, 1]]*xd3a[[2, 3]]*xd3a[[3, 5]] + xd3a[[1, 1]]*xd3a[[3, 1]]*xd3a[[2, 3]]*xd3a[[2, 5]] + xd3a[[2, 1]]^2*xd3a[[1, 3]]*xd3a[[3, 5]], -xd3a[[1, 1]]*xd3a[[2, 2]]*xd3a[[2, 3]]*xd3a[[3, 5]] + xd3a[[2, 1]]*xd3a[[3, 2]]*xd3a[[1, 3]]*xd3a[[2, 5]], xd3a[[1, 4]]^2*xd3a[[3, 5]]^2 - xd3a[[1, 4]]*xd3a[[3, 4]]*xd3a[[1, 5]]*xd3a[[3, 5]] + xd3a[[3, 4]]^2*xd3a[[1, 5]]^2, xd3a[[1, 4]]^2*xd3a[[2, 5]]*xd3a[[3, 5]] - xd3a[[1, 4]]*xd3a[[2, 4]]*xd3a[[1, 5]]*xd3a[[3, 5]] + xd3a[[2, 4]]*xd3a[[3, 4]]*xd3a[[1, 5]]^2, xd3a[[1, 4]]^2*xd3a[[2, 5]]^2 - xd3a[[1, 4]]*xd3a[[2, 4]]*xd3a[[1, 5]]*xd3a[[2, 5]] + xd3a[[2, 4]]^2*xd3a[[1, 5]]^2, xd3a[[1, 1]]*xd3a[[1, 4]]*xd3a[[2, 5]]^2 - xd3a[[1, 1]]*xd3a[[2, 4]]*xd3a[[1, 5]]*xd3a[[2, 5]] + xd3a[[2, 1]]*xd3a[[2, 4]]*xd3a[[1, 5]]^2, xd3a[[1, 1]]^2*xd3a[[2, 5]]^2 - xd3a[[1, 1]]*xd3a[[2, 1]]*xd3a[[1, 5]]*xd3a[[2, 5]] + xd3a[[2, 1]]^2*xd3a[[1, 5]]^2, -xd3a[[1, 1]]*xd3a[[1, 4]]*xd3a[[3, 4]]*xd3a[[3, 5]] + xd3a[[1, 1]]*xd3a[[3, 4]]^2*xd3a[[1, 5]] + xd3a[[3, 1]]*xd3a[[1, 4]]^2*xd3a[[3, 5]], -xd3a[[1, 3]]*xd3a[[1, 4]]*xd3a[[2, 4]]*xd3a[[3, 5]] + xd3a[[1, 3]]*xd3a[[2, 4]]*xd3a[[3, 4]]*xd3a[[1, 5]] + xd3a[[2, 3]]*xd3a[[1, 4]]^2*xd3a[[3, 5]], -xd3a[[1, 1]]*xd3a[[1, 4]]*xd3a[[2, 4]]*xd3a[[3, 5]] + xd3a[[1, 1]]*xd3a[[2, 4]]*xd3a[[3, 4]]*xd3a[[1, 5]] + xd3a[[3, 1]]*xd3a[[1, 4]]^2*xd3a[[2, 5]], -xd3a[[1, 3]]*xd3a[[1, 4]]*xd3a[[2, 4]]*xd3a[[2, 5]] + xd3a[[1, 3]]*xd3a[[2, 4]]^2*xd3a[[1, 5]] + xd3a[[2, 3]]*xd3a[[1, 4]]^2*xd3a[[2, 5]], -xd3a[[1, 1]]*xd3a[[1, 3]]*xd3a[[2, 4]]*xd3a[[2, 5]] + xd3a[[1, 1]]*xd3a[[2, 3]]*xd3a[[1, 4]]*xd3a[[2, 5]] + xd3a[[2, 1]]*xd3a[[1, 3]]*xd3a[[2, 4]]*xd3a[[1, 5]], xd3a[[1, 1]]^2*xd3a[[2, 3]]*xd3a[[2, 5]] - xd3a[[1, 1]]*xd3a[[2, 1]]*xd3a[[1, 3]]*xd3a[[2, 5]] + xd3a[[2, 1]]^2*xd3a[[1, 3]]*xd3a[[1, 5]], -xd3a[[2, 2]]*xd3a[[2, 3]]*xd3a[[1, 4]]*xd3a[[3, 4]] + xd3a[[3, 2]]*xd3a[[1, 3]]*xd3a[[2, 4]]^2, xd3a[[2, 2]]^2*xd3a[[3, 4]]^2 - xd3a[[2, 2]]*xd3a[[3, 2]]*xd3a[[2, 4]]*xd3a[[3, 4]] + xd3a[[3, 2]]^2*xd3a[[2, 4]]^2, xd3a[[2, 1]]*xd3a[[2, 2]]*xd3a[[3, 4]]^2 - xd3a[[2, 1]]*xd3a[[3, 2]]*xd3a[[2, 4]]*xd3a[[3, 4]] + xd3a[[3, 1]]*xd3a[[3, 2]]*xd3a[[2, 4]]^2, xd3a[[2, 1]]^2*xd3a[[3, 4]]^2 - xd3a[[2, 1]]*xd3a[[3, 1]]*xd3a[[2, 4]]*xd3a[[3, 4]] + xd3a[[3, 1]]^2*xd3a[[2, 4]]^2, xd3a[[1, 1]]*xd3a[[2, 2]]*xd3a[[3, 4]]^2 - xd3a[[1, 1]]*xd3a[[3, 2]]*xd3a[[2, 4]]*xd3a[[3, 4]] + xd3a[[3, 1]]*xd3a[[3, 2]]*xd3a[[1, 4]]*xd3a[[2, 4]], xd3a[[1, 1]]*xd3a[[2, 1]]*xd3a[[3, 4]]^2 - xd3a[[1, 1]]*xd3a[[3, 1]]*xd3a[[2, 4]]*xd3a[[3, 4]] + xd3a[[3, 1]]^2*xd3a[[1, 4]]*xd3a[[2, 4]], -xd3a[[1, 1]]*xd3a[[2, 2]]*xd3a[[2, 3]]*xd3a[[3, 4]] + xd3a[[1, 1]]*xd3a[[3, 2]]*xd3a[[2, 3]]*xd3a[[2, 4]] + xd3a[[2, 1]]*xd3a[[2, 2]]*xd3a[[1, 3]]*xd3a[[3, 4]], -xd3a[[1, 1]]*xd3a[[2, 1]]*xd3a[[2, 3]]*xd3a[[3, 4]] + xd3a[[1, 1]]*xd3a[[3, 1]]*xd3a[[2, 3]]*xd3a[[2, 4]] + xd3a[[2, 1]]^2*xd3a[[1, 3]]*xd3a[[3, 4]], -xd3a[[1, 1]]*xd3a[[3, 2]]*xd3a[[2, 3]]*xd3a[[2, 4]] - xd3a[[2, 1]]*xd3a[[2, 2]]*xd3a[[1, 3]]*xd3a[[3, 4]] + xd3a[[2, 1]]*xd3a[[3, 2]]*xd3a[[1, 3]]*xd3a[[2, 4]], xd3a[[1, 3]]^2*xd3a[[2, 4]]^2 - xd3a[[1, 3]]*xd3a[[2, 3]]*xd3a[[1, 4]]*xd3a[[2, 4]] + xd3a[[2, 3]]^2*xd3a[[1, 4]]^2, xd3a[[1, 1]]^2*xd3a[[3, 4]]^2 - xd3a[[1, 1]]*xd3a[[3, 1]]*xd3a[[1, 4]]*xd3a[[3, 4]] + xd3a[[3, 1]]^2*xd3a[[1, 4]]^2, -xd3a[[1, 1]]*xd3a[[1, 3]]*xd3a[[2, 3]]*xd3a[[2, 4]] + xd3a[[1, 1]]*xd3a[[2, 3]]^2*xd3a[[1, 4]] + xd3a[[2, 1]]*xd3a[[1, 3]]^2*xd3a[[2, 4]], xd3a[[1, 1]]^2*xd3a[[2, 3]]^2 - xd3a[[1, 1]]*xd3a[[2, 1]]*xd3a[[1, 3]]*xd3a[[2, 3]] + xd3a[[2, 1]]^2*xd3a[[1, 3]]^2])


        @test MS3a[2] == ideal(R3b, [-xd3b[[2, 4]]*xd3b[[3, 5]] + xd3b[[3, 4]]*xd3b[[2, 5]], -xd3b[[1, 1]]*xd3b[[2, 4]] + xd3b[[2, 1]]*xd3b[[1, 4]], -xd3b[[1, 1]]*xd3b[[3, 2]] + xd3b[[3, 1]]*xd3b[[1, 2]], xd3b[[1, 1]]*xd3b[[2, 3]]*xd3b[[3, 5]] - xd3b[[2, 1]]*xd3b[[1, 3]]*xd3b[[3, 5]] + xd3b[[3, 1]]*xd3b[[1, 3]]*xd3b[[2, 5]], xd3b[[1, 1]]*xd3b[[2, 3]]*xd3b[[3, 4]] - xd3b[[2, 1]]*xd3b[[1, 3]]*xd3b[[3, 4]] + xd3b[[3, 1]]*xd3b[[1, 3]]*xd3b[[2, 4]], -xd3b[[1, 2]]*xd3b[[2, 3]]*xd3b[[3, 4]] - xd3b[[3, 2]]*xd3b[[1, 3]]*xd3b[[2, 4]] + xd3b[[3, 2]]*xd3b[[2, 3]]*xd3b[[1, 4]], -xd3b[[2, 1]]*xd3b[[1, 3]]*xd3b[[3, 4]] + xd3b[[3, 1]]*xd3b[[2, 3]]*xd3b[[1, 4]], xd3b[[2, 1]]^2*xd3b[[3, 5]]^2 - xd3b[[2, 1]]*xd3b[[3, 1]]*xd3b[[2, 5]]*xd3b[[3, 5]] + xd3b[[3, 1]]^2*xd3b[[2, 5]]^2, xd3b[[2, 1]]^2*xd3b[[3, 4]]*xd3b[[3, 5]] - xd3b[[2, 1]]*xd3b[[3, 1]]*xd3b[[2, 4]]*xd3b[[3, 5]] + xd3b[[3, 1]]^2*xd3b[[2, 4]]*xd3b[[2, 5]], -xd3b[[1, 2]]*xd3b[[1, 3]]*xd3b[[2, 4]]*xd3b[[3, 5]] + xd3b[[1, 2]]*xd3b[[2, 3]]*xd3b[[1, 4]]*xd3b[[3, 5]] + xd3b[[3, 2]]*xd3b[[1, 3]]*xd3b[[1, 4]]*xd3b[[2, 5]], xd3b[[2, 1]]*xd3b[[1, 2]]*xd3b[[3, 4]]*xd3b[[3, 5]] - xd3b[[2, 1]]*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[3, 5]] + xd3b[[3, 1]]*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[2, 5]], xd3b[[1, 1]]*xd3b[[2, 1]]*xd3b[[3, 4]]*xd3b[[3, 5]] - xd3b[[2, 1]]*xd3b[[3, 1]]*xd3b[[1, 4]]*xd3b[[3, 5]] + xd3b[[3, 1]]^2*xd3b[[1, 4]]*xd3b[[2, 5]], xd3b[[1, 1]]*xd3b[[3, 1]]*xd3b[[2, 3]]*xd3b[[2, 5]] - xd3b[[2, 1]]^2*xd3b[[1, 3]]*xd3b[[3, 5]], xd3b[[2, 1]]*xd3b[[1, 2]]*xd3b[[2, 3]]*xd3b[[3, 5]] + xd3b[[2, 1]]*xd3b[[3, 2]]*xd3b[[1, 3]]*xd3b[[2, 5]] - xd3b[[3, 1]]*xd3b[[1, 2]]*xd3b[[2, 3]]*xd3b[[2, 5]], xd3b[[1, 1]]*xd3b[[1, 2]]*xd3b[[2, 3]]*xd3b[[3, 5]] + xd3b[[1, 1]]*xd3b[[3, 2]]*xd3b[[1, 3]]*xd3b[[2, 5]] - xd3b[[2, 1]]*xd3b[[1, 2]]*xd3b[[1, 3]]*xd3b[[3, 5]], xd3b[[2, 1]]^2*xd3b[[3, 4]]^2 - xd3b[[2, 1]]*xd3b[[3, 1]]*xd3b[[2, 4]]*xd3b[[3, 4]] + xd3b[[3, 1]]^2*xd3b[[2, 4]]^2, -xd3b[[1, 2]]*xd3b[[1, 3]]*xd3b[[2, 4]]*xd3b[[3, 4]] + xd3b[[1, 2]]*xd3b[[2, 3]]*xd3b[[1, 4]]*xd3b[[3, 4]] + xd3b[[3, 2]]*xd3b[[1, 3]]*xd3b[[1, 4]]*xd3b[[2, 4]], xd3b[[2, 1]]*xd3b[[1, 2]]*xd3b[[3, 4]]^2 - xd3b[[2, 1]]*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[3, 4]] + xd3b[[3, 1]]*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[2, 4]], xd3b[[1, 1]]*xd3b[[2, 1]]*xd3b[[3, 4]]^2 - xd3b[[1, 1]]*xd3b[[3, 1]]*xd3b[[2, 4]]*xd3b[[3, 4]] + xd3b[[3, 1]]^2*xd3b[[1, 4]]*xd3b[[2, 4]], xd3b[[1, 1]]*xd3b[[3, 1]]*xd3b[[2, 3]]*xd3b[[2, 4]] - xd3b[[2, 1]]^2*xd3b[[1, 3]]*xd3b[[3, 4]], -xd3b[[1, 1]]*xd3b[[3, 2]]*xd3b[[2, 3]]*xd3b[[2, 4]] + xd3b[[2, 1]]*xd3b[[1, 2]]*xd3b[[2, 3]]*xd3b[[3, 4]] + xd3b[[2, 1]]*xd3b[[3, 2]]*xd3b[[1, 3]]*xd3b[[2, 4]], xd3b[[1, 1]]*xd3b[[1, 2]]*xd3b[[2, 3]]*xd3b[[3, 4]] + xd3b[[1, 1]]*xd3b[[3, 2]]*xd3b[[1, 3]]*xd3b[[2, 4]] - xd3b[[2, 1]]*xd3b[[1, 2]]*xd3b[[1, 3]]*xd3b[[3, 4]], xd3b[[1, 3]]^2*xd3b[[2, 4]]^2 - xd3b[[1, 3]]*xd3b[[2, 3]]*xd3b[[1, 4]]*xd3b[[2, 4]] + xd3b[[2, 3]]^2*xd3b[[1, 4]]^2, xd3b[[1, 2]]^2*xd3b[[3, 4]]^2 - xd3b[[1, 2]]*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[3, 4]] + xd3b[[3, 2]]^2*xd3b[[1, 4]]^2, xd3b[[1, 1]]*xd3b[[1, 2]]*xd3b[[3, 4]]^2 - xd3b[[1, 1]]*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[3, 4]] + xd3b[[3, 1]]*xd3b[[3, 2]]*xd3b[[1, 4]]^2, xd3b[[1, 1]]^2*xd3b[[3, 4]]^2 - xd3b[[1, 1]]*xd3b[[3, 1]]*xd3b[[1, 4]]*xd3b[[3, 4]] + xd3b[[3, 1]]^2*xd3b[[1, 4]]^2, -xd3b[[1, 1]]*xd3b[[1, 3]]*xd3b[[2, 3]]*xd3b[[2, 4]] + xd3b[[1, 1]]*xd3b[[2, 3]]^2*xd3b[[1, 4]] + xd3b[[2, 1]]*xd3b[[1, 3]]^2*xd3b[[2, 4]], xd3b[[1, 1]]^2*xd3b[[2, 3]]^2 - xd3b[[1, 1]]*xd3b[[2, 1]]*xd3b[[1, 3]]*xd3b[[2, 3]] + xd3b[[2, 1]]^2*xd3b[[1, 3]]^2, xd3b[[1, 1]]*xd3b[[3, 1]]*xd3b[[3, 2]]*xd3b[[2, 5]]^2 + xd3b[[2, 1]]^2*xd3b[[1, 2]]*xd3b[[3, 5]]^2 - xd3b[[2, 1]]*xd3b[[3, 1]]*xd3b[[1, 2]]*xd3b[[2, 5]]*xd3b[[3, 5]], -xd3b[[1, 2]]*xd3b[[1, 3]]*xd3b[[2, 3]]*xd3b[[2, 4]]*xd3b[[3, 5]] + xd3b[[1, 2]]*xd3b[[1, 3]]*xd3b[[2, 3]]*xd3b[[3, 4]]*xd3b[[2, 5]] + xd3b[[1, 2]]*xd3b[[2, 3]]^2*xd3b[[1, 4]]*xd3b[[3, 5]] + xd3b[[3, 2]]*xd3b[[1, 3]]^2*xd3b[[2, 4]]*xd3b[[2, 5]], xd3b[[1, 1]]*xd3b[[3, 1]]*xd3b[[3, 2]]*xd3b[[2, 4]]*xd3b[[2, 5]] + xd3b[[2, 1]]^2*xd3b[[1, 2]]*xd3b[[3, 4]]*xd3b[[3, 5]] - xd3b[[2, 1]]^2*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[3, 5]], xd3b[[1, 1]]*xd3b[[3, 2]]^2*xd3b[[1, 4]]*xd3b[[2, 5]] + xd3b[[2, 1]]*xd3b[[1, 2]]^2*xd3b[[3, 4]]*xd3b[[3, 5]] - xd3b[[2, 1]]*xd3b[[1, 2]]*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[3, 5]], xd3b[[1, 1]]^2*xd3b[[3, 2]]*xd3b[[2, 3]]*xd3b[[2, 5]] - xd3b[[2, 1]]^2*xd3b[[1, 2]]*xd3b[[1, 3]]*xd3b[[3, 5]], xd3b[[1, 2]]*xd3b[[2, 3]]^2*xd3b[[1, 4]]*xd3b[[3, 4]] + xd3b[[3, 2]]*xd3b[[1, 3]]^2*xd3b[[2, 4]]^2, xd3b[[1, 1]]*xd3b[[3, 1]]*xd3b[[3, 2]]*xd3b[[2, 4]]^2 + xd3b[[2, 1]]^2*xd3b[[1, 2]]*xd3b[[3, 4]]^2 - xd3b[[2, 1]]^2*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[3, 4]], xd3b[[1, 1]]*xd3b[[3, 2]]^2*xd3b[[1, 4]]*xd3b[[2, 4]] + xd3b[[2, 1]]*xd3b[[1, 2]]^2*xd3b[[3, 4]]^2 - xd3b[[2, 1]]*xd3b[[1, 2]]*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[3, 4]], xd3b[[1, 1]]^2*xd3b[[3, 2]]*xd3b[[2, 3]]*xd3b[[2, 4]] - xd3b[[2, 1]]^2*xd3b[[1, 2]]*xd3b[[1, 3]]*xd3b[[3, 4]], xd3b[[1, 2]]^2*xd3b[[2, 3]]^2*xd3b[[3, 5]]^2 + xd3b[[1, 2]]*xd3b[[3, 2]]*xd3b[[1, 3]]*xd3b[[2, 3]]*xd3b[[2, 5]]*xd3b[[3, 5]] + xd3b[[3, 2]]^2*xd3b[[1, 3]]^2*xd3b[[2, 5]]^2, xd3b[[1, 1]]^2*xd3b[[3, 2]]^2*xd3b[[2, 5]]^2 + xd3b[[2, 1]]^2*xd3b[[1, 2]]^2*xd3b[[3, 5]]^2 - xd3b[[2, 1]]*xd3b[[3, 1]]*xd3b[[1, 2]]^2*xd3b[[2, 5]]*xd3b[[3, 5]], xd3b[[1, 1]]^2*xd3b[[3, 2]]^2*xd3b[[2, 4]]*xd3b[[2, 5]] + xd3b[[2, 1]]^2*xd3b[[1, 2]]^2*xd3b[[3, 4]]*xd3b[[3, 5]] - xd3b[[2, 1]]^2*xd3b[[1, 2]]*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[3, 5]], xd3b[[1, 1]]^2*xd3b[[3, 2]]^2*xd3b[[2, 4]]^2 + xd3b[[2, 1]]^2*xd3b[[1, 2]]^2*xd3b[[3, 4]]^2 - xd3b[[2, 1]]^2*xd3b[[1, 2]]*xd3b[[3, 2]]*xd3b[[1, 4]]*xd3b[[3, 4]]])


        @test MS4a[2] == ideal(R4a, [])

        @test MS4b[2] == ideal(R4b, [-xd4[[4, 1]]*xd4[[5, 2]] + xd4[[5, 1]]*xd4[[4, 2]]])
                
    end

    
end

