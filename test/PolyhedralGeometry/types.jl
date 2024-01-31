@testset "types" begin
    
  @testset "IncidenceMatrix" begin
        
    im = IncidenceMatrix([[1,2,3],[4,5,6]])
    @test nrows(im) == 2
    @test ncols(im) == 6
    @test row(im, 1) isa Set{Int}
    @test row(im, 1) == Set{Int}([1, 2, 3])
    @test column(im, 2) isa Set{Int}
    @test column(im, 2) == Set{Int}([1])
        
  end

  a = [1, 2, 3]
  b = [8, 6, 4]

  (ENF, _) = _prepare_scalar_types()[2]

  @testset "$T" for (T, fun) in ((PointVector, point_vector), (RayVector, ray_vector))

    @test fun(a) isa T{QQFieldElem}

    for f in (ZZ, QQ, ENF)

      U = elem_type(f)
      @testset "$T{$U}" begin

        @test T{U} <: AbstractVector
        @test T{U} <: AbstractVector{U}

        @test fun(f, a) isa T
        @test fun(f, a) isa T{U}

        @test fun(f, 7) == f.(zeros(Int, 7))

        A = fun(f, a)

        @test A[1] isa U
        @test A[2] == 2
        @test A[begin] == 1
        @test A[end] == 3
        A[3] = f(7)
        @test A[end] == 7
        A[3] = f(3)

        @test size(A) == (3,)

        @test_throws BoundsError A[0]
        @test_throws BoundsError A[4]

        for g in (ZZ, QQ, ENF)
          V = elem_type(g)

          B = fun(g, b)

          for op in [+, -]
            @test op(A, B) isa T

            @test op(A, B) == op(a, b)
          end

          @test *(g(3), A) isa T

          @test *(g(3), A) == 3 * a

        end

        for op in [+, -]
          @test op(A, b) isa T
          @test op(A, b) isa T{U}
          @test op(A, b) == op(a, b)
        end

        @test 3 * A isa T
        @test 3 * A isa T{U}
        @test 3 * A == 3 * a

        if f != ENF
          let h = Int
            Ah = h.(A)
            @test Ah isa Vector{Int}
            @test Ah == [1, 2, 3]
          end
          
          let h = ENF
            Ah = h.(A)
            @test Ah isa T{elem_type(ENF)}
            @test Ah == [1, 2, 3]
          end
        end

      end
    end

  end





  @testset "$T" for (T, f) in ((AffineHalfspace, affine_halfspace), (AffineHyperplane, affine_hyperplane))
        
    for p in [QQ, ENF]
      U = elem_type(p)
      @test f(p, a, 0) isa T{U}
      @test f(p, permutedims(a), 0) isa T{U}

      @test f(p, a, 0) == f(p, permutedims(a), 0) == f(p, a) == f(p, permutedims(a))

      A = f(p, a, 0)
      B = f(p, b, 2)

      @test A != B

      @test normal_vector(A) isa Vector{U}
      @test normal_vector(A) == a
      @test negbias(A) isa U
      @test negbias(A) == 0
            
            
      @test normal_vector(B) == b
      @test negbias(B) == 2
    end
        
    @test f(a, 0) isa T{QQFieldElem}
        
  end
    
  @testset "$T" for (T, f) in ((LinearHalfspace, linear_halfspace), (LinearHyperplane, linear_hyperplane))
        
    for p in [QQ, ENF]
      U = elem_type(p)
      @test f(p, a) isa T{U}
      @test f(p, permutedims(a)) isa T{U}

      @test f(p, a) == f(p, permutedims(a))

      A = f(p, a)
      B = f(p, b)

      @test A != B

      @test normal_vector(A) isa Vector{U}
      @test normal_vector(A) == a
      @test negbias(A) isa U
      @test negbias(A) == 0
            
            
      @test normal_vector(B) == b
      @test negbias(B) == 0
    end
        
    @test f(a) isa T{QQFieldElem}
        
  end
    
  for p in [QQ, ENF]
    U = elem_type(p)
    let A = linear_halfspace(p, a)
      @test invert(A) isa LinearHalfspace{U}
      Ai = invert(A)
      @test normal_vector(Ai) == -a
    end
    let A = affine_halfspace(p, a, 8)
      @test invert(A) isa AffineHalfspace{U}
      Ai = invert(A)
      @test normal_vector(Ai) == -a
      @test negbias(Ai) == -8
    end
  end
    
  @test halfspace(a) isa LinearHalfspace{QQFieldElem}
  @test hyperplane(a) isa LinearHyperplane{QQFieldElem}
  @test halfspace(a, 0) isa AffineHalfspace{QQFieldElem}
  @test hyperplane(a, 0) isa AffineHyperplane{QQFieldElem}

end
