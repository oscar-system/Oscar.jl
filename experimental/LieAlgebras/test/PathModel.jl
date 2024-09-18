@testset "LieAlgebras.PathModel" begin
  GAP.Packages.load("QuaGroup")

  @testset "StraightLinePathModel" begin
    #=@testset "falpha" begin
      R = root_system(:A, 2)
      P = path_model(R, [1, 1])
      dom = dominant_path(P)

      @test falpha(dom, 1) !== dom
      @test falpha(dom, 1, 1) !== dom
      @test falpha(dom, [1], [1]) !== dom
    end=#

    @testset "(P::LSPathModel)(wt::WeightLatticeElem)" begin
      R = root_system(:A, 3)
      P = straight_line_path_model(R, [2, 2, 2])
      L = lie_algebra(QQ, :A, 3)

      char = dominant_character(L, [2, 2, 2])
      for k in keys(char)
        @test char[k] == length(P(k))
      end
    end

    @test "falpha" begin
      # test that everything is copied
    end

    #=@testset "falpha!" begin
      R = root_system(:A, 3)
      P = straight_line_path_model(R, [2, 2, 2])
      g = dominant_path(P)

      gR = GAP.Globals.RootSystem(GAP.Obj("A"), 3)
      gp = GAP.Globals.DominantLSPath(gR, GAP.Obj([2, 2, 2]))

      w0 = longest_element(weyl_group(R))
      path = zeros(length(w0))

      n = 1
      while true
        ii = Int(w0[n])
        while Oscar.LieAlgebras.phi(g, ii) > 0
          path[n] += 1
          falpha!(g, ii)
          gp = GAP.Globals.Falpha(gp, ii)
          @test gp != GAP.Globals.fail

          ls = GAP.Globals.LSSequence(gp)
          for (j, s) in Iterators.enumerate(Iterators.accumulate(+, g.d))
            @test QQ(ls[2][j+1]) == s # GAP starts at 0
          end
          for (j, w) in Iterators.enumerate(g.s)
            wt = w*P.wt
            for k in 1:length(ls[1][j])
              @test wt[k] == ls[1][j][k]
            end
          end
        end

        @test iszero(falpha(g, ii))
        @test GAP.Globals.Falpha(gp, ii) == GAP.Globals.fail

        if n < length(w0)
          n += 1
        else
          n -= 1
        end
      end
    end =#
  end
end
