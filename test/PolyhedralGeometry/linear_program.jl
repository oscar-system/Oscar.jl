NF, sr2 = quadratic_field(2)
ENF, sre2 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
Qx, x = QQ["x"]
K, (a1, a2) = embedded_number_field([x^2 - 2, x^3 - 5], [(0, 2), (0, 2)])

for f in (QQ, ENF)

    T = elem_type(f)

    @testset "(MixedInteger)LinearProgram{$T}" begin

        pts = [1 0; 0 0; 0 1]
        Q0 = convex_hull(f, pts)
        Q1 = convex_hull(f, pts, [1 1])
        Q2 = convex_hull(f, pts, [1 1], [1 1])
        square = cube(f, 2)
        C1 = cube(f, 2, 0, 1)
        Pos = polyhedron(f, [-1 0 0; 0 -1 0; 0 0 -1], [0,0,0])
        L = polyhedron(f, [-1 0 0; 0 -1 0], [0,0])
        point = convex_hull(f, [0 1 0])
        # this is to make sure the order of some matrices below doesn't change
        Polymake.prefer("beneath_beyond") do
            affine_hull(point)
        end
        s = simplex(f, 2)
        rsquare = cube(f, 2, QQFieldElem(-3,2), QQFieldElem(3,2))

        @testset "linear programs" begin
            LP1 = linear_program(square,[1,3])
            LP2 = linear_program(square,[2,2]; k=3, convention = :min)
            LP3 = linear_program(Pos,[1,2,3])
            @test LP1 isa LinearProgram{T}
            @test LP2 isa LinearProgram{T}
            @test LP3 isa LinearProgram{T}
            
            @test solve_lp(LP1)==(4,[1,1])
            @test solve_lp(LP2)==(-1,[-1,-1])
            if T == QQFieldElem
                @test string(solve_lp(LP3))==string("(inf, nothing)")
            else
                @test string(solve_lp(LP3))==string("((inf), nothing)")
            end
        end
        
        if T == QQFieldElem
            
            @testset "mixed integer linear programs" begin
                MILP1 = mixed_integer_linear_program(rsquare, [1,3], integer_variables=[1])
                MILP2 = mixed_integer_linear_program(rsquare, [2,2]; k=3, convention = :min)
                MILP3 = mixed_integer_linear_program(Pos, [1,2,3])
                @test MILP1 isa MixedIntegerLinearProgram{T}
                @test MILP2 isa MixedIntegerLinearProgram{T}
                @test MILP3 isa MixedIntegerLinearProgram{T}
                
                @test solve_milp(MILP1)==(11//2,[1,3//2])
                @test solve_milp(MILP2)==(-1,[-1,-1])
                @test string(solve_milp(MILP3))==string("(inf, nothing)")
            end
            
            @testset "LinearProgram: lp and mps files" begin
                LP1 = linear_program(square,[1,3])
                LP2 = linear_program(rsquare,[2,2]; k=3, convention = :min)
                
                buffer = IOBuffer()
                @test save_lp(buffer, LP2) === nothing
                @test String(take!(buffer)) ==
                    """
                      MINIMIZE
                        obj: +2 x1 +2 x2 +3
                      Subject To
                        ie0: +2 x1 >= -3
                        ie1: -2 x1 >= -3
                        ie2: +2 x2 >= -3
                        ie3: -2 x2 >= -3
                      BOUNDS
                        x1 free
                        x2 free
                      END
                      """
                MILP1 = mixed_integer_linear_program(rsquare,[1,3], integer_variables=[1])
                MILP2 = mixed_integer_linear_program(rsquare,[2,2]; k=3, convention = :max)
                
                @test save_mps(buffer, MILP1) === nothing
                @test String(take!(buffer)) ==
                    """
                      * Class:	MIP
                      * Rows:		5
                      * Columns:	2
                      * Format:	MPS
                      *
                      NAME          unnamed#2
                      ROWS
                       N  C0000000
                       G  R0000000
                       G  R0000001
                       G  R0000002
                       G  R0000003
                      COLUMNS
                          M0000000  'MARKER'                 'INTORG'
                          x1        C0000000  1                        R0000000  1
                          x1        R0000001  -1                       
                          M0000000  'MARKER'                 'INTEND'
                          x2        C0000000  3                        R0000002  1
                          x2        R0000003  -1                       
                      RHS
                          B         R0000000  -1.5                     R0000001  -1.5
                          B         R0000002  -1.5                     R0000003  -1.5
                      BOUNDS
                       FR BND       x1  
                       FR BND       x2  
                      ENDATA
                      """
                
                for lp in (LP1, LP2, MILP1, MILP2)
                    mktempdir() do path
                        @test save_lp(joinpath(path,"lp.lp"), lp) === nothing
                        loaded = load_lp(joinpath(path,"lp.lp"))
                        @test typeof(loaded) == typeof(lp)
                        @test feasible_region(lp) == feasible_region(loaded)
                        @test objective_function(lp) == objective_function(loaded)
                        @test optimal_value(lp) == optimal_value(loaded)
                        if lp isa MixedIntegerLinearProgram
                            @test optimal_solution(lp) == optimal_solution(loaded)
                        else
                            @test optimal_vertex(lp) == optimal_vertex(loaded)
                        end
                        
                        @test save_mps(joinpath(path,"lp.mps"), lp) === nothing
                        loaded = load_mps(joinpath(path,"lp.mps"))
                        @test typeof(loaded) == typeof(lp)
                        @test feasible_region(lp) == feasible_region(loaded)
                        @test objective_function(lp) == objective_function(loaded)
                        if lp.convention === :max
                            # mps file don't store max / min
                            @test optimal_value(lp) == optimal_value(loaded)
                            if lp isa MixedIntegerLinearProgram
                                @test optimal_solution(lp) == optimal_solution(loaded)
                            else
                                @test optimal_vertex(lp) == optimal_vertex(loaded)
                            end
                        end
                    end
                end
            end
        end
    end
end
