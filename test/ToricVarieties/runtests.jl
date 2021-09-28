const pm = Polymake

using Oscar

# This enables a debug flag that will produce some more output
# during polymake's wrapper compilation (at run-time).
# We want to keep this active for CI to make it easier to debug future
# compilation errors, especially when they only appear during CI.
if (haskey(ENV, "GITHUB_ACTIONS"))
    Polymake.shell_execute(raw"$Verbose::cpp=3;")
end

@testset "Toric Varieties" begin
    
    square = cube(2)
    ntv_square = NormalToricVariety(square)

    @testset "core functionality" begin
        @test isnormal(ntv_square)
        @test iscomplete(ntv_square)
        @test isprojective(ntv_square)
        @test !isaffine(ntv_square)
        @test issimplicial(ntv_square)
        @test issmooth(ntv_square)
    end

end
