using Oscar.Oscardb.Mongoc

@testset verbose=true "oscarDB" begin
    # conditions from which the test database dump was generated.
    add_constraints_poly = ["DIM" => Dict("\$lte" => 3)]
    function _acp(a::Array)
        return append!(Array{Pair{String,Any}}(a), add_constraints_poly)
    end
    add_constraints_mat = [:("N_ELEMENTS" <= 4)]
    function _acm(a::Array)
        return append!(Array{Expr}(a), add_constraints_mat)
    end

    @testset verbose=true "Basic functionality" begin
        # Types
        @test Oscar.Oscardb.get_db() isa Oscar.Oscardb.Database
        db = Oscar.Oscardb.get_db()
        @test db["Surfaces"] isa Oscar.Oscardb.Collection
        try
            @test Mongoc.ping(db.mdb.client)["ok"] == 1
        catch
            @test "not" == "connected"
        end
        collection_bo = db["Surfaces"]
        query = Dict("_id" => "abelian_d10_pi6")

        @test length(collection_bo) == 48
        @test length(collection_bo, query) == 1

        @testset verbose=true "Iterator (Collection)" begin
            @test iterate(collection_bo) isa Tuple{Mongoc.BSON, Mongoc.Cursor{Mongoc.Collection}}
        end
    end
end
