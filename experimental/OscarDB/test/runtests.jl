using Oscar.OscarDB.Mongoc

@testset verbose=true "OscarDB" begin

  @testset verbose=true "Basic functionality" begin
    # Types
    @test Oscar.OscarDB.get_db() isa Oscar.OscarDB.Database
    db = Oscar.OscarDB.get_db()
    try
      @test Mongoc.ping(db.mdb.client)["ok"] == 1
    catch
      @test "not" == "connected"
    end
  end

  @testset "Collections" begin
    # can only test collections that are in the .github/oscardb_dump/oscar
    db = Oscar.OscarDB.get_db()
    @testset "TransitiveSimplicialComplexes" begin
      collection_tsc = db["TransitiveSimplicialComplexes"]

      @testset verbose=true "Iterator (Collection)" begin
        @test iterate(collection_tsc) isa Tuple{Mongoc.BSON, Mongoc.Cursor{Mongoc.Collection}}
      end

      @testset "Types" begin
        @test Oscar.OscarDB.get_db() isa Oscar.OscarDB.Database
        @test db["TransitiveSimplicialComplexes"] isa Oscar.OscarDB.Collection
        @test collection_tsc isa Oscar.OscarDB.Collection
      end

      @testset "Querying" begin
        tsc = Oscar.OscarDB.find_one(db["TransitiveSimplicialComplexes"],
                                     Dict("data.betti_numbers" => ["0", "0", "0", "1"]))
        @test betti_numbers(tsc) == [0, 0, 0, 1]
      end
    end

    @testset "AlgebraicStatistics.SmallTreeModels" begin
      collection_stm = db["AlgebraicStatistics.SmallTreeModels"]

      @testset verbose=true "Iterator (Collection)" begin
        @test iterate(collection_stm) isa Tuple{Mongoc.BSON, Mongoc.Cursor{Mongoc.Collection}}
      end

      @testset "Types" begin
        @test Oscar.OscarDB.get_db() isa Oscar.OscarDB.Database
        @test db["AlgebraicStatistics.SmallTreeModels"] isa Oscar.OscarDB.Collection
        @test collection_stm isa Oscar.OscarDB.Collection
      end

      @testset "Querying" begin
        @test length(collection_stm, Dict("data.model_type" => "JC")) == 6
      end
    end
  end

  # here to avoid
  @testset "Doctests" begin

  end
end
