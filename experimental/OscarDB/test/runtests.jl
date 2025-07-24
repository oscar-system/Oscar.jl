using Oscar.OscarDB.Mongoc

@testset "OscarDB" begin
  db = Oscar.OscarDB.get_db()
  collection_tsc = db["TransitiveSimplicialComplexes"]
  
  @testset "Types" begin
    @test Oscar.OscarDB.get_db() isa Oscar.OscarDB.Database
    @test db["TransitiveSimplicialComplexes"] isa Oscar.OscarDB.Collection
    @test collection_tsc isa Oscar.OscarDB.Collection
  end

  @testset "Querying" begin
    tsc = Oscar.OscarDB.find_one(db["TransitiveSimplicialComplexes"],
                                 Dict("data.betti_numbers" => ["0", "1", "0"]))
    @test betti_numbers(tsc) == [0, 1, 0]
  end
end
