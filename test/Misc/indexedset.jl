@testset "Basic IndexedSet usage" begin
  # Construction
  items = ["ball", "flower", "house"]
  s = IndexedSet(items)
  @test !isempty(s)
  @test collect(s) == items
  @test length(s) == 3

  # Membership and indexing.
  @test "ball" in s
  @test s[s("flower")] == "flower"

  # Missing elements.
  @test !("stone" in s)
  @test s("stone") == 0

  # Pushing new and redundant elements.
  push!(s, "ball", "stone")
  @test collect(s) == [items; "stone"]

  # Appending new items.
  items_2 = ["garden", "rock"]
  append!(s, items_2)
  @test collect(s) == [items; "stone"; items_2]
end
