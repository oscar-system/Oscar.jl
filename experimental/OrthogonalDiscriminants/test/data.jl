@testset "all_od_infos" begin
  iob = IOBuffer()
  show_OD_info("G2(3)", iob)
  str1 = String(take!(iob))

  all_entries = all_od_infos();
  @test all_entries isa Vector

  show_OD_info("G2(3)", iob)
  str2 = String(take!(iob))
  @test str1 == str2

  @test length(all_entries) == length(all_od_infos(is_simple)) +
                               length(all_od_infos(! is_simple))
  @test length(all_entries) == length(all_od_infos(is_sporadic_simple)) +
                               length(all_od_infos(! is_sporadic_simple))
  @test length(all_od_infos(is_sporadic_simple, is_simple)) ==
                               length(all_od_infos(is_sporadic_simple))
  @test length(all_od_infos(is_sporadic_simple, ! is_simple)) == 0
  A5_entries = all_od_infos(identifier => "A5")
  @test length(A5_entries) == 3
  @test length(all_od_infos(identifier => ["A5", "A6"])) ==
        length(A5_entries) +
        length(all_od_infos(identifier => "A6"))
  A_entries = all_od_infos(identifier => x -> x[1] == 'A')
  @test A5_entries == A_entries[1:3]
  @test A5_entries == all_od_infos(is_simple, identifier => "A5")
  @test all_od_infos(identifier => ["A5", "A6"], characteristic => [3, 5]) ==
        all_od_infos(identifier => ["A5", "A6"], characteristic => isodd)
  @test all_od_infos(characteristic => [0, 2]) ==
        all_od_infos(characteristic => iseven)
  @test length(all_od_infos(dim => isodd)) == 0
  A_dims = collect(Set([parse(Int, filter(isdigit, x[:charname])) for x in A_entries]))
  @test all_od_infos(identifier => x -> x[1] == 'A', dim => A_dims) == A_entries
  Oplusminus = all_od_infos(orthogonal_discriminant => ["O+", "O-"]);
  @test length(all_od_infos(orthogonal_discriminant => "O+")) +
        length(all_od_infos(orthogonal_discriminant => "O-")) ==
        length(Oplusminus)
  @test length(Oplusminus) <= length(all_od_infos(characteristic => is_positive))
  deg_1 = all_od_infos(degree => 1);
  deg_other = all_od_infos(degree => (x -> x > 1));
  @test length(all_entries) == length(deg_1) + length(deg_other)
  @test deg_1 == all_od_infos(degree => [1]);
end

@testset "orthogonal_discriminant(s)" begin
  for name in ["S4", "A5", "A6", "J1"]
    t = character_table(name)
    vals = orthogonal_discriminants(t)
    @test vals == [orthogonal_discriminant(chi) for chi in t]
    @test length(vals) == length(t)
    @test vals[1] == ""
    for p in vcat([i for (i, e) in factor(order(t))], [101])
      modt = mod(t, p)
      vals = orthogonal_discriminants(modt)
      @test vals == [orthogonal_discriminant(chi) for chi in modt]
      @test length(vals) == length(modt)
      @test vals[1] == ""
    end
  end
end

@testset "orthogonal_discriminants for central extensions" begin
  name = "2.L3(4)"
  t = character_table(name)
  @test orthogonal_discriminants(t) ==
  ["","21","","","","","","","","105","","","1","1","1","1","-1","-1"]
end


# Make sure that the tests start without cached character tables.
# (Running the tests will store information that changes some test outputs,
# thus running the tests twice needs these calls.)
GAP.Globals.UnloadCharacterTableData()
empty!(Oscar.character_tables_by_id)

Oscar.@_AuxDocTest "show and print character tables", (fix = false),
raw"""
a table with irrational ODs

```jldoctest show_with_ODs.test
julia> using Oscar

julia> t = character_table("L2(11)");

julia> Oscar.OrthogonalDiscriminants.show_with_ODs(t, stdout)
L2(11)

                2  2  2  1  .  .  1   .   .
                3  1  1  1  .  .  1   .   .
                5  1  .  .  1  1  .   .   .
               11  1  .  .  .  .  .   1   1
                                           
                  1a 2a 3a 5a 5b 6a 11a 11b
               2P 1a 1a 3a 5b 5a 3a 11b 11a
               3P 1a 2a 1a 5b 5a 2a 11a 11b
               5P 1a 2a 3a 1a 1a 6a 11a 11b
              11P 1a 2a 3a 5a 5b 6a  1a  1a
    d      OD   2                          
X_1 1           +  1  1  1  1  1  1   1   1
X_2 2           o  5  1 -1  .  .  1   B  /B
X_3 2           o  5  1 -1  .  .  1  /B   B
X_4 1     -11   + 10 -2  1  .  .  1  -1  -1
X_5 1     -11   + 10  2  1  .  . -1  -1  -1
X_6 1           + 11 -1 -1  1  1 -1   .   .
X_7 2 22-11b5   + 12  .  .  A A*  .   1   1
X_8 2 11b5+33   + 12  .  . A*  A  .   1   1

A = -z_5^3 - z_5^2 - 1
A* = z_5^3 + z_5^2
B = z_11^9 + z_11^5 + z_11^4 + z_11^3 + z_11
/B = -z_11^9 - z_11^5 - z_11^4 - z_11^3 - z_11 - 1
```

a table with unknown ODs (outside the scope of the database)

```jldoctest show_with_ODs.test
julia> t = character_table("D8");

julia> Oscar.OrthogonalDiscriminants.show_with_ODs(t, stdout)
D8

          2  3  3  2  2  2
                          
            1a 2a 4a 2b 2c
         2P 1a 1a 2a 1a 1a
    d OD  2               
X_1 1     +  1  1  1  1  1
X_2 1     +  1  1  1 -1 -1
X_3 1     +  1  1 -1  1 -1
X_4 1     +  1  1 -1 -1  1
X_5 1  ?  +  2 -2  .  .  .
```
"""
