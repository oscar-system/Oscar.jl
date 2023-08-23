@testset "all_od_infos" begin
  all_entries = all_od_infos();
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
  @test length(Oplusminus) <= length(all_od_infos(characteristic => ispositive))
end
