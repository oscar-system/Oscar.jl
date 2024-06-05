@testset "constructors for polycyclic groups" begin
  # Given a GAP group, `PcGroup` just wraps the given full pc group
  # or signals an error.
  G = small_group(12, 2)
  PCG = PcGroup(GapObj(G))
  @test PCG isa PcGroup
  @test GapObj(PCG) === GapObj(G)

  G = derived_subgroup(G)[1]
  @test_throws AssertionError PcGroup(GapObj(G))

  G = symmetric_group(4)
  @test_throws AssertionError PcGroup(GapObj(G))

  # Given a GAP group `G`, `pc_group` may create a new full pc group
  # if `G` is a proper subgroup or if its generators are not the defining ones.
  G = small_group(12, 2)
  PCG = pc_group(GapObj(G))
  @test PCG isa PcGroup
  @test GapObj(PCG) === GapObj(G)

  G = derived_subgroup(G)[1]
  PCG = pc_group(GapObj(G))
  @test PCG isa PcGroup
  @test GapObj(PCG) !== GapObj(G)

  G = symmetric_group(4)
  @test_throws ArgumentError pc_group(GapObj(G)) isa PcGroup
end

@testset "create polycyclic groups from collectors" begin

  # finite polycyclic groups
  c = collector(2, Int);
  set_relative_order!(c, 1, 2)
  set_relative_order!(c, 2, 3)
  set_power!(c, 1, [2 => 1])
  gg = pc_group(c)
  @test describe(gg) == "C6"

  c = collector(2, Int);
  set_relative_orders!(c, [2, 3])
  set_conjugate!(c, 2, 1, [2 => 2])
  gg = pc_group(c)
  @test describe(gg) == "S3"

  c = collector(2, Int);
  set_relative_orders!(c, [2, 3])
  set_commutator!(c, 2, 1, Pair{Int, Int}[])
  gg = pc_group(c)
  @test describe(gg) == "C6"

  c = collector(2, Int);
  set_relative_orders!(c, [2, 3])
  set_commutator!(c, 2, 1, [2 => 1])
  gg = pc_group(c)
  @test describe(gg) == "S3"

  c = collector(4, Int);
  set_relative_orders!(c, [3, 3, 3, 3])
  gg = pc_group(c)
  @test describe(gg) == "C3 x C3 x C3 x C3"

  # infinite polycyclic groups
  c = collector(2);
  set_relative_orders!(c, ZZRingElem[2, 0])
  set_conjugate!(c, 2, 1, [2 => ZZRingElem(-1)])
  gg = pc_group(c)
  @test describe(gg) == "an infinite group"

  # collectors created by hand are independent of collectors stored in groups
  x = GapObj(one(gg))
  cgg = GAP.getbangproperty(x, :collector)
  @test GAP.Globals.IsMutable(cgg)
  @test cgg !== c.X
end
