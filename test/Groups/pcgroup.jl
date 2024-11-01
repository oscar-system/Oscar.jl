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

  G = small_group(24, 12)
  H = sylow_subgroup(G, 2)[1]
  F, mp = full_group(H)
  @test F == G
  @test F !== small_group(24, 12)
  @test domain(mp) == H && codomain(mp) == G
  @test full_group(G)[1] == G
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

@testset "generate letters from polycyclic group element" begin

  # finite polycyclic groups
  c = collector(2, Int);
  set_relative_order!(c, 1, 2)
  set_relative_order!(c, 2, 3)
  set_power!(c, 1, [2 => 1])
  gg = pc_group(c)
  @test letters(gg[1]^5*gg[2]^-4) == [1, 2]
  @test letters(gg[1]^5*gg[2]^4) == [1] # all positive exp
  @test letters(gg[1]^-5*gg[2]^-7) == [1, 2, 2] # all negative exp
  @test letters(gg[1]^2*gg[2]^3) == [2] # both identity elements

  # finite polycyclic subgroup
  gg = pc_group(symmetric_group(4))
  G = derived_subgroup(gg)[1]
  @test letters(G[1]^2) == [2, 2]
  @test letters(G[1]^2*G[2]^3*G[3]^3) == [2, 2, 3, 4]
  @test letters(G[1]^-2*G[2]^-3*G[3]^-3) == [2, 3, 4]
end

@testset "create polycyclic group element from syllables" begin

  # finite polycyclic groups
  c = collector(2, Int);
  set_relative_order!(c, 1, 2)
  set_relative_order!(c, 2, 3)
  set_power!(c, 1, [2 => 1])
  gg = pc_group(c)

  element = gg[1]^5*gg[2]^-4
  sylls = syllables(element)
  @test sylls == [1 => ZZ(1), 2 => ZZ(1)] # check general usage
  @test gg(sylls) == element # this will pass the check

  sylls = [1 => ZZ(1), 2 => ZZ(2), 1 => ZZ(3)]
  @test_throws ArgumentError gg(sylls) # repeating generators

  sylls = [2 => ZZ(1), 1 => ZZ(2)]
  @test_throws ArgumentError gg(sylls) # not in ascending order

  sylls = [2 => ZZ(1), 1 => ZZ(2), 1 => ZZ(3)] # both conditions
  @test_throws ArgumentError gg(sylls)
end
