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
  set_relative_orders!(c, fmpz[2, 0])
  set_conjugate!(c, 2, 1, [2 => fmpz(-1)])
  gg = pc_group(c)
  @test describe(gg) == "an infinite group"

  # collectors created by hand are independent of collectors stored in groups
  x = one(gg).X
  cgg = GAP.getbangproperty(x, :collector)
  @test GAP.Globals.IsMutable(cgg)
  @test ! GAP.Globals.IsIdenticalObj(cgg, c.X)
end
