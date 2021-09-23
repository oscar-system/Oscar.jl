@testset "AlgebraicGroups" begin
  G = AlgGroup_GL(2,QQ)
  rep = LinearRepresentation(G)
  rep3 = representation_on_sym_power( rep, 3 )
end
