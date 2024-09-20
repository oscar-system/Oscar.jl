using Distributed

process_ids = addprocs(1)

@everywhere using Oscar

@testset "Interprocess Serialization" begin
  channels = Oscar.params_channels(Union{Ring, MatSpace})

  Qx, x = QQ["x"]
  F, a = number_field(x^2 + x + 1)
  MR = matrix_space(F, 2, 2)

  Oscar.put_params(channels, Qx)
  Oscar.put_params(channels, F)
  Oscar.put_params(channels, MR)
  
  c = [MR([a^i F(1); a a + 1]) for i in 1:5]
  dets = pmap(det, c)
  total = reduce(*, dets)

  @test total == F(4)
end

map(rmprocs, process_ids)

process_ids = addprocs(3)
@everywhere using Oscar

@testset "parallel smoothness test for schemes: proof-of-concept" begin
  channels = Oscar.params_channels(Union{AffineScheme, Ring})
  
  IP2 = projective_space(QQ, [:x, :y, :z])
  S = homogeneous_coordinate_ring(IP2)
  (x, y, z) = gens(S)
  
  I = ideal(S, y^2*z + x^3 + x^2*z)
  
  X, _ = sub(IP2, I)
  X_cov = covered_scheme(X)
  U = affine_charts(X_cov)
  for a in U
    Oscar.put_params(channels, ambient_coordinate_ring(a))
    Oscar.put_params(channels, OO(a))
  end
  results = pmap(is_smooth, U)
  @test results == [1, 1, 0]
end  

map(rmprocs, process_ids)


process_ids = addprocs(3)
@everywhere using Oscar

@testset "parallel smoothness test for schemes" begin
  X = rational_d9_pi6()
  #set_assertion_scope(:Parallelization, 5)
  #set_verbosity_scope(:Parallelization, 5)
  #IP4 = projective_space(QQ, [:x, :y, :z, :v, :w])
  #S = homogeneous_coordinate_ring(IP4)
  #(x, y, z, v, w) = gens(S)
  
  #I = ideal(S, sum(u^6 for u in gens(S)))
  #X, _ = sub(IP4, I)
  X_cov = covered_scheme(X)
  Oscar.is_smooth_parallel2(X_cov) # works!
end  

map(rmprocs, process_ids)

