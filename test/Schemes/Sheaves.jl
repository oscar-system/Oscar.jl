@testset "constant sheaf of integers on Spec" begin
  R, x = QQ["x", "y", "z"]
  X = Spec(R)

  is_open_func(U::AbsSpec, V::AbsSpec) = is_open_embedding(U, V)

  function production_func(U::AbsSpec)
    return ZZ
  end

  function restriction_func(V::AbsSpec, U::AbsSpec)
    return identity_map(ZZ)
  end

  const_sheaf_ZZ= SheafOnSpec(X, production_func, restriction_func,
                              OpenType=AbsSpec, OutputType=typeof(ZZ), 
                              RestrictionType=typeof(identity_map(ZZ)),
                              is_open_func=is_open_func
                             )

  U = hypersurface_complement(X, [x[1]])
  G = const_sheaf_ZZ(U)
  @test G === const_sheaf_ZZ(U)
  V = hypersurface_complement(U, [x[2]])
  @test G === const_sheaf_ZZ(V)
  rho = const_sheaf_ZZ(U, V)
  @test rho==identity_map(ZZ)
  @test domain(rho) === const_sheaf_ZZ(U)
  @test codomain(rho) === const_sheaf_ZZ(V)
end
