@testset "MultigradedImplicitization" begin
  R, (x11, x12, x21, x22, y11, y12, y21, y22) = QQ["x11", "x12", "x21", "x22", "y11", "y12", "y21", "y22"]

  S, (s1, s2, t1, t2, m1, m2, n1, n2) = QQ["s1", "s2", "t1", "t2", "m1", "m2", "n1", "n2"]

  phi = hom(R, S, [s1*t1, s1*t2, s2*t1, s2*t2, m1*n1, m1*n2, m2*n1, m2*n2])
  
  components = components_of_kernel(2, phi)
  G = parent(first(keys(components)))
  @test components[G([0, 0, 0, 0, 1, 1, 1, 1])] == [y11 * y22 - y12 * y21]
  @test components[G([1, 1, 1, 1, 0, 0, 0, 0])] == [x11 * x22 - x12 * x21]

  components = oscar_worker_pool(1) do wp
    components_of_kernel(2, phi; wp=wp)
  end

  G = parent(first(keys(components)))
  @test components[G([0, 0, 0, 0, 1, 1, 1, 1])] == [y11 * y22 - y12 * y21]
  @test components[G([1, 1, 1, 1, 0, 0, 0, 0])] == [x11 * x22 - x12 * x21]
end
