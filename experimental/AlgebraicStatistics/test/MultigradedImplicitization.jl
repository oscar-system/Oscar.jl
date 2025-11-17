using Distributed

oscar_worker_pool(1) do wp
  @testset "MultigradedImplicitization" begin
    R, (x11, x12, x21, x22, y11, y12, y21, y22) = QQ["x11", "x12", "x21", "x22", "y11", "y12", "y21", "y22"]
    S, (s1, s2, t1, t2, m1, m2, n1, n2) = QQ["s1", "s2", "t1", "t2", "m1", "m2", "n1", "n2"]

    phi = hom(R, S, [s1*t1, s1*t2, s2*t1, s2*t2, m1*n1, m1*n2, m2*n1, m2*n2])
    
    comps = components_of_kernel(2, phi)
    G = parent(first(keys(comps)))
    for (_, v) in comps
      @test all(is_zero, phi.(v))
    end

    @test isone(length(comps[G([0, 0, 0, 0, 1, 1, 1, 1])]))

    for (d, v) in comps
      @test all(parent(p) === R for p in v)
    end
  end

  @testset "K3P" begin
    T = graph_from_edges(Directed, [[6, 5], [5, 1], [5, 2], [6, 3], [6, 4]])
    M = kimura3_model(T)
    phi = parametrization(M)
    _, q = model_ring(M)
    comps = components_of_kernel(2, phi; wp=wp)
    for (_, v) in comps
      @test all(is_zero, phi.(v))
    end
  end
end
