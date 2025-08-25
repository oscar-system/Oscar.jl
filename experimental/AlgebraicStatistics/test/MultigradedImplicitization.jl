using Distributed

oscar_worker_pool(1) do wp
  @testset "MultigradedImplicitization" begin
    R, (x11, x12, x21, x22, y11, y12, y21, y22) = QQ["x11", "x12", "x21", "x22", "y11", "y12", "y21", "y22"]

    S, (s1, s2, t1, t2, m1, m2, n1, n2) = QQ["s1", "s2", "t1", "t2", "m1", "m2", "n1", "n2"]

    phi = hom(R, S, [s1*t1, s1*t2, s2*t1, s2*t2, m1*n1, m1*n2, m2*n1, m2*n2])
    
    components = components_of_kernel(2, phi)
    G = parent(first(keys(components)))
    @test components[G([0, 0, 0, 0, 1, 1, 1, 1])] == [y11 * y22 - y12 * y21]
    @test components[G([1, 1, 1, 1, 0, 0, 0, 0])] == [x11 * x22 - x12 * x21]


    @time components_of_kernel(2, phi; wp=wp)

    G = parent(first(keys(components)))
    @test components[G([0, 0, 0, 0, 1, 1, 1, 1])] == [y11 * y22 - y12 * y21]
    @test components[G([1, 1, 1, 1, 0, 0, 0, 0])] == [x11 * x22 - x12 * x21]
  end


  @testset "K3P" begin
    G = abelian_group(2, 2)
    elG = collect(G)
    T = graph_from_edges(Directed, [[6, 5], [5, 1], [5, 2], [6, 3], [6, 4]])
    M = kimura3_model(T)

    inds = [g for g in collect(Iterators.product([elG for _ in 1:4]...)) if sum(g) == elG[1]]
    S, p = polynomial_ring(QQ, "q" => inds)
    q = Dict()
    for i in 1:length(p)
      q[inds[i]] = p[i]
    end

    param_inds = [(i, g) for i in 1:8 for g in elG]
    R, b = polynomial_ring(QQ, "a" => param_inds)
    a = Dict()
    for i in 1:length(param_inds)
      ind = param_inds[i]
      a[ind[1], ind[2]] = b[i]
    end

    images = []

    for g in inds

      im = a[1, g[1]]*a[2, g[2]]*a[3, g[3]]*a[4, g[4]]*(a[5, g[1]]*a[6, g[1]+g[2]]*a[7, g[4]] + a[6, g[2]]*a[7, g[2]+g[3]]*a[8, g[1]])
      push!(images, im)
    end

    phi = hom(S, R, images);
    components_of_kernel(2, phi; wp=wp)
  end
end
