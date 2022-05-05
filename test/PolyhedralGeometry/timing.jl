@testset "Timing" begin
    using Oscar
    using Oscar.Graphs


    lp_provide = ["FACETS", "VERTICES", "VERTICES_IN_FACETS", "LATTICE", "BOUNDED"]

    function create_matching(n)
        g = complete_graph(n)
        p = Polymake.polytope.fractional_matching_polytope(g.pm_graph);
        Polymake.prefer("ppl") do
            for prop in lp_provide
                Polymake.give(p, prop)
            end
        end
        Polymake.setname!(p, "matching($n)")
        return Polyhedron(p)
    end


    function create_randbox(d, n)
        p = Polymake.polytope.rand_box(d, n, 5, seed=42)
        Polymake.prefer("ppl") do
            for prop in lp_provide
                Polymake.give(p, prop)
            end
        end
        Polymake.setname!(p, "rand_box($d,$n,5)")
        return Polyhedron(p)
    end


    function create_knapsack_lp(d, b)
        (first, second) = (2, 3)
        ineq = [b, -first, -second];
        for i in 2:d-1
            third = first + second;
            first = second;
            second = third;
            push!(ineq, -third)
        end
        p = Polymake.polytope.fractional_knapsack(ineq);
        Polymake.prefer("ppl") do
            for prop in lp_provide
                Polymake.give(p, prop)
            end
        end
        Polymake.setname!(p, "knapsack($d,$b)")
        return Polyhedron(p)
    end


    function create_knapsack_ch(d, b)
        k = create_knapsack_lp(d, b)
        points = lattice_points(k) |> Oscar.matrix_for_polymake |> Oscar.homogenize
        #p = convex_hull(points)
        #for prop in ["BOUNDED", "POINTED", "FEASIBLE"]
        #   Polymake.give(Oscar.pm_polytope(p), prop)
        #end
        p = Polymake.polytope.Polytope(POINTS=points, BOUNDED=true, POINTED=true)
        Polymake.setname!(p, "knapsack_pts($d,$b)")
        return Polyhedron(p)
    end


    function create_cutpoly(n)
        g = Graph{Undirected}(n+6)
        add_edge!(g,1,2)
        add_edge!(g,2,3)
        add_edge!(g,2,5)
        add_edge!(g,3,4)
        add_edge!(g,4,5)
        add_edge!(g,5,6)
        for i in 1:n
            add_edge!(g,5+i,6+i)
        end
        c = Polymake.polytope.fractional_cut_polytope(g.pm_graph)
        Polymake.give(c,"N_RAYS | N_INPUT_RAYS")
        Polymake.give(c,"POINTED")
        Polymake.setname!(c,"non-sym-cutpoly($n)")
        return Polyhedron(c);
    end


    function create_voronoi(d, n)
        c = cube(d-1);
        r = Polymake.polytope.rand_inner_points(Oscar.pm_polytope(c), n);
        p = Polymake.polytope.VoronoiPolyhedron(SITES=r.POINTS);
        Polymake.give(p,"FACETS")
        Polymake.setname!(p, "rand-voronoi($d,$n)");
        return Polyhedron(p);
    end
    

    function benchmark(poly, pref, fun, bound)
        # Since timings can be unstable, especially with the jit compiler, we
        # repeat every experiment and return the minimal running time.
        repeat = 5
        Polymake.prefer("$pref"; application="polytope") do
            result = 10
            for i in 1:repeat
                copy = deepcopy(poly)
                result = min(@elapsed fun(copy), result)
                if result <= bound
                    return true
                end
            end
            return false
        end
    end

    @testset "convex hull" begin
        c4 = create_cutpoly(4)
        @test benchmark(c4, "cdd", facets, 2.0)
        @test benchmark(c4, "libnormaliz", facets, 1.0)
        @test benchmark(c4, "ppl", facets, 1.0)

        ks_6_40 = create_knapsack_ch(6, 40)
        @test benchmark(ks_6_40, "beneath_beyond", facets, 3.0)
        @test benchmark(ks_6_40, "libnormaliz", facets, 1.0)
        @test benchmark(ks_6_40, "ppl", facets, 1.0)

        voronoi_3_1500 = create_voronoi(3, 1500)
        @test benchmark(voronoi_3_1500, "libnormaliz", vertices, 3.0)
    end


    @testset "lattice points" begin
        ks_10_60 = create_knapsack_lp(10, 60)
        @test benchmark(ks_10_60, "projection", lattice_points, 2.0)
        @test benchmark(ks_10_60, "libnormaliz", lattice_points, 1.0)

        m6 = create_matching(6)
        @test benchmark(m6, "_4ti2", lattice_points, 1.0)
        @test benchmark(m6, "projection", lattice_points, 1.0)
        @test benchmark(m6, "libnormaliz", lattice_points, 1.0)

        rbx = create_randbox(4, 70)
        @test benchmark(rbx, "bbox", lattice_points, 2.0)
        @test benchmark(rbx, "projection", lattice_points, 1.0)
        @test benchmark(rbx, "libnormaliz", lattice_points, 1.0)
    end
end
