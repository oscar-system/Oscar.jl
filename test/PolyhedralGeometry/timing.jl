@testset "Timing" begin
    using Oscar
    using Oscar.Graphs

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
        lp_provide = ["FACETS", "VERTICES", "VERTICES_IN_FACETS", "LATTICE", "BOUNDED"]
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
    

    function benchmark(poly, pref, fun)
        Polymake.prefer("$pref"; application="polytope") do
            copy = deepcopy(poly)
            return @elapsed fun(copy)
        end
    end

    @testset "convex hull" begin
        c4 = create_cutpoly(4)
        @test benchmark(c4, "cdd", facets) <= 1.0
        @test benchmark(c4, "libnormaliz", facets) <= 1.0
        @test benchmark(c4, "ppl", facets) <= 1.0

        ks_6_40 = create_knapsack_ch(6, 40)
        @test benchmark(ks_6_40, "beneath_beyond", facets) <= 1.0
        @test benchmark(ks_6_40, "libnormaliz", facets) <= 1.0
        @test benchmark(ks_6_40, "ppl", facets) <= 1.0

        voronoi_3_1500 = create_voronoi(3, 1500)
        @test benchmark(voronoi_3_1500, "libnormaliz", vertices) <= 1.0
    end

end
