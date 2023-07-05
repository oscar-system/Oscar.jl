@testset "Scalar types" begin

    # small helper for quantitative testing
    # returns the scalar product of the difference of the vertices of an edge x (which is encoded as a Polyhedron),
    # i.e. the square of the edge length
    function _edge_length_for_test(x::Polyhedron)
        v = vertices(x)
        s = v[1] - v[2]
        return sum([s[i]^2 for i in 1:3])
    end
    
    # this example relies on the same source as polymake's johnson_solid(84):
    # https://de.wikipedia.org/wiki/Trigondodekaeder
    # as of now, johnson solid J84 is realized as a Polyhedron over floats within polymake
    Qx, x = QQ["x"]
    K, r = number_field(x^3 - 3x^2 - 4x + 8, "r")
    Ky, y = K["y"]
    L, = number_field(y^2 - (2-r^2)//2, "q")
    Lz, z = L["z"]
    E, q = Hecke.embedded_field(L, real_embeddings(L)[2])
    pq1 = E(roots(z^2 - (3 + 2r - r^2))[2])
    pq2 = E(roots(z^2 - (3 - r^2))[1])
    p = (pq1 + pq2)//2
    r = E(r)
    
    v = [0 1 p; 0 -1 p; r 0 q; -r 0 q; 0 r -q; 0 -r -q; 1 0 -p; -1 0 -p]
    V = matrix(E, v)
    sd = convex_hull(E, v; non_redundant = true)
    
    @test sd isa Polyhedron{Hecke.EmbeddedNumFieldElem{Hecke.NfRelElem{nf_elem}}}
    @test point_matrix(vertices(sd)) == V
    # volume formula from source
    @test volume(sd) == 8r//12 * (6pq2 + 2*r*q)
    # the snub disphenoid consists of 12 triangles
    @test nvertices.(faces(sd, 2)) == repeat([3], 12)
    # and here the 18 edges are of length 2
    @test _edge_length_for_test.(faces(sd, 1)) == repeat([4], 18)
    # scaling the Polyhedron by 3 yields edge lengths of 6
    @test _edge_length_for_test.(faces(3*sd, 1)) == repeat([36], 18)
    let pc = polyhedral_complex(E, IncidenceMatrix(facets(sd)), vertices(sd); non_redundant = true)
        @test maximal_polyhedra(pc) == faces(sd, 2)
    end
    let c = convex_hull(E, permutedims([0]), permutedims([r]))
        ms = product(sd, c)
        @test recession_cone(ms) == positive_hull(E, [0 0 0 1])
    end
    nf = normal_fan(sd)
    nfc =  polyhedral_fan(E, rays(nf), maximal_cones(IncidenceMatrix,nf))
    @test is_regular(nfc)

    @testset "Scalar detection" begin
        let j = johnson_solid(12)
            @test j isa Polyhedron{QQFieldElem}
            s = vertices(j)[1][1]
            @test s isa QQFieldElem
        end
        let j = johnson_solid(64)
            @test j isa Polyhedron{Float64}
            s = vertices(j)[1][1]
            @test s isa Float64
        end
        let a = archimedean_solid("truncated_cube")
            @test a isa Polyhedron{Hecke.EmbeddedNumFieldElem{nf_elem}}
            s = vertices(a)[1][1]
            @test s isa Hecke.EmbeddedNumFieldElem{nf_elem}
            f = Oscar.get_parent_field(a)
            F = number_field(f)
            isq = Hecke.is_quadratic_type(F)
            @test isq[1]
            @test isq[2] == 2
        end
    end

end
