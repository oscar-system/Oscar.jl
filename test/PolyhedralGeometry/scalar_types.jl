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
    E = Hecke.embedded_field(L, real_embeddings(L)[2])
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
    ff = face_fan(sd)
    @test vector_matrix(rays(ff)) == V

end
