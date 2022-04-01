@testset "AtlasSLProgram" begin
    for s0 in ["", "\n# comment\n \n # comment  \n\n # comment",
               """ # \n  echo  "1 2  4" \necho a echo becho""" ]
        p0 = AtlasSLProgram(s0)
        @test p0.code == s0
        @test p0.ngens == 2
        @test isempty(p0.lines)
        @test p0.outputs == 1:2
        @test endswith(p0.code, "echo") == !isempty(p0.echo)
    end

    s1 = " inp 3 3 1  2 \n oup 2 2  2"
    p1 = AtlasSLProgram(s1)
    @test p1.code == s1
    @test p1.ngens == 3
    @test isempty(p1.lines)
    @test p1.outputs == [3, 3]

    p2 = AtlasSLProgram("inp 2 \n inp 1 4 \n oup 3 4 1 2")
    @test p2.ngens == 3
    @test p2.outputs == [3, 1, 2]

    p3 = AtlasSLProgram("oup 1 2")
    @test p3.ngens == 2
    @test p3.outputs == [2]

    @test_throws ArgumentError AtlasSLProgram("inp 1 2") # implicit "oup 1 2"

    p4 = AtlasSLProgram("inp 1 2 \n oup 2 2 2")
    @test p4.ngens == 1
    @test p4.outputs == [1, 1]

    for bad in ["inp 1 2", "oup 3", "oup 2 \n inp 2",
                "inp 2 \n inp 1", "oup 2 \n oup 1",
                "inp 2 1 2 3", "oup 2 1",
                "inp 3 \n oup 4", "oup 2 2 3",
                "mu 1 2 3 \n inp 3", "oup 2 \n mu 1 2 3"]
        @test_throws ArgumentError AtlasSLProgram(bad)
    end

    for bad in ["cjr 1", "cjr 1 2 3",
                "cj 1 2", "cj 1 2 3 4",
                "com 1 2", "com 1 2 3 4",
                "mu 1 2", "mu 1 2 3 4",
                "pwr 1 2", "pwr 1 2 3 4",
                "iv 1", "iv 1 2 3",
                "cp 1", "cp 1 2 3"]
        @test_throws ArgumentError AtlasSLProgram(bad)
    end
    q1 = AtlasSLProgram("""
           cjr 1 2
           cj 1 2 3
           com 2 3 4
           mu 3 4 5
           pwr 4 5 6
           iv 6 7
           cp 7 8
           oup 3 6 7 8
         """)
    @test q1.ngens == 2
    @test q1.outputs == [6, 7, 8]
    @test q1.lines == [AtlasLine(:cj, 1, 1, 2),
                       AtlasLine(:cj, 3, 1, 2),
                       AtlasLine(:com, 4, 2, 3),
                       AtlasLine(:mu, 5, 3, 4),
                       AtlasLine(:pwr, 6, 4, 5),
                       AtlasLine(:iv, 7, 6),
                       AtlasLine(:cp, 8, 7)
                       ]
    @test_throws ArgumentError AtlasSLProgram("unknown 1 2")
    @test_throws ArgumentError AtlasSLProgram("chor 1 2")

    d1 = AtlasSLDecision("echo something\ninp 2\nchor 1 2\nchor 2 3\nmu 1 2 3\nchor 3 5")
    @test d1.lines ==  [AtlasLine(:chor, 0, 1, 2),
                        AtlasLine(:chor, 0, 2, 3),
                        AtlasLine(:mu, 3, 1, 2),
                        AtlasLine(:chor, 0, 3, 5)]
    @test d1.ngens == 2
    @test d1.echo == ["something"]
    @test_throws ArgumentError AtlasSLDecision("chor 1 2\noup 1")
end

@testset "AtlasSLProgram evaluate" begin
    ab = LazyRec[Gen(:a), Gen(:b)]
    a, b = ab
    res = empty(ab)

    p = AtlasSLProgram("cjr 2 1 \noup 1 2")
    @test SLP.evaluate(p, ab) == [a^-1 * b * a]
    @test res === SLP.evaluate!(res, p, ab) == [a^-1 * b * a]
    g = SLP.compile(GAPSLProgram, p)
    @test g isa GAPSLProgram
    @test SLP.evaluate(g, ab) == [a^-1 * b * a]
    sl = SLP.compile(SLProgram, p)
    @test sl isa SLProgram
    @test SLP.evaluate(sl, ab) == [a^-1 * b * a]
    sl = SLP.compile(p)
    @test sl isa SLProgram
    @test SLP.evaluate(sl, ab) == [a^-1 * b * a]

    p = AtlasSLProgram("cj 1 2 3\noup 2 1 3")
    @test SLP.evaluate(p, ab) == [a, b^-1 * a * b]
    @test res === SLP.evaluate!(res, p, ab) == [a, b^-1 * a * b]
    g = SLP.compile(GAPSLProgram, p)
    @test SLP.evaluate(g, ab) == [a, b^-1 * a * b]

    p = AtlasSLProgram("com 1 2 3 \noup 1 3")
    @test SLP.evaluate(p, ab) == [a^-1 * b^-1 * a * b]
    g = SLP.compile(GAPSLProgram, p)
    @test SLP.evaluate(g, ab) == [a^-1 * b^-1 * a * b]

    p = AtlasSLProgram("iv 1 2")
    @test SLP.evaluate(p, ab) == [a, a^-1]
    g = SLP.compile(GAPSLProgram, p)
    @test SLP.evaluate(g, ab) == [a, a^-1]

    p = AtlasSLProgram("mu 1 2 1")
    @test SLP.evaluate(p, ab) == [a * b, b]
    g = SLP.compile(GAPSLProgram, p)
    @test SLP.evaluate(g, ab) == [a * b, b]

    p = AtlasSLProgram("inp 1 \n pwr 4 1 2")
    @test SLP.evaluate(p, ab) == [a, a^4]
    g = SLP.compile(GAPSLProgram, p)
    @test SLP.evaluate(g, ab) == [a, a^4]

    p = AtlasSLProgram("cp 1 3 \n cp 2 4 \n oup 4")
    @test SLP.evaluate(p, ab) == [a, b, a, b]
    g = SLP.compile(GAPSLProgram, p)
    @test SLP.evaluate(g, ab) == [a, b, a, b]

    @test_throws ArgumentError SLP.compile(GAPSLDecision, p)

    # bug with overwriting results
    p = AtlasSLProgram("oup 2 2 1")
    @test SLP.evaluate(p, ab) == [b, a]
    g = SLP.compile(GAPSLProgram, p)
    @test SLP.evaluate(g, ab) == [b, a]
end

@testset "AtlasSLDecision evaluate /compile" begin
    p = perm"(1, 2)(3, 4)"
    q = Perm([1, 4, 2, 3])
    pq = [p, q]

    # same test as in tests for GAPSLDecision
    for (i, j, k, r) in [(2, 3, 5, false),
                         (2, 3, 3, true),
                         (2, 2, 3, false),
                         (1, 3, 3, false)]
        d = AtlasSLDecision(
            """
            mu 1 2 3
            chor 1 $i
            chor 2 $j
            chor 3 $k
            """)

        @test SLP.evaluate(d, pq) == r
        @test_throws ArgumentError SLP.compile(GAPSLProgram, d)
    end

    for (p, q, r) in eachcol(rand(parent(p), 3, 50))
        pqr = [p, q, r]
        for (i, j, k) in eachcol(rand(1:6, 3, 10))
            ad = AtlasSLDecision(
                 """
                 inp 3
                 chor 1 $i
                 chor 2 $j
                 chor 3 $k
                 """)
            gd = GAPSLDecision([
                ["Order", 1, i],
                ["Order", 2, j],
                ["Order", 3, k]],
                               3)
            @test SLP.evaluate(ad, pqr) == SLP.evaluate(gd, pqr)
            @test SLP.compile(GAPSLDecision, ad).lines == gd.lines
        end
    end
end
