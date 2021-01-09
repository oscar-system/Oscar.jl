@testset "mpoly-graded" begin

    Qx, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    t = gen(Hecke.Globals.Qx)
    k1 , l= number_field(t^5+t^2+2)
    NFx = PolynomialRing(k1, ["x", "y", "z"])[1]
    k2 = Nemo.GF(23)
    GFx = PolynomialRing(k2, ["x", "y", "z"])[1]
    RNmodx=PolynomialRing(Nemo.ResidueRing(ZZ,17), :x => 1:2)[1]
    Rings= [Qx, NFx, GFx, RNmodx]

    A = abelian_group([4 3 0 1;   0 0 0 3])
    GrpElems = [A([convert(fmpz,x) for x= [1,0,2,1]]), A([convert(fmpz,x) for x= [1,1,0,0]]), A([convert(fmpz,x) for x= [2,2,1,0]])]

    T = grade(Qx)
    @test sprint(show, "text/plain", T(x)//T(1)) isa String

    Rings_dec=[]
    v= [1,1,2]
    for R in Rings
        temp = []
        push!(Rings_dec, [decorate(R), grade(R, [v[i] for i=1:ngens(R)]), filtrate(R, [v[i] for i=1:ngens(R)]), filtrate(R, [GrpElems[i] for i =1:ngens(R)],(x,y) -> x[1]+x[2]+x[3]+x[4] < y[1]+y[2]+y[3]+y[4]), grade(R, [GrpElems[i] for i =1:ngens(R)])])
    end

    function rmPols(i::Int64, j::Int64)
        Pols = Array{elem_type(Rings_dec[i][j]),1}()
        if i == 2
            coeff = fill(zero(k1),3,4)
            for k = 1:4
                for j = 1:3
                    coeff[j,k] = dot(rand(-10:10,5),[one(k1),l,l^2,l^3,l^4])
                end
            end
        else   
            init = rand(0:22,4,4)
            coeff = [(base_ring(Rings_dec[i][j]))(x) for x = init]
        end
        for t = 1:4
            f = MPolyBuildCtx(Rings_dec[i][j])
            g = MPolyBuildCtx(Rings_dec[i][j].R)
            for z = 1:3
                e = rand(0:6, ngens(Rings_dec[i][j].R))
                push_term!(f, coeff[z,t], e)
                push_term!(g, coeff[z,t], e)
            end
            f = finish(f)
            @test Rings_dec[i][j](finish(g)) == Rings_dec[i][j](f)
            push!(Pols, f)
        end
        return(Pols)
    end

    function homPols(Polynomials::Array{<:MPolyElem,1})
        R = parent(Polynomials[1])
        D = Array{Any}(undef,4)
        Monomials = [zero(R)]
        for i = 1:4
            D[i] = homogenous_components(Polynomials[i])
            Monomials = vcat(Monomials, collect(Oscar.monomials(Polynomials[i])))
        end
        homPolys = []
        for deg in unique([degree(mon) for mon = Monomials])
            g = R(0)
            for i = 1:4
                if haskey(D[i], deg)
                    temp = get(D[i], deg, 'x')
                    @test homogenous_component(Polynomials[i], deg) == temp
                    g += temp
                end
            end
        @test ishomogenous(g)
        push!(homPolys, g)
        end    
        return homPolys
    end

    for i = 1:4
        @test !Oscar.isgraded(Rings_dec[i][1])
        @test Oscar.isgraded(Rings_dec[i][2])
        @test !Oscar.isgraded(Rings_dec[i][3])
        @test !Oscar.isgraded(Rings_dec[i][4])
        @test Oscar.isgraded(Rings_dec[i][5])
        @test Oscar.isfiltrated(Rings_dec[i][1])
        @test !Oscar.isfiltrated(Rings_dec[i][2])
        @test Oscar.isfiltrated(Rings_dec[i][3])
        @test Oscar.isfiltrated(Rings_dec[i][4])
        @test !Oscar.isfiltrated(Rings_dec[i][5])
    end
    d_Elems = Array{Any, 1}()
    for i= 1:3
        push!(d_Elems, [4*Rings_dec[i][1].D[1], 5*Rings_dec[i][2].D[1], 3*Rings_dec[i][3].D[1], Rings_dec[i][4].D([2,2,1,0]), Rings_dec[i][5].D([2,2,1,0])])
    end
    Dimensions = [15, 12, 6, 1, 1]
    Polys = Array{Any}(undef,4,5)
    for i = 1:4, j=1:5
        base_ring(Rings_dec[i][j])
        @test ngens(Rings_dec[i][j]) == length(gens(Rings_dec[i][j]))
        @test gen(Rings_dec[i][j], 1) == Base.getindex(Rings_dec[i][j], 1)

        Polys[i,j] = rmPols(i,j)
        @test one(Rings_dec[i][j]) == Rings_dec[i][j](1)
        @test zero(Rings_dec[i][j])== Rings_dec[i][j](0)
        @test iszero(zero(Rings_dec[i][j]))
        @test !iszero(one(Rings_dec[i][j]))
        @test isone(one(Rings_dec[i][j]))
        @test !isone(zero(Rings_dec[i][j]))
        @test divexact(one(Rings_dec[i][j]), one(Rings_dec[i][j])) == one(Rings_dec[i][j])
        @test (Polys[i,j][1] + Polys[i,j][2])^2 == Polys[i,j][1]^2 + 2*Polys[i,j][1]*Polys[i,j][2] + Polys[i,j][2]^2
        @test (Polys[i,j][3] - Polys[i,j][4])^2 == Polys[i,j][3]^2 + 2*(-Polys[i,j][3])*Polys[i,j][4] + Polys[i,j][4]^2
        @test Polys[i,j][2] * (Polys[i,j][3] + Polys[i,j][4]) == Oscar.addeq!(Oscar.mul!(Polys[i,j][1], Polys[i,j][2], Polys[i,j][3]), Oscar.mul!(Polys[i,j][1], Polys[i,j][2], Polys[i,j][4]))
        @test parent(Polys[i,j][1]) == Rings_dec[i][j]    
        for k= 1:Oscar.length(Polys[i,j][4])
            @test Oscar.coeff(Polys[i,j][4],k) * Oscar.monomial(Polys[i,j][4], k) == finish(push_term!(MPolyBuildCtx(Rings_dec[i][j]), collect(Oscar.MPolyCoeffs(Polys[i,j][4]))[k], collect(Oscar.MPolyExponentVectors(Polys[i,j][4]))[k]))
        end
        homogenous_Polys = homPols(Polys[i,j])
        I = ideal(homogenous_Polys)
        R_quo = Oscar.MPolyQuo(Rings_dec[i][j], I)
        @test base_ring(R_quo) == Rings_dec[i][j]
        @test modulus(R_quo) == I
        f = R_quo(Polys[i,j][2])
        D = homogenous_components(f)
        for deg in [degree(R_quo(mon)) for mon  = collect(Oscar.monomials(f.f))]
            h = get(D, deg, 'x')
            @test ishomogenous(R_quo(h))
            @test h == homogenous_component(f, deg)        
        end
        if j == 1 || j== 3 || j==4
            @test Oscar.isfiltrated(R_quo)
        else
            @test !Oscar.isfiltrated(R_quo)
        end
        if j == 2 || j==5
            @test Oscar.isgraded(R_quo)
        else
            @test !Oscar.isgraded(R_quo)
        end
        @test decoration(R_quo) == decoration(Rings_dec[i][j])
        if i!= 4
            d_GrpElems = d_Elems[i]
            H = homogenous_component(Rings_dec[i][j], d_GrpElems[j])
            @test Oscar.hasrelshp(H[1], Rings_dec[i][j]) !== nothing
            for g in gens(H[1])
                @test degree(H[2](g)) == d_GrpElems[j]
                @test (H[2].g)(Rings_dec[i][j](g)) == g
            end
            @test dim(H[1]) == Dimensions[j]
            #H_quo = homogenous_component(R_quo, d_GrpElems[j])
            #Oscar.hasrelshp(H_quo[1], R_quo) !== nothing
            #for g in gens(H_quo[1])
            #    degree(H_quo[2](g)) == d_GrpElems[j]
            #    (H_quo[2].g)(R_quo(g)) == g
            #end
        end
    end      
        
    l_dec=Rings_dec[2][1](l)    
end
