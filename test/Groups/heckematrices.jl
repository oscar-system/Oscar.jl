n=2
q=3
F = GF(q,1)[1]
G = GL(n,q)
S = SL(n,q)

@testset "Elements of GL(n,q)" begin
   x = matrix(F,2,2,[2,0,0,1])
   y = matrix(F,2,2,[-1,1,-1,0])
   gx = G(x)
   gy = G(y)

   @test F==base_ring(G)
   @test F==base_ring(gx)
   @test x in G
   @test x != gx
   @test parent(x) != G
   @test parent(gx)==G
   @test x[1,1]==F(2)
   @test gx[1,1]==F(2)
   @test order(x)==q-1
   @test !(x in S)
   @test y in S
   @test S(y) in S
   @test_throws AssertionError S(x)

end

@testset "Functions and operations with matrices" begin
   x = matrix(F,2,2,[2,0,0,1])
   y = matrix(F,2,2,[-1,1,-1,0])
   gx = G(x)
   gy = G(y)

   @test gx*gy==G(gx*y)
   @test gx*gy==G(x*gy)
   @test gx*gy==G(x*y)
   @test G(x^-1)==gx^-1
   @test gx^gy==G(x^y)
   @test gx^gy==G(x^gy)
   @test gx^gy==gx^y
   @test G(comm(x,y))==G(x^-1*y^-1*x*y)
   
   cc = conjugacy_class(G,x)
   @test order(cc)==index(G,centralizer(G,x)[1])
   @test !isconjugate(G,x,y)[1]
   @test !isconjugate(G,x,gy)[1]
   @test !isconjugate(G,gx,y)[1]
   z=rand(G)
   @test isconjugate(G,x,x^z)[1]
   w = isconjugate(G,x,x^z)[2]
   @test x^w==x^z

   H = sub(G,[y])[1]
   cc = conjugacy_class(G,H)
   @test order(cc)==index(G,normalizer(G,y)[1])
   @test quo(G,[y])==quo(G,S)
   z = rand(G)
   K=H^z
   @test issubgroup(G,K)[1]
   @test isconjugate(G,H,K)[1]

end

@testset "Homomorphisms with Hecke matrices" begin
   x = matrix(F,2,2,[2,0,0,1])
   y = matrix(F,2,2,[-1,1,-1,0])
   gx = G(x)
   gy = G(y)

   A = automorphism_group(G)
   f = hom(G,G,l -> l^x)
   af = A(f)
   @test f(x)==gx
   @test af(x)==gx
   @test f(y)==G(y^x)
   @test (f^-1)(x)==gx
   @test x^af==gx
   @test x^f==gx
   @test f==hom(G,G,[x,y],[x,y^x])
end
