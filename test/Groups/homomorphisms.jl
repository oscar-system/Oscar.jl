
n = 6
@testset "Homomorphism in Sym($n)" begin
   G = symmetric_group(n)
   x = rand(G)
   while x != one(G)
      x = rand(G)
   end

   g = hom(G, G, gens(G), [y^x for y in gens(G)])
   @test typeof(g) == Oscar.GAPGroupHomomorphism{PermGroup, PermGroup}
   @test domain(g) == G
   @test codomain(g) == G
   @test image(g)[1] == G
   @test g(one(G)) == one(G)
   @test g(x)==x
   @test Set(elements(centralizer(G,x)[1])) == Set([y for y in G if g(y)==y])
   z = rand(G)
   @test g(z) == z^x
   @test isinjective(g)
   @test issurjective(g)
   @test isinvertible(g)
   @test isbijective(g)
   og = order(g)
   @test og isa Integer
   @test og == order(x)
   @test g^og == id_hom(G)
   @test g^(og+1) == g
   @test g^(1-og) == g

end

@testset "Operations on homomorphism in Sym($n)" begin
   G = symmetric_group(n)
   x = cperm(G,[1,2,3])
   Hx,fx = sub(G,[x])

   @test isinjective(fx)
   @test !issurjective(fx)
   @test fx(Hx[1])==x
   @test image(fx)[1] == Hx
   @test image(fx)[2] == fx

   y = cperm(G,[1,3,4])
   z = cperm(G,[1,2,5])
   Hy,fy = sub(G,[y])
   Hz,fz = sub(G,[z])
   f = hom(Hx,Hy,gens(Hx),gens(Hy))
   g = hom(Hy,Hz,gens(Hy),gens(Hz))
   @test isbijective(f)
   @test domain(f)==Hx
   @test codomain(f)==Hy
   @test image(f)[1]==Hy
   @test image(f)[2]==id_hom(Hy)
   @test f(x)==y
   @test inv(f)==f^-1
   @test domain(f^-1)==codomain(f)
   @test codomain(f^-1)==domain(f)
   @test image(f^-1)[1]==Hx
   @test image(f^-1)[2]==id_hom(Hx)
   @test (f^-1)(y)==x
   @test g(f(x))==z
   @test (f*g)(x)==z
   @test (f*g)^-1 == g^-1*f^-1
   @test_throws AssertionError g*f
   ty = trivial_morphism(Hy,Hy)
   @test f*ty==trivial_morphism(Hx,Hy)
   @test ty*g==trivial_morphism(Hy,Hz)
   @test f*trivial_morphism(Hy,Hz) == trivial_morphism(Hx,Hz)
   @test trivial_morphism(Hz,Hy)*f^-1 == trivial_morphism(Hz,Hx)
   @test order(kernel(f)[1])==1
   @test kernel(ty)[1] == Hy
   @test kernel(ty)[2] == id_hom(Hy)
   @test kernel(f*ty)[1]==Hx
   
end

@testset "Isomorphic groups" begin
   @testset "Dihedral_as_permutation" for n in 4:10
      G = symmetric_group(n)
      D = dihedral_group(PermGroup,2*n)
      H = sub(G, [G(D[1]),G(D[2])])[1]
   
      @test order(H)==2*n
      @test H == D
   end

   @testset "Abelian_as_permutation" for n in 15:20
      G = symmetric_group(n)
      for j in 2:n-2
         D = abelian_group(PermGroup,[j,n-j])
         H = sub(G, [G(D[1]),G(D[2])])[1]

         @test H == D
         @test order(D) == j*(n-j)
         @test (gcd(j,n-j)==1 && iscyclic(D)) || (gcd(j,n-j)>1 && !iscyclic(D))
      end
   end

end

