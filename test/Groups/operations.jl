
@testset "Operations on PermGroup" begin
   @testset "Elements of Sym($i)" for i in 8:16
      if i>1
      G=symmetric_group(i)
      x,y = @inferred rand(G, 2)
      #z=cperm(1:i)
      z=G(vcat([j for j in 2:i],[1]))
      if i>3
         w=cperm([1,2],[j for j=3:i])
      else
         w=cperm([1,2])
      end
      ox=order(x)
      oy=order(y)
      oz=order(z)

      @test x isa PermGroupElem
      @test ox isa fmpz
      @test inv(x) in G
      @test mul(x,y) == x*y
      @test mul!(x,x,y) == x*y
      @test inv(x)==x^-1
      @test inv!(x,x) == x^-1
      @test sign(z)==(-1)^(i-1)
      @test sign(x*y)==sign(x)*sign(y)
      @test parent(x)==G
      @test parent(z)==G
      @test x/y == x*y^-1
      @test (x*y).X == x.X*y.X
      @test x^2 == x*x
      @test x*one(G)==x
      @test oz==i
      @test x*x^-1 == one(G)
      @test isone(x*x^-1)
      @test x^(ox-1)==x^-1
      @test x^w == w^-1*x*w
      @test z^y == y^-1*z*y
      @test conj(x,y) == y^-1*x*y
      @test comm(x,y) == x^-1*conj(x,y)
      @test comm(x,y) == div_left(x^y,x)
      @test div_right(x,y) == x/y
      @test z(1)==2
      @test (z^-1)(i)==i-1
      @test x(y(1))==(y*x)(1)
      @test x>y || x==y || x<y
      @test isequal(x,y) || isless(x,y) || isless(y,x)
      @test (isless(x,y) && x<y) || (isequal(x,y) && x==y) || (isless(y,x) && x>y)
      @test y==G(Vector(y))
      end
   end
end

@testset "Matrix manipulation" begin
   F = GF(5, 1)
   
   I = identity_matrix(F,6)
   x = matrix(F,2,2,[2,3,4,0])
   I1 = deepcopy(I)
   I1[3:4,3:4] = x
   @test I==identity_matrix(F,6)
   I[3:4,3:4] = x
   @test I==I1
   @test I[3:4,3:4]==x
   V = VectorSpace(F,6)
   @test matrix([V[i] for i in 1:6])==identity_matrix(F,6)
   L = [1,4,6,2,3,5]
   P = permutation_matrix(F,L)
   @testset for i in 1:6
      @test V[i]*P == V[L[i]]
   end
   p = symmetric_group(6)(L)
   @test permutation_matrix(F,p)==permutation_matrix(F,L)
   @test upper_triangular_matrix(F.([2,3,1,1,0,1]))==matrix(F,3,3,[2,3,1,0,1,0,0,0,1])
   @test lower_triangular_matrix(F.([2,3,1,1,0,1]))==matrix(F,3,3,[2,0,0,3,1,0,1,0,1])
   @test_throws ArgumentError lower_triangular_matrix(F.([2,3,1,1]))

   R,t=PolynomialRing(F,"t")
   f = t^4+2*t^3+4*t+1
   @test f(identity_matrix(F,6))==f(1)*identity_matrix(F,6)
   @test_throws ArgumentError conjugate_transpose(x)
   @test is_symmetric(P+transpose(P))
   @test is_skewsymmetric_matrix(P-transpose(P))

   F,z = FiniteField(2,2)
   x=matrix(F,4,4,[1,z,0,0,0,1,z^2,z,z,0,0,1,0,0,z+1,0])
   y=x+transpose(x)
   @test is_symmetric(y)
   @test is_hermitian_matrix(x+conjugate_transpose(x))
   @test is_skewsymmetric_matrix(y)
   y[1,1]=1
   @test !is_skewsymmetric_matrix(y)
   @test conjugate_transpose(x)==transpose(matrix(F,4,4,[1,z+1,0,0,0,1,z,z+1,z+1,0,0,1,0,0,z,0]))

end

@testset "Operations with vector spaces" begin
   F= GF(7, 1)
   V=VectorSpace(F,5)

   @test V([1,0,2,0,6])==V[1]+2*V[3]-V[5]
   U = sub(V,[V[1],V[3]])[1]
   W = complement(V,U)[1]
   @test dim(intersect(U,W)[1])==0
   @test dim(W)==3
   W0 = sub(V,[])[1]
   @test complement(V,W0)[1]==V
   @test complement(V,sub(V,gens(V))[1])[1]==W0

   v1=V([1,2,3,4,5])
   v2=V([1,6,0,5,2])
   G = GL(5,F)
   B=matrix(F,5,5,[1,2,3,1,0,4,5,2,0,1,3,2,5,4,0,1,6,4,3,5,2,0,4,1,1])
   @test v1*B == V([ sum([v1[i]*B[i,j] for i in 1:5]) for j in 1:5 ])
   @test V(transpose(B*v2))==V([ sum([v2[i]*B[j,i] for i in 1:5]) for j in 1:5 ])
   B = G(B)
   @test v1*B == V([ sum([v1[i]*B[i,j] for i in 1:5]) for j in 1:5 ])
   @test V(transpose(B*v2))==V([ sum([v2[i]*B[j,i] for i in 1:5]) for j in 1:5 ])
   @test map(x->x+1,v1)==V([2,3,4,5,6])

   # see the discussion of pull/1368 and issues/872
   @test_throws MethodError v1*v2
end

# from file matrices/stuff_field_gen.jl
@testset "Stuff on fields" begin
   F = GF(3)
   R,t = PolynomialRing(F,"t")
   f = t^2+1
   f1 = Oscar._change_type(f)
   @test collect(coefficients(f1))==collect(coefficients(f))
   F1 = GF(3, 1)
   @test base_ring(f1)==F1
   x = Oscar._centralizer(f1)(companion_matrix(f1))
   @test order(GL(2,F1)(x))==8
   K,z = FiniteField(f,"z")
   @test z^4==1
   @test (change_base_ring(K,Oscar._centralizer(f1))(z))^4 !=1

   F = GF(17, 1)
   a = F(3)
   b = F(13)
   @test a^Oscar.disc_log(a,b)==b
   @test_throws String Oscar.disc_log(b,a)
   @test Oscar.disc_log(F(16),F(1))==0
   @test_throws String Oscar.disc_log(b,F(0))
end
