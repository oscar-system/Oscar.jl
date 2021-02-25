@testset "Definition forms" begin
   T,t = PolynomialRing(FiniteField(3),"t")
   F,z = FiniteField(t^2+1,"z")

   B = matrix(F,4,4,[0 1 0 0; 2 0 0 0; 0 0 0 z+2; 0 0 1-z 0])
   @test isskewsymmetric_matrix(B)
   f = alternating_form(B)
   @test f isa SesquilinearForm
   @test gram_matrix(f)==B
   @test f==alternating_form(gram_matrix(f))
   @test base_ring(B)==F
   @test !issymmetric(B)
   @test !ishermitian_matrix(B)
   @test isalternating_form(f)
   @test !isquadratic_form(f)
   @test !issymmetric_form(f)
   @test !ishermitian_form(f)
   @test_throws AssertionError f = symmetric_form(B)
   @test_throws AssertionError f = hermitian_form(B)
   @test hermitian_form(B; check=false) isa SesquilinearForm

   B = matrix(F,4,4,[0 1 0 0; 1 0 0 0; 0 0 0 z+2; 0 0 -1-z 0])
   @test ishermitian_matrix(B)
   f = hermitian_form(B)
   @test f isa SesquilinearForm
   @test gram_matrix(f)==B
   @test ishermitian_form(f)
   @test f.mat_iso isa Oscar.GenMatIso
   @test f.X isa GapObj
   @test_throws AssertionError f = symmetric_form(B)
   @test_throws AssertionError f = alternating_form(B)
   @test_throws ArgumentError corresponding_quadratic_form(f)
   f = SesquilinearForm(B,:unitary)
   @test_throws ErrorException f.X

   B = matrix(F,4,4,[0 1 0 0; 1 0 0 0; 0 0 0 z+2; 0 0 z+2 0])
   @test issymmetric(B)
   f = symmetric_form(B)
   @test f isa SesquilinearForm
   @test gram_matrix(f)==B
   @test issymmetric_form(f)
   @test !ishermitian_form(f)
   @test_throws AssertionError f = alternating_form(B)
   Qf = corresponding_quadratic_form(f)
   R = PolynomialRing(F,4)[1]
   p = R[1]*R[2]+(z+2)*R[3]*R[4]
   Q = quadratic_form(p)
   @test Q==Qf
   @test corresponding_quadratic_form(corresponding_bilinear_form(Q))==Q
   @test corresponding_bilinear_form(Q)==f
   B1 = matrix(F,4,4,[0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 z+2 0])
   Q1 = quadratic_form(B1)
   @test Q1==Q
   @test gram_matrix(Q1)!=B1
   @test defining_polynomial(Q1) isa AbstractAlgebra.Generic.MPolyElem
   pf = defining_polynomial(Q1)
   @test defining_polynomial(Q1)==parent(pf)[1]*parent(pf)[2]+(z+2)*parent(pf)[3]*parent(pf)[4]
# I can't test simply pf==p, because it returns FALSE. The line
 #      PolynomialRing(F,4)[1]==PolynomialRing(F,4)[1]
# returns FALSE.
   @test_throws ArgumentError corresponding_quadratic_form(Q)
   @test_throws ArgumentError corresponding_bilinear_form(f)

   R,x = PolynomialRing(F,"x")
   p = x^2*z
   Q = quadratic_form(p)
   @test isquadratic_form(Q)
   f = corresponding_bilinear_form(Q)
   @test issymmetric_form(f)
   @test gram_matrix(f)==matrix(F,1,1,[-z])

   T,t = PolynomialRing(FiniteField(2),"t")
   F,z = FiniteField(t^2+t+1,"z")
   R = PolynomialRing(F,4)[1]
   p = R[1]*R[2]+z*R[3]*R[4]
   Q = quadratic_form(p)
   @test isquadratic_form(Q)
   @test gram_matrix(Q)==matrix(F,4,4,[0 1 0 0; 0 0 0 0; 0 0 0 z; 0 0 0 0])
   f = corresponding_bilinear_form(Q)
   @test isalternating_form(f)
   @test gram_matrix(f)==matrix(F,4,4,[0 1 0 0; 1 0 0 0; 0 0 0 z; 0 0 z 0])
   @test_throws ArgumentError corresponding_quadratic_form(f)


end

@testset "Evaluating forms" begin
   F,z=GF(3,2,"z")
   V=VectorSpace(F,6)

   x = matrix(F,6,6,[1,0,0,0,z+1,0,0,0,0,2,1+2*z,1,0,0,1,0,0,z,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
   f = alternating_form(x-transpose(x))
   for i in 1:6
   for j in i:6
      @test f(V[i],V[j])==f.matrix[i,j]
   end
   end

   f = hermitian_form(x+conjugate_transpose(x))
   for i in 1:6
   for j in i:6
      @test f(V[i],V[j])==f(V[j],V[i])^3         # x -> x^3 is the field automorphism
   end
   end

   Q = quadratic_form(x)
   f = corresponding_bilinear_form(Q)
   @test f.matrix==x+transpose(x)
   for i in 1:6
   for j in i:6
      @test f(V[i],V[j])==Q(V[i]+V[j])-Q(V[i])-Q(V[j])
   end
   end

   @test_throws AssertionError Q(V[1],V[5])
   @test_throws AssertionError f(V[2])

   g = rand(GL(6,F))
   v = rand(V)
   w = rand(V)
   @test (f^g)(v,w)==f(v*g^-1,w*g^-1)
   @test (Q^g)(v)==Q(v*g^-1)
end

@testset "Methods with forms" begin
   F = GF(5,1)[1]
   V = VectorSpace(F,6)
   x = zero_matrix(F,6,6)
   f = symmetric_form(x)
   @test radical(f)[1]==sub(V,gens(V))[1]
   x[3,5]=1; x[4,6]=1;
   f = symmetric_form(x+transpose(x))
   @test radical(f)[1]==sub(V,[V[1],V[2]])[1]
   @test f*F(3)==F(3)*f;
   @test issymmetric_form(f*F(3))
   @test gram_matrix(f*F(3))==F(3)*(x+transpose(x))
   @test isdegenerate(f)
   x[1,2]=1;

   f = symmetric_form(x+transpose(x))
   @test radical(f)[1]==sub(V,[])[1]
   Q = corresponding_quadratic_form(f)
   @test radical(Q)[1]==radical(f)[1]
   @test witt_index(f)==3
   @test witt_index(Q)==3
   @test !isdegenerate(f)

   F = GF(2,1)[1]
   x = matrix(F,2,2,[1,0,0,0])
   Q = quadratic_form(x)
   @test dim(radical(Q)[1])==1
   f = corresponding_bilinear_form(Q)
   @test dim(radical(f)[1])==2

end
