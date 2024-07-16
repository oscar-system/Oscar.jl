@testset "Definition forms" begin
   T,t = polynomial_ring(GF(3),"t")
   F,z = finite_field(t^2+1,"z")

   B = matrix(F,4,4,[0 1 0 0; 2 0 0 0; 0 0 0 z+2; 0 0 1-z 0])
   @test is_alternating(B)
   f = alternating_form(B)
   @test f isa SesquilinearForm
   @test gram_matrix(f)==B
   @test f==alternating_form(gram_matrix(f))
   @test base_ring(B)==F
   @test !is_symmetric(B)
   @test !is_hermitian(B)
   @test is_alternating(f)
   @test !is_quadratic(f)
   @test !is_symmetric(f)
   @test !is_hermitian(f)
   @test_throws AssertionError f = symmetric_form(B)
   @test_throws AssertionError f = hermitian_form(B)

   B = matrix(F,4,4,[0 1 0 0; 1 0 0 0; 0 0 0 z+2; 0 0 -1-z 0])
   @test is_hermitian(B)
   f = hermitian_form(B)
   @test f isa SesquilinearForm
   @test gram_matrix(f)==B
   @test is_hermitian(f)
   @test GAP.Obj(f) isa GapObj
   @test_throws AssertionError f = symmetric_form(B)
   @test_throws AssertionError f = alternating_form(B)
   @test_throws ArgumentError corresponding_quadratic_form(f)
   @test_throws ErrorException f = SesquilinearForm(B,:unitary)

   B = matrix(F,4,4,[0 1 0 0; 1 0 0 0; 0 0 0 z+2; 0 0 z+2 0])
   @test is_symmetric(B)
   f = symmetric_form(B)
   @test f isa SesquilinearForm
   @test gram_matrix(f)==B
   @test is_symmetric(f)
   @test !is_hermitian(f)
   @test_throws AssertionError f = alternating_form(B)
   Qf = corresponding_quadratic_form(f)
   R = polynomial_ring(F,4)[1]
   p = R[1]*R[2]+(z+2)*R[3]*R[4]
   Q = quadratic_form(p)
   @test Q==Qf
   @test corresponding_quadratic_form(corresponding_bilinear_form(Q))==Q
   @test corresponding_bilinear_form(Q)==f
   B1 = matrix(F,4,4,[0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 z+2 0])
   Q1 = quadratic_form(B1)
   @test Q1==Q
   @test gram_matrix(Q1)!=B1
   @test defining_polynomial(Q1) isa MPolyRingElem
   pf = defining_polynomial(Q1)
   @test defining_polynomial(Q1)==parent(pf)[1]*parent(pf)[2]+(z+2)*parent(pf)[3]*parent(pf)[4]
# I can't test simply pf==p, because it returns FALSE. The line
 #      polynomial_ring(F,4)[1]==polynomial_ring(F,4)[1]
# returns FALSE.
   @test_throws ArgumentError corresponding_quadratic_form(Q)
   @test_throws ArgumentError corresponding_bilinear_form(f)

   R,x = polynomial_ring(F,"x")
   p = x^2*z
   Q = quadratic_form(p)
   @test is_quadratic(Q)
   f = corresponding_bilinear_form(Q)
   @test is_symmetric(f)
   @test gram_matrix(f)==matrix(F,1,1,[-z])

   T,t = polynomial_ring(GF(2),"t")
   F,z = finite_field(t^2+t+1,"z")
   R = polynomial_ring(F,4)[1]
   p = R[1]*R[2]+z*R[3]*R[4]
   Q = quadratic_form(p)
   @test is_quadratic(Q)
   @test gram_matrix(Q)==matrix(F,4,4,[0 1 0 0; 0 0 0 0; 0 0 0 z; 0 0 0 0])
   f = corresponding_bilinear_form(Q)
   @test is_alternating(f)
   @test gram_matrix(f)==matrix(F,4,4,[0 1 0 0; 1 0 0 0; 0 0 0 z; 0 0 z 0])
   @test_throws ArgumentError corresponding_quadratic_form(f)


end

@testset "Evaluating forms" begin
   F,z = finite_field(3,2,"z")
   V=vector_space(F,6)

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

   @test_throws ArgumentError Q(V[1],V[5])
   @test_throws ArgumentError f(V[2])

   g = rand(GL(6,F))
   v = rand(V)
   w = rand(V)
   @test (f^g)(v,w)==f(v*g^-1,w*g^-1)
   @test (Q^g)(v)==Q(v*g^-1)
end

@testset "Methods with forms" begin
   F = GF(5,1)
   V = vector_space(F,6)
   x = zero_matrix(F,6,6)
   f = symmetric_form(x)
   @test radical(f)[1]==sub(V,gens(V))[1]
   x[3,5]=1; x[4,6]=1;
   f = symmetric_form(x+transpose(x))
   @test radical(f)[1]==sub(V,[V[1],V[2]])[1]
   @test f*F(3)==F(3)*f;
   @test is_symmetric(f*F(3))
   @test gram_matrix(f*F(3))==F(3)*(x+transpose(x))
   @test is_degenerate(f)
   x[1,2]=1;

   f = symmetric_form(x+transpose(x))
   @test radical(f)[1]==sub(V,[])[1]
   Q = corresponding_quadratic_form(f)
   @test radical(Q)[1]==radical(f)[1]
   @test witt_index(f)==3
   @test witt_index(Q)==3
   @test !is_degenerate(f)

   F = GF(2,1)
   x = matrix(F,2,2,[1,0,0,0])
   Q = quadratic_form(x)
   @test dim(radical(Q)[1])==1
   f = corresponding_bilinear_form(Q)
   @test dim(radical(f)[1])==2

end

@testset "TransformForm" begin
   # symmetric
   F = GF(3,1)
   x = zero_matrix(F,6,6)
   x[4,5]=1; x[5,6]=1; x=x+transpose(x)
   y = zero_matrix(F,6,6)
   y[1,3]=2;y[3,4]=1; y=y+transpose(y)
   f = symmetric_form(x); g = symmetric_form(y)
   is_true,z = is_congruent(f,g)
   @test is_true
   @test f^z == g

   F = GF(7,1)
   x = diagonal_matrix(F.([1,4,2,3,6,5,4]))
   y = diagonal_matrix(F.([3,1,5,6,4,2,1]))
   f = symmetric_form(x); g = symmetric_form(y)
   is_true,z = is_congruent(f,g)
   @test is_true
   @test f^z == g
   y = diagonal_matrix(F.([3,1,5,6,4,3,1]))
   f = symmetric_form(x); g = symmetric_form(y)
   is_true,z = is_congruent(f,g)
   @test !is_true
   @test z===nothing

   T,t = polynomial_ring(GF(3),"t")
   F,a = finite_field(t^2+1,"a")
   x = zero_matrix(F,6,6)
   x[1,2]=1+2*a; x[3,4]=a; x[5,6]=1; x=x+transpose(x)
   y = diagonal_matrix(F.([a,1,1,a+1,2,2*a+2]))
   f = symmetric_form(x); g = symmetric_form(y)
   is_true,z = is_congruent(f,g)
   @test is_true
   @test f^z == g
   y = diagonal_matrix(F.([a,1,1,a+1,2,2*a]))
   g = symmetric_form(y)
   is_true,z = is_congruent(f,g)
   @test !is_true

   #alternating
   F = GF(3,1)
   x = zero_matrix(F,6,6)
   x[4,5]=1; x[5,6]=1; x=x-transpose(x)
   y = zero_matrix(F,6,6)
   y[1,3]=2;y[3,4]=1; y=y-transpose(y)
   f = alternating_form(x); g = alternating_form(y)
   is_true,z = is_congruent(f,g)
   @test is_true
   @test f^z == g

   F,a = finite_field(2,3,"a")
   x = zero_matrix(F,6,6)
   x[1,2]=a; x[2,3]=a^2+1; x[3,4]=1; x[1,5]=a^2+a+1; x[5,6]=1; x=x-transpose(x)
   y = zero_matrix(F,6,6)
   y[1,6]=1; y[2,5]=a; y[3,4]=a^2+1; y = y-transpose(y)
   f = alternating_form(x); g = alternating_form(y)
   is_true,z = is_congruent(f,g)
   @test is_true
   @test f^z == g
   x = zero_matrix(F,8,8)
   x[1,2]=a; x[2,3]=a^2-1; x[3,4]=1; x[1,5]=a^2-a-1; x[5,6]=1; x[7,8]=-1; x=x-transpose(x)
   y = zero_matrix(F,8,8)
   y[1,8]=1; y[2,7]=a; y[3,6]=a^2+1; y[4,5] = -a^2-a-1; y = y-transpose(y)
   f = alternating_form(x); g = alternating_form(y)
   is_true,z = is_congruent(f,g)
   @test is_true
   @test f^z == g
   y = zero_matrix(F,8,8)
   y[1,8]=1; y[2,7]=a; y[3,6]=a^2+1; y = y-transpose(y)
   f = alternating_form(x); g = alternating_form(y)
   is_true,z = is_congruent(f,g)
   @test !is_true

   y = zero_matrix(F,6,6)
   y[1,6]=1; y[2,5]=a; y[3,4]=a^2+1; y = y-transpose(y)
   g = alternating_form(y)
   @test_throws ArgumentError is_congruent(f,g)

   F = GF(3,1)
   y = zero_matrix(F,8,8)
   y[1,3]=2;y[3,4]=1; y=y-transpose(y)
   g = alternating_form(y)
   @test_throws ArgumentError is_congruent(f,g)

   #hermitian
   F,a = finite_field(3,2,"a")
   x = zero_matrix(F,6,6)
   x[4,5]=1; x[5,6]=a-1; x=x+conjugate_transpose(x)
   y = zero_matrix(F,6,6)
   y[1,2]=2; y=y+conjugate_transpose(y)
   f = hermitian_form(x); g = hermitian_form(y)
   is_true,z = is_congruent(f,g)
   @test is_true
   @test f^z == g
   x = diagonal_matrix(F.([1,2,1,1,1]))
   y = diagonal_matrix(F.([2,1,0,0,2]))
   y[3,4]=a; y[4,3]=a^3; y[4,5]=1+a; y[5,4]=(1+a)^3;
   f = hermitian_form(x); g = hermitian_form(y)
   is_true,z = is_congruent(f,g)
   @test is_true
   @test f^z == g

   F,a = finite_field(2,2,"a")
   x = zero_matrix(F,6,6)
   x[4,5]=1; x[5,6]=a+1; x=x+conjugate_transpose(x)
   y = zero_matrix(F,6,6)
   y[1,2]=1; y=y+conjugate_transpose(y)
   f = hermitian_form(x); g = hermitian_form(y)
   is_true,z = is_congruent(f,g)
   @test is_true
   @test f^z == g
   x = diagonal_matrix(F.([1,1,1,1,1]))
   y = diagonal_matrix(F.([1,1,0,0,1]))
   y[3,4]=a; y[4,3]=a^2; y[4,5]=1+a; y[5,4]=(1+a)^2;
   f = hermitian_form(x); g = hermitian_form(y)
   is_true,z = is_congruent(f,g)
   @test is_true
   @test f^z == g

   #quadratic
   F = GF(5,1)
   R = polynomial_ring(F,6)[1]
   p1 = R[1]*R[2]
   p2 = R[4]*R[5]+3*R[4]^2
   Q1 = quadratic_form(p1)
   Q2 = quadratic_form(p2)
   is_true,z = is_congruent(Q1,Q2)
   @test is_true
   @test Q1^z == Q2
   p2 = R[4]*R[5]+3*R[4]^2+2*R[5]^2
   Q2 = quadratic_form(p2)
   is_true,z = is_congruent(Q1,Q2)
   @test !is_true

   p1 = R[1]^2+R[2]^2+R[3]^2+R[4]^2+R[5]*R[6]
   p2 = R[1]*R[6]+R[2]*R[5]+R[3]*R[4]
   Q1 = quadratic_form(p1)
   Q2 = quadratic_form(p2)
   is_true,z = is_congruent(Q1,Q2)
   @test is_true
   @test Q1^z == Q2
   p1 = R[1]^2+R[2]^2+R[3]^2+R[4]^2+R[5]^2+R[5]*R[6]+2*R[6]^2
   Q1 = quadratic_form(p1)
   is_true,z = is_congruent(Q1,Q2)
   @test !is_true

   F,a = finite_field(2,2,"a")
   R = polynomial_ring(F,6)[1]
   p1 = R[1]*R[2]+R[3]*R[4]+R[5]^2+R[5]*R[6]+R[6]^2
   p2 = R[1]*R[6]+a*R[2]*R[5]+R[3]*R[4]
   Q1 = quadratic_form(p1)
   Q2 = quadratic_form(p2)
   is_true,z = is_congruent(Q1,Q2)
   @test is_true
   @test Q1^z == Q2
   p1 = R[1]*R[2]+R[3]*R[4]+R[5]^2+R[5]*R[6]+a*R[6]^2
   Q1 = quadratic_form(p1)
   is_true,z = is_congruent(Q1,Q2)
   @test !is_true
   p1 = R[1]*R[2]
   p2 = R[3]^2+(a+1)*R[3]*R[5]
   Q1 = quadratic_form(p1)
   Q2 = quadratic_form(p2)
   is_true,z = is_congruent(Q1,Q2)
   @test is_true
   @test Q1^z == Q2
   p2 = R[3]^2+R[3]*R[5]+(a+1)*R[5]^2
   Q2 = quadratic_form(p2)
   is_true,z = is_congruent(Q1,Q2)
   @test !is_true
end

# TODO
@testset "Relationship group - forms" begin
   G = GL(6,3)
   Op = GO(1,6,3)
   Om = GO(-1,6,3)
   F = base_ring(G)
   B = zero_matrix(F,6,6)
   for i in 1:6 B[i,7-i]=1 end
   f = symmetric_form(B)
   H = isometry_group(f)
   @testset for x in gens(H)
      @test f^x==f
   end
   @test order(H)==order(Op)
   B[3:4,3:4] = identity_matrix(F,2)
   f = symmetric_form(B)
   H = isometry_group(f)
   @testset for x in gens(H)
      @test f^x==f
   end
   @test order(H)==order(Om)

   G = GL(4,2)
   Op = GO(1,4,2)
   Om = GO(-1,4,2)
   F = base_ring(G)
   B = zero_matrix(F,4,4)
   B[1,4]=1; B[2,3]=1
   f = quadratic_form(B)
   H = isometry_group(f)
   @testset for x in gens(H)
      @test f^x==f
   end
   @test is_conjugate(G, H, Op)
   B[2,2]=1; B[3,3]=1;
   f = quadratic_form(B)
   H = isometry_group(f)
   @testset for x in gens(H)
      @test f^x==f
   end
   @test is_conjugate( G, H, Om)
   Q = f
   f = corresponding_bilinear_form(Q)
   H = isometry_group(f)
   @testset for x in gens(H)
      @test f^x==f
   end
   @test is_conjugate(G, H, Sp(4,2))

   G = GL(5,9)
   F = base_ring(G)
   B = diagonal_matrix(F.([1,1,1,2,2]))
   B[4,5] = gen(F)
   B[5,4] = gen(F)^3
   f = hermitian_form(B)
   H = isometry_group(f)
   @testset for x in gens(H)
      @test f^x==f
   end
   @test order(H)==order(GU(5,3))

   G = GU(2,3)
   L = invariant_sesquilinear_forms(G)
   @testset for f in L
       for g in gens(G)
          @test matrix(g)*f*conjugate_transpose(matrix(g))==f
       end
   end
   G = GO(-1,4,3)
   L = invariant_bilinear_forms(G)
   @testset for f in L
       for g in gens(G)
          @test matrix(g)*f*transpose(matrix(g))==f
       end
   end
   L = invariant_quadratic_forms(G)
   @testset for m in L
       Q = quadratic_form(m)
       for g in gens(G)
          @test Q^g==Q
       end
   end

   G = GO(0,5,3)
   L = preserved_quadratic_forms(G)
   @testset for f in L
       for g in gens(G)
          @test f^g==f
       end
   end
   L = invariant_quadratic_forms(G)
   L = [quadratic_form(f) for f in L]
   @testset for f in L
       for g in gens(G)
          @test f^g==f
       end
   end

   G = GO(1,6,4)
   L = preserved_quadratic_forms(G)
   @testset for f in L
       for g in gens(G)
          @test f^g==f
       end
   end
   L = invariant_quadratic_forms(G)
   L = [quadratic_form(f) for f in L]
   @testset for f in L
       for g in gens(G)
          @test f^g==f
       end
   end
   Q = quadratic_form(Oscar.invariant_quadratic_form(G))
   @testset for g in gens(G)
      Q^g==Q
   end
   

   G = Sp(4,4)
   L = preserved_sesquilinear_forms(G)
   @testset for f in L
      for g in gens(G)
         @test f^g==f
      end
   end
   B = Oscar.invariant_bilinear_form(G)
   @testset for g in gens(G)
      @test matrix(g)*B*transpose(matrix(g))==B
   end
   L = preserved_sesquilinear_forms(G)
   @testset for f in L
      for g in gens(G)
         @test f^g==f
      end
   end
   

   G = GO(-1,6,3)
   L = preserved_sesquilinear_forms(G)
   @testset for f in L
      for g in gens(G)
         @test f^g==f
      end
   end
   B = Oscar.invariant_bilinear_form(G)
   @testset for g in gens(G)
      @test matrix(g)*B*transpose(matrix(g))==B
   end
   B = Oscar.invariant_quadratic_form(G)
   @testset for g in gens(G)
      @test is_alternating(matrix(g)*B*transpose(matrix(g))-B)
   end

   G = GU(4,5)
   L = preserved_sesquilinear_forms(G)
   @testset for f in L
      for g in gens(G)
         @test f^g==f
      end
   end
   B = Oscar.invariant_sesquilinear_form(G)
   @testset for g in gens(G)
      @test matrix(g)*B*conjugate_transpose(matrix(g))==B
   end

   G = general_linear_group(2, 3)
   @test_throws ArgumentError Oscar.invariant_sesquilinear_form(G)

   @testset for q in [2,3,4,5]
      G = GU(5,q)
      x = rand(GL(5,base_ring(G)))
      G = G^x
      L = invariant_hermitian_forms(G)
      @testset for m in L
         f = hermitian_form(m)
         for g in gens(G)
            @test f^g==f
         end
      end
   end
   @test length(invariant_hermitian_forms(GL(2,5)))==0

   @testset for q in [3,4,5]
      G = Sp(6,q)^rand(GL(6,q))
      L = invariant_alternating_forms(G)
      @testset for m in L
         f = alternating_form(m)
         for g in gens(G)
            @test f^g==f
         end
      end
   end

   G = GO(1,6,3)^rand(GL(6,3))
   L = invariant_symmetric_forms(G)
   for m in L
      f = symmetric_form(m)
      for g in gens(G)
         @test f^g==f
      end
   end
   G = omega_group(-1,6,9)^rand(GL(6,9))
   L = invariant_symmetric_forms(G)
   for m in L
      f = symmetric_form(m)
      for g in gens(G)
         @test f^g==f
      end
   end
   G = SO(0,5,7)^rand(GL(5,7))
   L = invariant_symmetric_forms(G)
   for m in L
      f = symmetric_form(m)
      for g in gens(G)
         @test f^g==f
      end
   end

   G = GO(1,6,8)^rand(GL(6,8))
   L = invariant_quadratic_forms(G)
   for m in L
      f = quadratic_form(m)
      for g in gens(G)
         @test f^g==f
      end
   end
   G = omega_group(-1,6,2)^rand(GL(6,2))
   L = invariant_quadratic_forms(G)
   for m in L
      f = quadratic_form(m)
      for g in gens(G)
         @test f^g==f
      end
   end
   G = GO(0,5,4)^rand(GL(5,4))
   L = invariant_quadratic_forms(G)
   for m in L
      f = quadratic_form(m)
      for g in gens(G)
         @test f^g==f
      end
   end

   G = sub(GL(6,3),[one(GL(6,3))])[1]
   @test length(invariant_symmetric_forms(G))==21
   @test length(invariant_alternating_forms(G))==15
   @test length(invariant_quadratic_forms(G))==21

   # degenerate forms
   F = GF(3,1)
   B = diagonal_matrix(F.([1,1,1,0,0,0]))
   Q = quadratic_form(B)
   H = isometry_group(Q)
   @testset for h in gens(H)
      @test Q^h==Q
   end
   @test order(H)==order(GO(3,3))*order(GL(3,F))*3^9

   F = GF(2,1)
   B = zero_matrix(F,6,6)
   B[1,2]=1;B[3,4]=1;
   Q = quadratic_form(B)
   H = isometry_group(Q)
   @testset for h in gens(H)
      @test Q^h==Q
   end
   @test order(H)==order(GO(1,4,2))*order(GL(2,F))*2^8

   F = GF(3,2)
   B = zero_matrix(F,4,4)
   B[3,4]=gen(F);B[4,3]=gen(F)^3
   f = hermitian_form(B)
   H = isometry_group(f)
   @testset for h in gens(H)
      @test f^h==f
   end
   @test order(H)==order(GU(2,3))*order(GL(2,9))*9^4

   @testset "orthogonal sign" begin
      @test orthogonal_sign(orthogonal_group(+1, 4, 2)) == +1
      @test orthogonal_sign(orthogonal_group(-1, 4, 2)) == -1
      @test orthogonal_sign(orthogonal_group(+1, 4, 3)) == +1
      @test orthogonal_sign(orthogonal_group(-1, 4, 3)) == -1
      @test orthogonal_sign(orthogonal_group(0, 5, 3)) == 0
      # Odd dimensional orthogonal groups in char. 2 are not irreducible.
      @test_throws ErrorException orthogonal_sign(orthogonal_group(0, 5, 2))
      @test orthogonal_sign(general_linear_group(4, 2)) === nothing
      # If the abs. irred. module preserves an antisymmetric invariant
      # bilinear form then there is no nondegenerate quadratic form.
      F = GF(7)
      mats = [matrix(F, [0 1 0 0 0 0; 1 0 0 0 0 0; 0 0 0 1 0 0;
                         0 0 1 0 0 0; 5 5 2 2 6 0; 3 3 4 4 0 6]),
              matrix(F, [4 0 0 0 0 0; 0 0 1 0 0 0; 0 6 1 0 0 0;
                         0 0 0 0 1 0; 0 0 0 0 0 1; 6 0 0 2 4 3])]
      G = matrix_group(mats)
      @test describe(G) == "PSU(3,3)"
      @test orthogonal_sign(G) === nothing
   end
end

@testset "Orthogonal groups of ZZ-lattices" begin
  N1 = root_lattice(:A, 2)
  N2 = rescale(N1, 4)
  N,_ = direct_sum(N1,N2)
  @test order(orthogonal_group(N))==144

  L = integer_lattice(; gram=QQ[4 0 0 0 0; 0 16 4 10 8; 0 4 2 3 2; 0 10 3 10 5; 0 8 2 5 34])
  G = orthogonal_group(L)
  @test order(G)==32
  @test order(Oscar._isometry_group_via_decomposition(L, closed=false)[1]) == 32
  @test order(Oscar._isometry_group_via_decomposition(L, closed=false, direct=false)[1]) == 32
  @test order(Oscar._isometry_group_via_decomposition(L, closed=true, direct=false)[1]) == 32
  L = integer_lattice(; gram=QQ[4 0 0 0 0; 0 16 4 10 8; 0 4 2 3 2; 0 10 3 10 5; 0 8 2 5 34]) # avoid caching
  @test order(isometry_group(L, algorithm = :decomposition, depth = 1, bacher_depth = 0)) == 32

  gram = ZZ[2 1 -1 -1 -1 1 1 -1 0 0 0 0 0 0 0 0; 1 2 -1 -1 -1 1 1 -1 0 0 0 0 0 0 0 0; -1 -1 2 0 1 0 -1 1 0 0 0 0 0 0 0 0; -1 -1 0 2 1 -1 0 0 0 0 0 0 0 0 0 0; -1 -1 1 1 2 0 -1 0 0 0 0 0 0 0 0 0; 1 1 0 -1 0 2 0 -1 0 0 0 0 0 0 0 0; 1 1 -1 0 -1 0 2 -1 0 0 0 0 0 0 0 0; -1 -1 1 0 0 -1 -1 2 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 2 1 1 0 1 1 1 0; 0 0 0 0 0 0 0 0 1 2 1 0 1 1 0 0; 0 0 0 0 0 0 0 0 1 1 2 0 0 0 1 0; 0 0 0 0 0 0 0 0 0 0 0 2 1 0 -1 0; 0 0 0 0 0 0 0 0 1 1 0 1 4 1 0 1; 0 0 0 0 0 0 0 0 1 1 0 0 1 4 0 0; 0 0 0 0 0 0 0 0 1 0 1 -1 0 0 8 1; 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 18]
  L = integer_lattice(; gram)
  @test order(orthogonal_group(L; algorithm = :decomposition)) == 267544166400

  H = hyperbolic_plane_lattice()
  G = orthogonal_group(H)
  @test order(G) == 4
end
