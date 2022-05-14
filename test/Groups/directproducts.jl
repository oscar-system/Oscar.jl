@testset "Directproducts" begin
   S = symmetric_group(4)
   C = abelian_group(PcGroup, [2,2])
   G = direct_product(S,C)

   @test G isa DirectProductGroup
   @test order(G)==order(S)*order(C)
   @test exponent(G)==lcm(exponent(S),exponent(C))
   @test rand(G) isa elem_type(DirectProductGroup)
   @test factor_of_direct_product(G,1)==S
   @test factor_of_direct_product(G,2)==C
   @test_throws ArgumentError factor_of_direct_product(G,3)
   @test number_of_factors(G)==2
   @test isfull_direct_product(G)
   @test PermGroup(G) isa PermGroup

   G,emb,proj = direct_product(S,C; morphisms=true)
   @test G==direct_product(S,C)
   @test emb[1]==embedding(G,1)
   @test emb[2]==embedding(G,2)
   @test proj[1]==projection(G,1)
   @test proj[2]==projection(G,2)

   x = rand(G)
   @test x in G
   @test projection(G,1)(x) in S
   @test projection(G,2)(x) in C
   @test embedding(G,1)(rand(S)) in G
   @test embedding(G,2)(rand(C)) in G
   @test x==G(projection(G,1)(x), projection(G,2)(x))
   @test x==embedding(G,1)(projection(G,1)(x))*embedding(G,2)(projection(G,2)(x))
   S1 = image(embedding(G,1))[1]
   C1 = image(embedding(G,2))[1]
   @test intersect(G,S1)[1]==S1
   @test is_isomorphic(S1,S)
   @test is_isomorphic(quo(G,S1)[1],C)
   @test is_isomorphic(quo(G,C1)[1],S1)

   x = G(cperm([1,2]),C[1])
   @test x==G([cperm([1,2]),C[1]])
   H = sub(G,[x])[1]
   @test H==sub(x,x)[1]
   @test x==H(cperm([1,2]),C[1])
   @test_throws ArgumentError H(cperm([2,3]),C[1])
   @test order(H)==2
   @test index(G,H) isa Oscar.fmpz
   @test_throws ArgumentError write_as_full(H)
   @test G==write_as_full(G)
   @test intersect(G,H)[1]==H
   x = G(cperm([1,2,3]),C[1])
   H = sub(G,[x])[1]
   @test order(H)==6
   @test H==write_as_full(H)
   @test intersect(G,H)[1]==H

   P1=sylow_subgroup(S,2)[1]
   P2=sylow_subgroup(C,2)[1]
   P=sylow_subgroup(G,2)[1]
   PP=direct_product(P1,P2)
   @test order(P)==order(P1)*order(P2)
   @test representative_action(G,P,PP)[1]
   x = representative_action(G,P,PP)[2]
   @test P^x==PP

   @test_throws ArgumentError as_perm_group(G)
   @test_throws ArgumentError as_polycyclic_group(G)
   G = direct_product(S,S)
   @test G isa DirectProductGroup
   @test as_perm_group(G) isa PermGroup
   G = direct_product(C,C)
   @test as_polycyclic_group(G) isa PcGroup

   @testset "Inner direct products" begin
      C2 = cyclic_group(2)
      C3 = cyclic_group(3)
      G = inner_direct_product(C2,C2,C3,C3)
      @test G isa PcGroup
      @test is_isomorphic(G,direct_product(C2,C2,C3,C3))
      A = alternating_group(3)
      G = inner_direct_product(A,A,A)
      @test G isa PermGroup
      G1 = G
      @test degree(G)==9
      L = [A,A,A]
      @test G==inner_direct_product(L)

      G,Le,Lp = inner_direct_product(A,A,A; morphisms=true)
      @test G==G1
      @test is_injective(Le[1])
      @test is_surjective(Lp[1])
      @testset for i in 1:3
         @test domain(Le[i])==L[i]
         @test codomain(Le[i])==G
         @test domain(Lp[i])==G
         @test codomain(Lp[i])==L[i]
         @test Le[i]*Lp[i]==id_hom(L[i])
      end

      G1,Le1,Lp1 = inner_cartesian_power(A,3; morphisms=true)         
      @test G1==G
      @test Le1==Le
      @test Lp1==Lp
   end
      


   @testset "Cartesian Power" begin
      C = cyclic_group(3)
      G = cartesian_power(C,5)
      @test G isa DirectProductGroup
      @test order(G)==243
      @test number_of_factors(G)==5
      @test is_abelian(G)
      @test factor_of_direct_product(G,5)==C
      @test is_isomorphic(G,abelian_group(PcGroup,[3,3,3,3,3]))
      @test PermGroup(G) isa PermGroup
      x1 = G(C[1],one(C),one(C),one(C),one(C))
      x2 = G(one(C),C[1],one(C),one(C),one(C))
      x3 = G(one(C),one(C),C[1],one(C),one(C))
      H = sub(G,[x1,x2,x3])[1]
      K = sub(G,[x1*x2*x3])[1]
      @test !isfull_direct_product(H)
      @test !isfull_direct_product(K)
      Hf = write_as_full(H)
      @test isfull_direct_product(Hf)
      @test_throws ArgumentError write_as_full(K)
   end
end

@testset "Semidirectproducts" begin
   Q=quaternion_group(8)
   C=cyclic_group(2)
   A=automorphism_group(Q)
   au=A(hom(Q,Q,[Q[1],Q[2]],[Q[1]^3,Q[2]^3]))
   f = hom(C,A,[C[1]],[au])

   G = semidirect_product(Q,f,C)
   @test G isa SemidirectProductGroup{PcGroup,PcGroup}
   @test normal_subgroup(G)==Q
   @test is_isomorphic(acting_subgroup(G),C)
   @test homomorphism_of_semidirect_product(G)==f
   @test order(G)==16
   @test isfull_semidirect_product(G)
   @test PermGroup(G) isa PermGroup
   x = G(Q[1]*Q[2],C[1])
   @test parent(x)==G
   H = sub(G,[x])[1]
   @test H(Q[1]*Q[2],C[1])==x
   @test_throws ArgumentError H(Q[1],C[1])
   @test is_subgroup(G,H)[1]
   @test index(G,H)==4
   @test !isfull_semidirect_product(H)
   @test projection(G)(x)==projection(H)(x)
   @test H==center(G)[1]
   y=G(Q[1],one(C))
   K = sub(G,gens(G))[1]
#   @test isfull_semidirect_product(K)
   K = sub(G,[y])[1]
   @test K==sub(y)[1]
   @test K==sub(y)[1]
   @test y == embedding(G,1)(Q[1])
   @test_throws ArgumentError embedding(G,3)(Q[1])
   @test codomain(projection(K))==C
   @test order(image(projection(K))[1])==1
   @test embedding(G,2)*projection(G)==id_hom(C)
   @test image(embedding(G,1))[1]==kernel(projection(G))[1]
end


@testset "Wreathproducts" begin
   C = cyclic_group(2)
   H = symmetric_group(3)
   W = wreath_product(C,H)
   @test W isa WreathProductGroup
   @test order(W)==48
   @test isfull_wreath_product(W)
   @test W==sub(W,gens(W))[1]
   @test PermGroup(W) isa PermGroup

   H = sub(cperm([1,2,4]))[1]
   W = wreath_product(C,H)

   @test W isa WreathProductGroup
   @test order(W)==2^4*3
   @test !is_abelian(W)
   @test rand(W) isa elem_type(WreathProductGroup)
   f1 = C[1]
   x = W(f1,one(C),f1,one(C),cperm([1,4,2]))
   @test embedding(W,1)(f1)==W(f1,one(C),one(C),one(C),one(H))
   @test embedding(W,4)(f1)==W(one(C),one(C),one(C),f1,one(H))
   @test_throws ArgumentError embedding(W,7)(f1) in W
   @test embedding(W,5)(cperm([1,4,2]))==W(one(C),one(C),one(C),one(C),cperm([1,4,2]))
   @test projection(W)(x)==cperm([1,4,2])
   @test codomain(projection(W))==H
   @test domain(embedding(W,2))==C
   K = sub(W,gens(W))[1]
#   @test isfull_wreath_product(K)
   K = sub(W,[x])[1]
   @test K==sub(x)[1]
   @test K==sub(x)[1]
   @test K isa WreathProductGroup
   @test order(K)==6
   @test is_cyclic(K)
   @test index(W,K)==8
   @test_throws ArgumentError embedding(K,1)(f1) in K

   C = cyclic_group(3)
   a = hom(C,H,[C[1]],[cperm([1,2,4])])
   W = wreath_product(C,C,a)
   @test W isa WreathProductGroup
   @test order(W)==81
   @test is_isomorphic(normal_subgroup(W),C)
   @test is_isomorphic(acting_subgroup(W),C)
   @test homomorphism_of_wreath_product(W)==a
   @test embedding(W,1)(C[1]) in W
   @test W(C[1],one(C),one(C),C[1]) isa elem_type(WreathProductGroup)
   @test image(projection(W))[1]==C
   @test_throws ArgumentError embedding(W,5)(C[1])
end
