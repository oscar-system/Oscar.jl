@testset "Directproducts" begin
   S = symmetric_group(4)
   C = abelian_group(PcGroup, [2,2])
   G = direct_product(S,C)

   @test G isa DirectProductOfGroups
   @test order(G)==order(S)*order(C)
   @test exponent(G)==lcm(exponent(S),exponent(C))
   @test typeof(rand(G))==Oscar.GAPGroupElem{DirectProductOfGroups}
   @test factorofdirectproduct(G,1)==S
   @test factorofdirectproduct(G,2)==C
   @test number_of_factors(G)==2
   @test isfull_direct_product(G)

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
   @test isisomorphic(S1,S)[1]
   @test isisomorphic(quo(G,S1)[1],C)[1]
   @test isisomorphic(quo(G,C1)[1],S1)[1]

   x = G(cperm([1,2]),C[1])
   H = sub(G,[x])[1]
   @test order(H)==2
   @test index(G,H) isa Integer
   @test_throws ArgumentError write_as_full(H)
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
   @test isconjugate(G,P,PP)[1]
   x = isconjugate(G,P,PP)[2]
   @test P^x==PP

   @test_throws ArgumentError as_perm_group(G)
   G = direct_product(S,S)
   @test G isa DirectProductOfGroups
   @test as_perm_group(G) isa PermGroup
   G = direct_product(C,C)
   @test as_polycyclic_group(G) isa PcGroup
   S=SL(2,3)
   G = as_matrix_group(direct_product(S,S))
   @test G isa MatrixGroup
   @test dim(G)==4

   @testset "Cartesian Power" begin
      C = cyclic_group(3)
      G = cartesian_power(C,5)
      @test G isa DirectProductOfGroups
      @test order(G)==243
      @test number_of_factors(G)==5
      @test isabelian(G)
      @test factorofdirectproduct(G,5)==C
      @test isisomorphic(G,abelian_group(PcGroup,[3,3,3,3,3]))[1]
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
   @test G isa SemidirectProductOfGroups{PcGroup,PcGroup}
   @test homomorphism_of_semidirect_product(G)==f
   @test order(G)==16
   @test isfull_semidirect_product(G)
   x = G(Q[1]*Q[2],C[1])
   @test parent(x)==G
   H = sub(G,[x])[1]
   @test issubgroup(G,H)[1]
   @test index(G,H)==4
   @test !isfull_semidirect_product(H)
   @test projection(G)(x)==projection(H)(x)
   @test H==centre(G)[1]
   y=G(Q[1],one(C))
   K = sub(G,[y])[1]
   @test y == embedding(G,1)(Q[1])
   @test_throws ArgumentError embedding(G,3)(Q[1])
   @test codomain(projection(K))==C
   @test order(image(projection(K))[1])==1
   @test embedding(G,2)*projection(G)==id_hom(C)
   @test image(embedding(G,1))[1]==kernel(projection(G))[1]
end


@testset "Wreathproducts" begin
   C = cyclic_group(2)
   H = sub(cperm([1,2,4]))[1]
   W = wreath_product(C,H)

   @test typeof(W)==WreathProduct
   @test order(W)==2^4*3
   @test !isabelian(W)
   @test typeof(rand(W))==Oscar.GAPGroupElem{WreathProduct}
   f1 = C[1]
   x = W(f1,one(C),f1,one(C),cperm([1,4,2]))
   @test embedding(W,1)(f1)==W(f1,one(C),one(C),one(C),one(H))
   @test embedding(W,4)(f1)==W(one(C),one(C),one(C),f1,one(H))
   @test_throws ArgumentError embedding(W,7)(f1) in W
   @test embedding(W,5)(cperm([1,4,2]))==W(one(C),one(C),one(C),one(C),cperm([1,4,2]))
   @test projection(W)(x)==cperm([1,4,2])
   @test codomain(projection(W))==H
   @test domain(embedding(W,2))==C
   K = sub(W,[x])[1]
   @test typeof(K)==WreathProduct
   @test order(K)==6
   @test iscyclic(K)
   @test index(W,K)==8
   @test_throws ArgumentError embedding(K,1)(f1) in K
end
