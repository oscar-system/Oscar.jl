L = [ alternating_group(5), cyclic_group(18), SL(3,3), free_group(0), free_group(1), free_group(2) ]

import Oscar.AbstractAlgebra.GroupsCore

#include(joinpath(dirname(pathof(GroupsCore)), "..", "test", "conformance_test.jl"))

@testset "GAPGroups_interface_conformance for $(G)" for G in L

   # TODO: enable GroupsCore conformance tests; this requires a new GroupsCore
   # release with https://github.com/kalmarek/GroupsCore.jl/pull/41 merged.
   #test_Group_interface(G)
   #test_GroupElement_interface(rand(G, 2)...)

   g, h = rand(G,2)

   @test parent(g) isa typeof(G)
   @test parent(h) isa typeof(G)
   @test parent(g) == G
   @test parent(h) == G

   @testset "Parent methods" begin
      @test elem_type(typeof(G)) == typeof(g)
      @test elem_type(G) == typeof(g)
      @test parent_type(typeof(g)) == typeof(G)
      @test parent_type(g) == typeof(G)
      @test one(G) isa typeof(g)
      @test one(G)==one(g)==one(h)

      @test isfinite(G) isa Bool
#      @test hasorder(G) isa Bool
#      @test hasgens(G) isa Bool
      @test ngens(G) isa Int
      @test gens(G) isa Vector{typeof(g)}

      if isfinite(G)
         @test order(G) isa fmpz
         @test order(G) > 0
         @test is_trivial(G) == (order(G) == 1)
      else
        @test_throws GroupsCore.InfiniteOrder{typeof(G)} order(G)
      end
   end

   @testset "Comparison methods" begin
      if typeof(G) isa PermGroup
      @test (g==h) isa Bool
      @test isequal(g,h) isa Bool
      @test g == g
      @test isequal(h,h)
      @test (g<h) isa Bool
      @test isless(g,h) isa Bool
      @test g>h || g==h || g<h
      @test isequal(g,h) || isless(g,h) || isless(h,g)
      end
   end

   @testset "Group operations" begin
#      g1,h1 = deepcopy(g), deepcopy(h)
      g1,h1 = g,h

      @test inv(g) isa typeof(g)
      @test (g,h) == (g1,h1)
      @test g*h isa typeof(g)
      @test (g,h) == (g1,h1)
      @test g^2 == g*g
      @test (g,h) == (g1,h1)
      @test g^-3 == inv(g)*inv(g)*inv(g)
      @test (g,h) == (g1,h1)
      @test (g*h)^-1 == inv(h)*inv(g)
      @test (g,h) == (g1,h1)
      @test conj(g,h) == inv(h)*g*h
      @test (g,h) == (g1,h1)
      @test comm(g,h) == g^-1*h^-1*g*h
      @test (g,h) == (g1,h1)
      @test isone(g*inv(g)) && isone(inv(g)*g)
   end

   @testset "In-place operations" begin
#      g1,h1 = deepcopy(g), deepcopy(h)
      g1,h1 = g,h
      out = rand(G)  #to be replaced by out=similar(g)

      @test isone(one!(g))
#      g = deepcopy(g1)
      g = g1

      @testset "mul!" begin
         @test mul!(out,g,h) == g1*h1
         @test (g,h) == (g1,h1)

         @test mul!(out,g,h) == g1*h1
         @test (g,h) == (g1,h1)

         @test mul!(g,g,h) == g1*h1
         @test h==h1
#         g = deepcopy(g1)
         g = g1

         @test mul!(h,g,h) == g1*h1
         @test g == g1
#         h = deepcopy(h1)
         h = h1

         @test mul!(g,g,g) == g1*g1
#         g = deepcopy(g1)
         g = g1
      end

      @testset "conj!" begin
         res = h1^-1*g1*h1
         @test conj!(out,g,h) == res
         @test (g,h) == (g1,h1)

         @test conj!(g,g,h) == res
         @test h == h1
#         g = deepcopy(g1)
         g = g1

         @test conj!(h,g,h) == res
         @test g == g1
#         h = deepcopy(h1)
         h = h1

         @test conj!(g,g,g) == g1
#         g = deepcopy(g1)
         g = g1
      end

      @testset "comm!" begin
         res = g1^-1*h1^-1*g*h

         @test comm!(out,g,h) == res
         @test (g,h) == (g1,h1)

         @test comm!(g,g,h) == res
         @test h == h1
#         g = deepcopy(g1)
         g = g1

         @test comm!(h,g,h) == res
         @test g == g1
#         h = deepcopy(h1)
         h = h1
      end

      @testset "div_[left|right]!" begin
         res = g*h^-1

         @test div_right!(out,g,h) == res
         @test (g,h) == (g1,h1)

         @test div_right!(g,g,h) == res
         @test h == h1
#         g = deepcopy(g1)
         g = g1

         @test div_right!(h,g,h) == res
         @test g == g1
#         h = deepcopy(h1)
         h = h1

         @test div_right!(g,g,g) == one(g)
#         g = deepcopy(g1)
         g = g1

         res = h^-1*g

         @test div_left!(out,g,h) == res
         @test (g,h) == (g1,h1)

         @test div_left!(g,g,h) == res
         @test h == h1
#         g = deepcopy(g1)
         g = g1

         @test div_left!(h,g,h) == res
         @test g == g1
#         h = deepcopy(h1)
         h = h1

         @test div_left!(g,g,g) == one(g)
#         g = deepcopy(g1)
         g = g1
      end
   end
end


@testset "Iteration" begin
  for n = 4:6
    G = symmetric_group(n)
    L = [x for x in G]
    @test L isa Vector{PermGroupElem}
    @test length(L) == factorial(degree(G))
    @test length(unique(L)) == factorial(degree(G))
    @test rand(G) isa PermGroupElem
    @test rand(G) in G
    A = PermGroupElem[]
    for x in G
      push!(A, x)
    end
    @test length(A) == factorial(degree(G))
    s = 0         # check if the number of (n-1)-cycles is correct
    for x in G 
      if order(x) == (n-1)
        s+=1
      end
    end
    @test s == factorial(n-2)*n
  end
end


