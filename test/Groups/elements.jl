Oscar.@_AuxDocTest "show and print fp group elements", (fix = false),
raw"""
common setup

```jldoctest FPGroupElem.show
julia> using Oscar

julia> F = free_group('a':'z');

julia> long_word = prod(gens(F))*prod(inv.(gens(F)));
```

default `show` without unicode

```jldoctest FPGroupElem.show
julia> old = allow_unicode(false; temporary=true);

julia> withenv("LINES" => 30, "COLUMNS" => 80) do
         show(stdout, long_word)
       end
a*b*c*d*e*f*g*h*i*j*k*l*m*n*o*p*q*r*s*t*u*v*w*x*y*z*a^-1*b^-1*c^-1*d^-1*e^-1*f^-1*g^-1*h^-1*i^-1*j^-1*k^-1*l^-1*m^-1*n^-1*o^-1*p^-1*q^-1*r^-1*s^-1*t^-1*u^-1*v^-1*w^-1*x^-1*y^-1*z^-1

julia> withenv("LINES" => 30, "COLUMNS" => 80) do
         show(stdout, MIME("text/plain"), long_word)
       end
a*b*c*d*e*f*g*h*i*j*k*l*m*n*o*p*q*r*s*t*u*v*w*x*y*z*a^-1*b^-1*c^-1*d^-1*e^-1*f^-1*g^-1*h^-1*i^-1*j^-1*k^-1*l^-1*m^-1*n^-1*o^-1*p^-1*q^-1*r^-1*s^-1*t^-1*u^-1*v^-1*w^-1*x^-1*y^-1*z^-1

julia> withenv("LINES" => 30, "COLUMNS" => 80) do
         show(IOContext(stdout, :limit => true), long_word)
       end
a*b*c*d*e*f*g*h*i*j*k*l*m*n*o*p*q*r*s*t*u*v*w*x*y*z*a^-1*b^-1*c^-1*d^-1*e^...

julia> withenv("LINES" => 30, "COLUMNS" => 80) do
         show(IOContext(stdout, :limit => true), MIME("text/plain"), long_word)
       end
a*b*c*d*e*f*g*h*i*j*k*l*m*n*o*p*q*r*s*t*u*v*w*x*y*z*a^-1*b^-1*c^-1*d^-1*e^...

julia> allow_unicode(old; temporary=true);
```

default `show` with unicode

```jldoctest FPGroupElem.show
julia> old = allow_unicode(true; temporary=true);

julia> withenv("LINES" => 30, "COLUMNS" => 80) do
         show(stdout, long_word)
       end
a*b*c*d*e*f*g*h*i*j*k*l*m*n*o*p*q*r*s*t*u*v*w*x*y*z*a^-1*b^-1*c^-1*d^-1*e^-1*f^-1*g^-1*h^-1*i^-1*j^-1*k^-1*l^-1*m^-1*n^-1*o^-1*p^-1*q^-1*r^-1*s^-1*t^-1*u^-1*v^-1*w^-1*x^-1*y^-1*z^-1

julia> withenv("LINES" => 30, "COLUMNS" => 80) do
         show(stdout, MIME("text/plain"), long_word)
       end
a*b*c*d*e*f*g*h*i*j*k*l*m*n*o*p*q*r*s*t*u*v*w*x*y*z*a^-1*b^-1*c^-1*d^-1*e^-1*f^-1*g^-1*h^-1*i^-1*j^-1*k^-1*l^-1*m^-1*n^-1*o^-1*p^-1*q^-1*r^-1*s^-1*t^-1*u^-1*v^-1*w^-1*x^-1*y^-1*z^-1

julia> withenv("LINES" => 30, "COLUMNS" => 80) do
         show(IOContext(stdout, :limit => true), long_word)
       end
a*b*c*d*e*f*g*h*i*j*k*l*m*n*o*p*q*r*s*t*u*v*w*x*y*z*a^-1*b^-1*c^-1*d^-1*e^-1…

julia> withenv("LINES" => 30, "COLUMNS" => 80) do
         show(IOContext(stdout, :limit => true), MIME("text/plain"), long_word)
       end
a*b*c*d*e*f*g*h*i*j*k*l*m*n*o*p*q*r*s*t*u*v*w*x*y*z*a^-1*b^-1*c^-1*d^-1*e^-1…

julia> allow_unicode(old; temporary=true);
```
"""

@testset "Permutations" begin
  for n = 10:13
    G=symmetric_group(n)
    x=cperm(1:5,6:8,9:n)
    A=vcat([i for i in 10:n],[9])
    A=vcat([2,3,4,5,1,7,8,6],A)
    y=G(A)

    @test x==y
    @test A==Vector(y)
    @test typeof(Vector(y))==Vector{Int64}
    @test typeof(Vector{ZZRingElem}(y))==Vector{ZZRingElem}
    @test x==G(Vector(x))
    @test is_finite_order(x)
    @test order(x) == lcm(15,n-8)
    for T in [Int, BigInt, ZZRingElem]
      @test order(T, x) == lcm(15,n-8)
      @test order(T, x) isa T
    end

    @test x(3)==4
    @test x(8)==6
    @test x(n)==9

    @test 3^x == 4
    @test n^x == 9

    for T in [Int32, Int, BigInt, ZZRingElem]
      @test x(T(n)) == T(9)
      @test typeof(x(T(n))) == T
      @test T(n)^x == T(9)
      @test typeof(T(n)^x) == T
    end
  end

  G=symmetric_group(6)
  x=perm([2,3,4,5,6,1])
  @test x==perm(Int8[2,3,4,5,6,1])
  @test x==perm(ZZRingElem[2,3,4,5,6,1])
  @test x==perm(G,[2,3,4,5,6,1])
  @test x==perm(G,Int8[2,3,4,5,6,1])
  @test x==perm(G,ZZRingElem[2,3,4,5,6,1])
  @test cperm(G,Int[])==one(G)
  @test x==cperm(G,1:6)
  @test x==cperm(G,[1,2,3,4,5,6])
  @test one(G)==perm(1:6)
  @test_throws ArgumentError G(perm([2,3,4,5,6,7,1]))
  @test_throws ArgumentError G([2,3,1,4,6,5,7])
  @test G(perm([2,3,1,4,6,5,7]))==perm([2,3,1,4,6,5])
  @test_throws ArgumentError perm(G,[2,3,4,5,6,7,1])
  @test_throws ArgumentError perm(G, [1,1])
  @test one(G)==cperm(G,Int64[])

  G=alternating_group(6)
  @test_throws ArgumentError cperm(G, [1,2])
  @test cperm(G, [1,2],[3,4]) == perm(G,[2,1,4,3,5,6])
  @test cperm(G, [1,2],[2,3]) == cperm(G, [1,3,2])
  @test cperm(G, [1,2],[2,3]) == perm(G,[3,1,2,4,5,6])

end

@testset "Change of parent" begin
  for n = 10:12
    H = symmetric_group(n+3)
    K = symmetric_group(n-1)
    x = cperm(H, 1:5,6:8,9:n)
    @test parent(x)===H
    @test_throws ArgumentError cperm(K,1:5,6:8,9:n)
    y = H(x)
    @test parent(y) == H
    @test parent(y) === H
    @test parent(x) === H
    @test_throws ArgumentError K(x)

#    z = perm(vcat(2:(n-2),[1]))
#    @test parent(z) == symmetric_group(n-2)
    z = perm(symmetric_group(n),vcat(2:(n-2),[1]))
    @test parent(z) == symmetric_group(n)
    @test_throws MethodError perm(symmetric_group(n-3),z)

    @test cperm(K,Int64[]) == one(K)
  end

  for G in [
            PermGroup(small_group(2, 1))
            PcGroup(small_group(2, 1))
            FPGroup(small_group(2, 1))
            GL(2,2)
           ]
    H = rand(rand(subgroup_classes(G)))
    @test parent(one(H)) === H
    @test parent(G(one(H))) === G
  end


end

@testset "Eltypes" begin
   @test eltype(PermGroup)==PermGroupElem
   @test eltype(PcGroup)==PcGroupElem
   @test eltype(FPGroup)==FPGroupElem
   @test eltype(GL(2,3))==MatrixGroupElem{elem_type(typeof(GF(2))),dense_matrix_type(GF(2))}
   @test eltype(DirectProductGroup)==Oscar.BasicGAPGroupElem{DirectProductGroup}
   @test eltype(direct_product(symmetric_group(3),cyclic_group(2)))==Oscar.BasicGAPGroupElem{DirectProductGroup}
   @test eltype(SemidirectProductGroup)==Oscar.BasicGAPGroupElem{SemidirectProductGroup}
   @test eltype(WreathProductGroup)==Oscar.BasicGAPGroupElem{WreathProductGroup}
   @test eltype(AutomorphismGroup{PcGroup})==Oscar.BasicGAPGroupElem{AutomorphismGroup{PcGroup}}

   G = symmetric_group(5)
   x = cperm([1,4,2,5])
   H = sub(G,[x])[1]
   y = cperm(G,[2,3,4])
   w = cperm(G,[1,4])
   K = sub(G,[w])[1]
   cc = conjugacy_class(G,y)
   cs = conjugacy_class(G,H)
   lc = x*H
   rc = H*x
   dc = K*x*H
   @test [z for z in G] == @inferred collect(G)
   @test [z for z in cc] == @inferred collect(cc)
   @test [z for z in cs] == @inferred collect(cs)
   @test [z for z in lc] == @inferred collect(lc)
   @test [z for z in rc] == @inferred collect(rc)
   @test [z for z in dc] == @inferred collect(dc)
   @test typeof(collect(G))==Vector{typeof(x)}
   @test typeof(collect(lc))==Vector{typeof(x)}
   @test typeof(collect(rc))==Vector{typeof(x)}
   @test typeof(collect(dc))==Vector{typeof(x)}
   @test typeof(collect(cc))==Vector{typeof(x)}
   @test typeof(collect(cs))==Vector{typeof(H)}
   @test eltype(cc)==typeof(y)
   @test eltype(cs)==typeof(H)
   @test eltype(lc)==typeof(y)
   @test eltype(rc)==typeof(y)
   @test eltype(dc)==typeof(y)
end

@testset "Generators" begin
   L=[symmetric_group(4), cyclic_group(5), free_group(3), symplectic_group(4,3)]
   for G in L
      K=gens(G)
      @test length(K) == ngens(G)
      for i in 1:length(K)
         @test K[i] == G[i]
         @test K[i] == gen(G,i)
      end
      @test G[0] == gen(G, 0)
      @test G[0] == one(G)
      @test_throws BoundsError K[0]
   end

   G = free_group(2)
   @test_throws ArgumentError gen(G, 3)
   @test_throws ArgumentError gen(G, -3)
end

@testset "deepcopy" begin
   for g in [symmetric_group(5), free_group(2), small_group(8, 1),
             automorphism_group(alternating_group(4))]
     m = Oscar.BasicGAPGroupElem(g, GAP.Obj(gen(g, 1)))
     @test isdefined(m, :X)
     c = deepcopy(m);
     @test isdefined(c, :X)
     @test GAP.Obj(c) == GAP.Obj(m)

     @test deepcopy([one(g)]) == [one(g)]
   end
end

@testset "compatibility of parents" begin
   G = symmetric_group(4)
   g = symmetric_group(3)
   a = automorphism_group(g)
   L = [G,
        automorphism_group(alternating_group(4)),
        general_linear_group(2, 3),
        direct_product(cyclic_group(2), cyclic_group(3)),
        semidirect_product(g, id_hom(a), a),
        wreath_product(symmetric_group(2), symmetric_group(3))]
   for T in [FPGroup, SubFPGroup, PcGroup, SubPcGroup]
     push!(L, codomain(isomorphism(T, G)))
   end
   @testset for g in L
     s2 = sylow_subgroup(g, 2)[1]
     s3 = sylow_subgroup(g, 3)[1]
     @test parent(one(s2) * one(g)) == g
     @test parent(one(g) * one(s2)) == g
     @test parent(one(s2) * one(s3)) == g
     @test parent(one(s2) * one(s2)) == s2

     x = gen(s2, 1)
     y = gen(s3, 1)
     z = gen(g, 1)
     c = conj(x, y)
     @test c == inv(y) * x * y
     @test c == x^y
     @test parent(c) == parent(inv(y) * x * y)
     c = conj(x, z)
     @test c == inv(z) * x * z
     @test c == x^z
     @test parent(c) == parent(inv(z) * x * z)
   end

   @test_throws MethodError one(L[1]) * one(L[2])  # no generic method for different types
   @test_throws ArgumentError one(small_group(12, 1)) * one(small_group(12,2))
   g1 = codomain(isomorphism(FPGroup, L[1]))
   g2 = codomain(isomorphism(FPGroup, L[2]))
   @test_throws ArgumentError one(g1) * one(g2)

   g = free_group(2)
   x, y = gens(g)
   f, epi = quo(g, [x, y^2])
   @test free_group(g) === g
   @test free_group(f) == g
   @test underlying_word(x) === x
   @test underlying_word(epi(x)) == x
end
