
@testset "Operations on PermGroup" begin
   @testset "Elements of Sym($i)" for i in 8:16
      if i>1
      G=symmetric_group(i)
      x=@inferred rand(G)
      y=rand(G)
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
      @test y==G(listperm(y))
      end
   end
end
