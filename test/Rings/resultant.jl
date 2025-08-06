@testset "resultant" begin
  QQt, t = QQ[:t]
  R, (x1, x2, x0) = QQt[:x1, :x2, :x0]
  F = [x2^2 - t*x0*x1 , x2^2-t*x0*x2, x1^2+x2^2-x0^2]
  # this example requires coordinate change
  @test resultant(F) == -2*t^6 + t^4

  @test resultant([x1, x1, x1]) == 0

  R, (x, y, z) = QQ[:x, :y, :z]
  F0 = x^3 + y^2 * z
  F1 = x*y + y^2 + x*z + y * z
  F2 = y^4 + z^4
  r = resultant([F0, F1, F2])
  @test parent(r) === QQ && r == 16

  @test_throws ArgumentError resultant(typeof(F0)[])
  @test_throws ArgumentError resultant([F0])
  @test_throws ArgumentError resultant([F0, F0, F0, F0])
  @test_throws ArgumentError resultant([x, y^2, z^3 + x])
  _R, (x1, x2) = residue_ring(ZZ, 6)[1]["x1", "x2"]
  @test_throws ArgumentError resultant([x1, x2])

  # example from the Macaulay2 documentation:
   QQab, (a, b) = QQ[:a, :b]
   R, (x, y, z, w) = QQab[:x, :y, :z, :w]
   F = [(7//3)*x+(7//2)*y+z+2*w,
        ((10//7)*a+b)*x^2+(a+(5//4)*b)*x*y+(2*a+(1//2)*b)*y^2+((7//8)*a+(7//5)*b)*x*z+((3//4)*a+b)*y*z+((7//8)*a+(1//7)*b)*z^2+((5//7)*a+(4//3)*b)*x*w+(9*a+10*b)*y*w+((7//5)*a+(3//4)*b)*z*w+((4//3)*a+5*b)*w^2,
        ((1//2)*a+(7//5)*b)*x^3+((1//2)*a+10*b)*x^2*y+((8//9)*a+(3//5)*b)*x*y^2+(a+(7//6)*b)*y^3+((3//7)*a+(3//4)*b)*x^2*z+((1//3)*a+(9//10)*b)*x*y*z+((9//4)*a+b)*y^2*z+((1//6)*a+(1//5)*b)*x*z^2+(3*a+(5//2)*b)*y*z^2+((5//3)*a+(3//7)*b)*z^3+(a+b)*x^2*w+((4//5)*a+(5//4)*b)*x*y*w+((5//3)*a+(5//8)*b)*y^2*w+((3//2)*a+(1//6)*b)*x*z*w+((1//3)*a+(4//5)*b)*y*z*w+(9*a+(1//3)*b)*z^2*w+((7//3)*a+(5//4)*b)*x*w^2+(a+(3//4)*b)*y*w^2+((9//8)*a+(7//8)*b)*z*w^2+((9//7)*a+2*b)*w^3,
        2*x+(1//4)*y+(8//3)*z+(4//5)*w]
  @test resultant(F) == -(21002161660529014459938925799//2222549728809984000000)*a^5-(2085933800619238998825958079203//12700284164628480000000)*a^4*b-(348237304382147063838108483692249//889019891523993600000000)*a^3*b^2-(38379949248928909714532254698073//35278567123968000000000)*a^2*b^3-(1146977327343523453866040839029//1119954511872000000000)*a*b^4-(194441910898734675845094443//895963609497600000)*b^5
end
