@testset "Rings.ZZ.constructors" begin
   @test ZZ(123) isa Oscar.fmpz
   @test ZZ(2348732479832498732749823) isa Oscar.fmpz
   @test ZZ(49874359874359874359874387598374587984375897435) isa Oscar.fmpz
   @test ZZ() isa Oscar.fmpz
   @test ZZ(0) isa Oscar.fmpz
   @test ZZ(-123) isa Oscar.fmpz

   @test zero(ZZ) isa Oscar.fmpz
   @test zero(ZZ) == 0
   @test one(ZZ) isa Oscar.fmpz
   @test one(ZZ) == 1

   @test ZZ(typemin(Int)) == typemin(Int)
   @test ZZ(typemin(Int)) < 0
   @test ZZ(typemax(Int)) == typemax(Int)
   @test ZZ(typemax(Int)) > 0
end

@testset "Rings.ZZ.properties" begin
   @test iszero(zero(ZZ))
   @test isone(one(ZZ))
   @test is_unit(one(ZZ))
   @test is_unit(-one(ZZ))
   @test !is_unit(zero(ZZ))
   @test sign(ZZ(2)) == 1
   @test sign(ZZ(0)) == 0
   @test sign(-ZZ(2)) == -1
   @test sign(ZZ(2)) isa Oscar.fmpz
   @test sign(ZZ(0)) isa Oscar.fmpz
   @test sign(ZZ(-2)) isa Oscar.fmpz
   @test abs(ZZ(2)) == 2
   @test abs(ZZ(-2)) == 2
   @test abs(ZZ(0)) == 0
   @test isodd(ZZ(3))
   @test isodd(ZZ(-1))
   @test !isodd(ZZ(2))
   @test !isodd(ZZ(0))
   @test !isodd(ZZ(-2))
   @test iseven(ZZ(2))
   @test iseven(ZZ(0))
   @test iseven(ZZ(-2))
   @test !iseven(ZZ(-1))
   @test is_square(ZZ(0))
   @test is_square(ZZ(1))
   @test !is_square(ZZ(2))
   @test is_square(ZZ(4))
   @test !is_square(ZZ(-1))
   @test !is_square(ZZ(-4))
   @test is_prime(ZZ(2))
   @test !is_prime(ZZ(-2))
   @test !is_prime(ZZ(-1))
   @test !is_prime(ZZ(0))
   @test !is_prime(ZZ(1))
   @test is_probable_prime(ZZ(2))
   @test !is_probable_prime(ZZ(-2))
   @test !is_probable_prime(ZZ(-1))
   @test !is_probable_prime(ZZ(0))
   @test !is_probable_prime(ZZ(1))

   @test parent(ZZ(2)) == ZZ
   @test parent(ZZ()) == ZZ
end

@testset "Rings.ZZ.arithmetic" begin
   @test ZZ(2) + 3 isa Oscar.fmpz
   @test 3 + ZZ(2) isa Oscar.fmpz
   @test ZZ(2) - 3 isa Oscar.fmpz
   @test 3 - ZZ(2) isa Oscar.fmpz
   @test ZZ(2)*3 isa Oscar.fmpz
   @test 3*ZZ(2) isa Oscar.fmpz
   @test 0*ZZ(2) isa Oscar.fmpz
   @test ZZ(2)*0 isa Oscar.fmpz

   a = ZZ(2)

   @test a !== a + 0
   @test a !== a + BigInt(0)
   @test a !== a + ZZ(0)
   @test a !== 0 + a
   @test a !== BigInt(0) + a
   @test a !== ZZ(0) + a
   
   @test a !== a - 0
   @test a !== a - BigInt(0)
   @test a !== a - ZZ(0)

   @test a !== a*1
   @test a !== a*BigInt(1)
   @test a !== a*ZZ(1)
   @test a !== 1*a
   @test a !== BigInt(1)*a
   @test a !== ZZ(1)*a
end

@testset "Rings.ZZ.comparison" begin
   @test ZZ(2) == ZZ(2)
   @test ZZ(0) == ZZ(0)
   @test ZZ(-2) == ZZ(-2)
   @test ZZ(2) == 2
   @test ZZ(0) == 0
   @test ZZ(-2) == -2
   @test 2 == ZZ(2)
   @test 0 == ZZ(0)
   @test -2 == ZZ(-2)

   @test ZZ(3) > ZZ(2)
   @test ZZ(2) >= ZZ(2)
   @test ZZ(2) > ZZ(0)
   @test ZZ(2) < ZZ(3)
   @test ZZ(2) <= ZZ(2)
   @test ZZ(0) < ZZ(2)
   @test ZZ(2) > ZZ(-2)
   @test ZZ(-2) < ZZ(2)

   @test ZZ(3) > 2
   @test ZZ(2) >= 2
   @test ZZ(2) > 0
   @test ZZ(2) < 3
   @test ZZ(2) <= 2
   @test ZZ(0) < 2
   @test ZZ(2) > -2
   @test ZZ(-2) < 2

   @test 3 > ZZ(2)
   @test 2 >= ZZ(2)
   @test 2 > ZZ(0)
   @test 2 < ZZ(3)
   @test 2 <= ZZ(2)
   @test 0 < ZZ(2)
   @test 2 > ZZ(-2)
   @test -2 < ZZ(2)

   @test ZZ(2) < 32783627361723381726834
   @test ZZ(2) < 4243746873264873264873264876327486832468732
   @test ZZ(2) > -32783627361723381726834
   @test ZZ(2) > -4243746873264873264873264876327486832468732
   @test 32783627361723381726834 > ZZ(2)
   @test 4243746873264873264873264876327486832468732 > ZZ(2)
   @test -32783627361723381726834 < ZZ(2)
   @test -4243746873264873264873264876327486832468732 < ZZ(2)

   @test ZZ(32783627361723381726834) == 32783627361723381726834
   @test ZZ(4243746873264873264873264876327486832468732) == 4243746873264873264873264876327486832468732

   @test ZZ(2) != ZZ(3)
   @test ZZ(2) != 3
   @test 2 != ZZ(3)
   @test ZZ(2) != 32783627361723381726834
   @test ZZ(2) != 4243746873264873264873264876327486832468732
   @test 32783627361723381726834 != ZZ(2)
   @test 4243746873264873264873264876327486832468732 != ZZ(2)
end

@testset "Rings.ZZ.divexact" begin
   @test divexact(ZZ(6), ZZ(2)) == 3
   @test divexact(ZZ(6), ZZ(-2)) == -3
   @test divexact(ZZ(-6), ZZ(-2)) == 3
   @test divexact(ZZ(-6), ZZ(2)) == -3

   @test_throws DivideError divexact(ZZ(2), ZZ(0))
   @test_throws DivideError divexact(ZZ(0), ZZ(0))
   @test_throws DivideError divexact(ZZ(2), 0)
   @test_throws DivideError divexact(ZZ(0), 0)

   @test_throws ArgumentError divexact(ZZ(2), ZZ(3))
   @test_throws ArgumentError divexact(ZZ(2), 3)
   @test_throws ArgumentError divexact(2, ZZ(3))

   @test divexact(ZZ(6), 2) isa Oscar.fmpz
   @test divexact(6, ZZ(2)) isa Oscar.fmpz
   @test divexact(ZZ(32783627361723381726834), 32783627361723381726834) isa Oscar.fmpz
   @test divexact(32783627361723381726834, ZZ(32783627361723381726834)) isa Oscar.fmpz
   @test divexact(ZZ(4243746873264873264873264876327486832468732), 4243746873264873264873264876327486832468732) isa Oscar.fmpz
   @test divexact(32783627361723381726834, ZZ(32783627361723381726834)) isa Oscar.fmpz

   @test divexact(ZZ(6), 2) == 3
   @test divexact(6, ZZ(2)) == 3
   @test divexact(ZZ(32783627361723381726834), 32783627361723381726834) == 1
   @test divexact(32783627361723381726834, ZZ(32783627361723381726834)) == 1
   @test divexact(4243746873264873264873264876327486832468732, ZZ(4243746873264873264873264876327486832468732)) == 1

   @test divexact(ZZ(0), ZZ(2)) == 0
   @test divexact(0, ZZ(2)) isa Oscar.fmpz
   @test divexact(0, ZZ(2)) == 0
end

@testset "Rings.ZZ.powering" begin
   @test ZZ(2)^2 isa Oscar.fmpz
   @test ZZ(2)^1 isa Oscar.fmpz
   @test ZZ(2)^0 isa Oscar.fmpz
   @test ZZ(1)^-1 isa Oscar.fmpz
   @test ZZ(1)^-2 isa Oscar.fmpz
   @test ZZ(0)^2 isa Oscar.fmpz
   @test ZZ(0)^1 isa Oscar.fmpz
   @test ZZ(0)^0 isa Oscar.fmpz

   @test ZZ(2)^2 == 4
   @test ZZ(2)^1 == 2
   @test ZZ(2)^0 == 1
   @test ZZ(1)^-1 == 1
   @test ZZ(1)^-2 == 1
   @test ZZ(-1)^-1 == -1
   @test ZZ(-1)^-2 == 1
   @test ZZ(0)^2 == 0
   @test ZZ(0)^1 == 0
   
   @test ZZ(0)^0 == 1 # for convenience only

   @test_throws DomainError ZZ(0)^-1
   @test_throws DomainError ZZ(0)^-2
   @test_throws DomainError ZZ(-2)^-1
   @test_throws DomainError ZZ(2)^-1
   @test_throws DomainError ZZ(-2)^-2
   @test_throws DomainError ZZ(2)^-1

   @test ZZ(-3)^3 == -27
   @test ZZ(-1)^3 == -1
   @test ZZ(-1)^2 == 1
   @test ZZ(-1)^1 == -1

   @test ZZ(2)^62 isa Oscar.fmpz
   @test ZZ(2)^63 isa Oscar.fmpz
   @test ZZ(2)^64 isa Oscar.fmpz
   @test ZZ(2)^62 == BigInt(2)^62
   @test ZZ(2)^63 == BigInt(2)^63
   @test ZZ(2)^64 == BigInt(2)^64
end

@testset "Rings.ZZ.euclidean_division" begin
   @test mod(ZZ(2), ZZ(3)) == mod(2, 3)
   @test mod(ZZ(2), ZZ(-3)) == mod(2, -3)
   @test mod(ZZ(-2), ZZ(3)) == mod(-2, 3)
   @test mod(ZZ(-2), ZZ(-3)) == mod(-2, -3)

   @test rem(ZZ(2), ZZ(3)) == rem(2, 3)
   @test rem(ZZ(2), ZZ(-3)) == rem(2, -3)
   @test rem(ZZ(-2), ZZ(3)) == rem(-2, 3)
   @test rem(ZZ(-2), ZZ(-3)) == rem(-2, -3)

   @test div(ZZ(2), ZZ(3)) == div(2, 3)
   @test div(ZZ(2), ZZ(-3)) == div(2, -3)
   @test div(ZZ(-2), ZZ(3)) == div(-2, 3)
   @test div(ZZ(-2), ZZ(-3)) == div(-2, -3)

   @test divrem(ZZ(2), ZZ(3)) == divrem(2, 3)
   @test divrem(ZZ(2), ZZ(-3)) == divrem(2, -3)
   @test divrem(ZZ(-2), ZZ(3)) == divrem(-2, 3)
   @test divrem(ZZ(-2), ZZ(-3)) == divrem(-2, -3)

   @test mod(ZZ(2), 3) == mod(2, 3)
   @test mod(ZZ(2), -3) == mod(2, -3)
   @test mod(ZZ(-2), 3) == mod(-2, 3)
   @test mod(ZZ(-2), -3) == mod(-2, -3)

   @test rem(ZZ(2), 3) == rem(2, 3)
   @test rem(ZZ(2), -3) == rem(2, -3)
   @test rem(ZZ(-2), 3) == rem(-2, 3)
   @test rem(ZZ(-2), -3) == rem(-2, -3)

   @test div(ZZ(2), 3) == div(2, 3)
   @test div(ZZ(2), -3) == div(2, -3)
   @test div(ZZ(-2), 3) == div(-2, 3)
   @test div(ZZ(-2), -3) == div(-2, -3)

   @test divrem(ZZ(2), 3) == divrem(2, 3)
   @test divrem(ZZ(2), -3) == divrem(2, -3)
   @test divrem(ZZ(-2), 3) == divrem(-2, 3)
   @test divrem(ZZ(-2), -3) == divrem(-2, -3)

   @test mod(2, ZZ(3)) == mod(2, 3)
   @test mod(2, ZZ(-3)) == mod(2, -3)
   @test mod(-2, ZZ(3)) == mod(-2, 3)
   @test mod(-2, ZZ(-3)) == mod(-2, -3)

   @test rem(2, ZZ(3)) == rem(2, 3)
   @test rem(2, ZZ(-3)) == rem(2, -3)
   @test rem(-2, ZZ(3)) == rem(-2, 3)
   @test rem(-2, ZZ(-3)) == rem(-2, -3)

   @test div(2, ZZ(3)) == div(2, 3)
   @test div(2, ZZ(-3)) == div(2, -3)
   @test div(-2, ZZ(3)) == div(-2, 3)
   @test div(-2, ZZ(-3)) == div(-2, -3)

   @test divrem(2, ZZ(3)) == divrem(2, 3)
   @test divrem(2, ZZ(-3)) == divrem(2, -3)
   @test divrem(-2, ZZ(3)) == divrem(-2, 3)
   @test divrem(-2, ZZ(-3)) == divrem(-2, -3)

   @test mod(2, ZZ(3)) isa Oscar.fmpz
   @test rem(2, ZZ(3)) isa Oscar.fmpz
   @test div(2, ZZ(3)) isa Oscar.fmpz
   @test divrem(2, ZZ(3)) isa Tuple{Oscar.fmpz, Oscar.fmpz}

   @test mod(ZZ(2), 3) isa Oscar.fmpz
   @test rem(ZZ(2), 3) isa Oscar.fmpz
   @test div(ZZ(2), 3) isa Oscar.fmpz
   @test divrem(ZZ(2), 3) isa Tuple{Oscar.fmpz, Oscar.fmpz}

   @test mod(Int128(2), ZZ(3)) isa Oscar.fmpz
   @test rem(Int128(2), ZZ(3)) isa Oscar.fmpz
   @test div(Int128(2), ZZ(3)) isa Oscar.fmpz
   @test divrem(Int128(2), ZZ(3)) isa Tuple{Oscar.fmpz, Oscar.fmpz}

   @test mod(ZZ(2), Int128(3)) isa Oscar.fmpz
   @test rem(ZZ(2), Int128(3)) isa Oscar.fmpz
   @test div(ZZ(2), Int128(3)) isa Oscar.fmpz
   @test divrem(ZZ(2), Int128(3)) isa Tuple{Oscar.fmpz, Oscar.fmpz}

   @test mod(BigInt(2), ZZ(3)) isa Oscar.fmpz
   @test rem(BigInt(2), ZZ(3)) isa Oscar.fmpz
   @test div(BigInt(2), ZZ(3)) isa Oscar.fmpz
   @test divrem(BigInt(2), ZZ(3)) isa Tuple{Oscar.fmpz, Oscar.fmpz}

   @test mod(ZZ(2), BigInt(3)) isa Oscar.fmpz
   @test rem(ZZ(2), BigInt(3)) isa Oscar.fmpz
   @test div(ZZ(2), BigInt(3)) isa Oscar.fmpz
   @test divrem(ZZ(2), BigInt(3)) isa Tuple{Oscar.fmpz, Oscar.fmpz}

   @test_throws DivideError mod(ZZ(2), ZZ(0))
   @test_throws DivideError mod(ZZ(0), ZZ(0))
   @test_throws DivideError mod(ZZ(2), 0)
   @test_throws DivideError mod(ZZ(0), 0)
   @test_throws DivideError mod(2, ZZ(0))
   @test_throws DivideError mod(0, ZZ(0))

   @test_throws DivideError rem(ZZ(2), ZZ(0))
   @test_throws DivideError rem(ZZ(0), ZZ(0))
   @test_throws DivideError rem(ZZ(2), 0)
   @test_throws DivideError rem(ZZ(0), 0)
   @test_throws DivideError rem(2, ZZ(0))
   @test_throws DivideError rem(0, ZZ(0))

   @test_throws DivideError div(ZZ(2), ZZ(0))
   @test_throws DivideError div(ZZ(0), ZZ(0))
   @test_throws DivideError div(ZZ(2), 0)
   @test_throws DivideError div(ZZ(0), 0)
   @test_throws DivideError div(2, ZZ(0))
   @test_throws DivideError div(0, ZZ(0))

   @test_throws DivideError divrem(ZZ(2), ZZ(0))
   @test_throws DivideError divrem(ZZ(0), ZZ(0))
   @test_throws DivideError divrem(ZZ(2), 0)
   @test_throws DivideError divrem(ZZ(0), 0)
   @test_throws DivideError divrem(2, ZZ(0))
   @test_throws DivideError divrem(0, ZZ(0))
end

@testset "Rings.ZZ.conversion" begin
   @test Int(ZZ(123)) == 123
   @test Int(ZZ(0)) == 0
   @test Int(ZZ(-123)) == -123
   @test Int(ZZ(typemin(Int))) == typemin(Int)
   @test Int(ZZ(typemax(Int))) == typemax(Int)

   @test_throws InexactError Int(ZZ(typemin(Int)) - 1)
   @test_throws InexactError Int(ZZ(typemax(Int)) + 1)

   @test fits(Int, ZZ(123))
   @test fits(Int, ZZ(typemin(Int)))
   @test fits(Int, ZZ(typemax(Int)))

   @test Int(ZZ(123)) isa Int
   @test Int(ZZ(0)) isa Int

   @test BigInt(ZZ(123)) == 123
   @test BigInt(ZZ(0)) == 0
   @test BigInt(ZZ(1234498543743857934594389)) == 1234498543743857934594389

   @test BigInt(ZZ(123)) isa BigInt
   @test BigInt(ZZ(0)) isa BigInt
end

@testset "Rings.ZZ.gcd" begin
   @test gcd(ZZ(2), ZZ(3)) == 1
   @test gcd(ZZ(2), ZZ(-3)) == 1
   @test gcd(ZZ(-2), ZZ(3)) == 1
   @test gcd(ZZ(-2), ZZ(-3)) == 1

   @test gcd(ZZ(2), ZZ(0)) == 2
   @test gcd(ZZ(0), ZZ(2)) == 2
   @test gcd(ZZ(0), ZZ(0)) == 0

   @test gcd(ZZ(2), 0) == 2
   @test gcd(0, ZZ(2)) == 2
   @test gcd(ZZ(0), 0) == 0
   @test gcd(0, ZZ(0)) == 0

   @test gcd(ZZ(2), 3) isa Oscar.fmpz
   @test gcd(2, ZZ(3)) isa Oscar.fmpz
   @test gcd(ZZ(2), 0) isa Oscar.fmpz
   @test gcd(0, ZZ(2)) isa Oscar.fmpz
   @test gcd(ZZ(0), 0) isa Oscar.fmpz
   @test gcd(0, ZZ(0)) isa Oscar.fmpz

   @test gcd(ZZ(2), Int128(3)) isa Oscar.fmpz
   @test gcd(ZZ(2), BigInt(3)) isa Oscar.fmpz
   @test gcd(Int128(2), ZZ(3)) isa Oscar.fmpz
   @test gcd(BigInt(2), ZZ(3)) isa Oscar.fmpz

   @test lcm(ZZ(2), ZZ(3)) == 6
   @test lcm(ZZ(2), ZZ(-3)) == 6
   @test lcm(ZZ(-2), ZZ(3)) == 6
   @test lcm(ZZ(-2), ZZ(-3)) == 6

   @test lcm(ZZ(2), ZZ(0)) == 0
   @test lcm(ZZ(0), ZZ(2)) == 0
   @test lcm(ZZ(0), ZZ(0)) == 0

   @test lcm(ZZ(2), 0) == 0
   @test lcm(0, ZZ(2)) == 0
   @test lcm(ZZ(0), 0) == 0
   @test lcm(0, ZZ(0)) == 0

   @test lcm(ZZ(2), 3) isa Oscar.fmpz
   @test lcm(2, ZZ(3)) isa Oscar.fmpz
   @test lcm(ZZ(2), 0) isa Oscar.fmpz
   @test lcm(0, ZZ(2)) isa Oscar.fmpz
   @test lcm(ZZ(0), 0) isa Oscar.fmpz
   @test lcm(0, ZZ(0)) isa Oscar.fmpz

   for a = -10:10
      for b = -10:10
         g, s, t = gcdx(ZZ(a), ZZ(b))
         gcdx(ZZ(a), b) == g, s, t
         gcdx(a, ZZ(b)) == g, s, t
         @test g == a*s + b*t
         if abs(a) == abs(b)
            @test t == sign(b)
         elseif b == 0 || abs(b) == 2*g
            @test s == sign(a)
         elseif a == 0 || abs(a) == 2*g
            @test t == sign(b)
         else
            @test 2*g*abs(s) < abs(b)
            @test 2*g*abs(t) < abs(a)
         end
      end
   end

   @test gcdx(ZZ(2), ZZ(3)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(2), ZZ(0)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(0), ZZ(3)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(0), ZZ(0)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(2), ZZ(1)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(1), ZZ(3)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(-1), ZZ(-1)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}

   @test gcdx(ZZ(2), 3) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(2), 0) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(0), 3) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(0), 0) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(2), 1) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(1), 3) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(ZZ(-1), -1) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}

   @test gcdx(2, ZZ(3)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(2, ZZ(0)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(0, ZZ(3)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(0, ZZ(0)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(2, ZZ(1)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(1, ZZ(3)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
   @test gcdx(-1, ZZ(-1)) isa Tuple{Oscar.fmpz, Oscar.fmpz, Oscar.fmpz}
end

@testset "Rings.ZZ.roots" begin
   @test isqrt(ZZ(16)) == 4
   @test isqrt(ZZ(5)) == 2
   @test isqrt(ZZ(0)) == 0
   @test isqrt(ZZ(1)) == 1

   @test_throws DomainError isqrt(ZZ(-1))
   @test_throws DomainError isqrt(ZZ(-4))

   @test isqrtrem(ZZ(16)) == (4, 0)
   @test isqrtrem(ZZ(5)) == (2, 1)
   @test isqrtrem(ZZ(0)) == (0, 0)
   @test isqrtrem(ZZ(1)) == (1, 0)

   @test_throws DomainError isqrtrem(ZZ(-1))
   @test_throws DomainError isqrtrem(ZZ(-4))

   @test root(ZZ(16), 2) == 4
   @test root(ZZ(5), 3) == 1
   @test root(ZZ(-8), 3) == -2
   @test root(ZZ(0), 2) == 0
   @test root(ZZ(0), 3) == 0
   @test root(ZZ(1), 1) == 1
   @test root(ZZ(-1), 1) == -1
   @test root(ZZ(1), 2) == 1

   @test_throws DomainError root(ZZ(2), -3)
   @test_throws DomainError root(ZZ(2), 0)
   @test_throws DomainError root(ZZ(1), -1)
   @test_throws DomainError root(ZZ(1), 0)
   @test_throws DomainError root(ZZ(0), 0)
   @test_throws DomainError root(ZZ(-2), 2)
   @test_throws DomainError root(ZZ(-2), 0)
end

@testset "Rings.ZZ.factorisation" begin
   F = factor(ZZ(-60))

   @test unit(F) == -1
   @test unit(F) isa Oscar.fmpz

   @test length(collect(F)) == 3
   @test 2 in F
   @test ZZ(3) in F
   @test F[2] == 2
   @test F[ZZ(3)] == 1
   @test !(ZZ(-1) in F)
   @test !(-1 in F)
   @test !(ZZ(0) in F)
   @test !(0 in F)
   @test !(7 in F)
   
   F = factor(ZZ(60))

   @test unit(F) == 1

   @test_throws ArgumentError factor(ZZ(0))
   @test_throws ErrorException F[0]
   @test_throws ErrorException F[7]
   @test_throws ErrorException F[ZZ(0)]
   @test_throws ErrorException F[ZZ(7)]
end

@testset "Rings.ZZ.combinatorial" begin
   @test factorial(ZZ(0)) == 1
   @test factorial(ZZ(1)) == 1
   @test factorial(ZZ(3)) == 6

   @test factorial(ZZ(0)) isa Oscar.fmpz
   @test factorial(ZZ(3)) isa Oscar.fmpz

   @test_throws DomainError factorial(ZZ(-1))

   @test rising_factorial(3, 2) == 12
   @test rising_factorial(3, 0) == 1
   @test rising_factorial(0, 3) == 0
   @test rising_factorial(-3, 3) == -6

   @test rising_factorial(2, 3) isa Int
   @test rising_factorial(2, 0) isa Int

   @test_throws DomainError rising_factorial(0, -1)
   @test_throws DomainError rising_factorial(-3, -5)

   @test rising_factorial(ZZ(3), 2) == 12
   @test rising_factorial(ZZ(3), 0) == 1
   @test rising_factorial(ZZ(0), 3) == 0
   @test rising_factorial(ZZ(-3), 3) == -6

   @test rising_factorial(ZZ(2), 3) isa Oscar.fmpz
   @test rising_factorial(ZZ(2), 0) isa Oscar.fmpz

   @test_throws DomainError rising_factorial(ZZ(0), -1)
   @test_throws DomainError rising_factorial(ZZ(-3), -5)

   @test rising_factorial(ZZ(3), ZZ(2)) == 12
   @test rising_factorial(ZZ(3), ZZ(0)) == 1
   @test rising_factorial(ZZ(0), ZZ(3)) == 0
   @test rising_factorial(ZZ(-3), ZZ(3)) == -6

   @test rising_factorial(ZZ(2), ZZ(3)) isa Oscar.fmpz
   @test rising_factorial(ZZ(2), ZZ(0)) isa Oscar.fmpz

   @test_throws DomainError rising_factorial(ZZ(0), ZZ(-1))
   @test_throws DomainError rising_factorial(ZZ(-3), ZZ(-5))

   @test primorial(6) == 30
   @test primorial(5) == 30
   @test primorial(2) == 2
   @test primorial(1) == 1
   @test primorial(0) == 1

   @test primorial(0) isa Int
   @test primorial(1) isa Int
   @test primorial(2) isa Int

   @test_throws DomainError primorial(-1)

   @test primorial(ZZ(6)) == 30
   @test primorial(ZZ(5)) == 30
   @test primorial(ZZ(2)) == 2
   @test primorial(ZZ(1)) == 1
   @test primorial(ZZ(0)) == 1

   @test primorial(ZZ(0)) isa fmpz
   @test primorial(ZZ(1)) isa fmpz
   @test primorial(ZZ(2)) isa fmpz

   @test_throws DomainError primorial(ZZ(-1))

   @test bell(0) == 1
   @test bell(1) == 1
   @test bell(3) == 5

   @test bell(0) isa Int
   @test bell(1) isa Int
   @test bell(3) isa Int

   @test_throws DomainError bell(-1)

   @test bell(ZZ(0)) == 1
   @test bell(ZZ(1)) == 1
   @test bell(ZZ(3)) == 5

   @test bell(ZZ(0)) isa Oscar.fmpz
   @test bell(ZZ(1)) isa Oscar.fmpz
   @test bell(ZZ(3)) isa Oscar.fmpz

   @test_throws DomainError bell(ZZ(-1))

   @test binomial(ZZ(3), ZZ(2)) == 3
   @test binomial(ZZ(3), ZZ(0)) == 1
   @test binomial(ZZ(3), ZZ(3)) == 1
   @test binomial(ZZ(0), ZZ(0)) == 1
   @test binomial(ZZ(1), ZZ(0)) == 1
   @test binomial(ZZ(1), ZZ(1)) == 1

   @test binomial(ZZ(2), ZZ(-3)) == 0
   @test binomial(ZZ(0), ZZ(-3)) == 0
   @test binomial(ZZ(-1), ZZ(-1)) == 0
   @test binomial(ZZ(3), ZZ(5)) == 0

   @test binomial(ZZ(3), ZZ(2)) isa Oscar.fmpz
   @test binomial(ZZ(0), ZZ(0)) isa Oscar.fmpz
   @test binomial(ZZ(-3), ZZ(2)) isa Oscar.fmpz
   @test binomial(ZZ(2), ZZ(-3)) isa Oscar.fmpz
   @test binomial(ZZ(3), ZZ(5)) isa Oscar.fmpz

   @test number_of_partitions(0) == 1
   @test number_of_partitions(ZZ(0)) == 1
   @test number_of_partitions(3) == 3
   @test number_of_partitions(ZZ(3)) == 3
   @test number_of_partitions(-1) == 0
   @test number_of_partitions(ZZ(-1)) == 0

   @test number_of_partitions(0) isa Int
   @test number_of_partitions(ZZ(0)) isa Oscar.fmpz
   @test number_of_partitions(3) isa Int
   @test number_of_partitions(ZZ(3)) isa Oscar.fmpz
   @test number_of_partitions(-1) isa Int
   @test number_of_partitions(ZZ(-1)) isa Oscar.fmpz
end

@testset "Rings.ZZ.number_theoretical" begin
   @test fibonacci(1) == 1
   @test fibonacci(2) == 1
   @test fibonacci(3) == 2
   @test fibonacci(0) == 0
   @test fibonacci(-1) == 1
   @test fibonacci(-2) == -1

   @test fibonacci(1) isa Int
   @test fibonacci(2) isa Int
   @test fibonacci(3) isa Int
   @test fibonacci(0) isa Int
   @test fibonacci(-1) isa Int

   @test fibonacci(ZZ(1)) == 1
   @test fibonacci(ZZ(2)) == 1
   @test fibonacci(ZZ(3)) == 2
   @test fibonacci(ZZ(0)) == 0
   @test fibonacci(ZZ(-1)) == 1
   @test fibonacci(ZZ(-2)) == -1

   @test fibonacci(ZZ(1)) isa Oscar.fmpz
   @test fibonacci(ZZ(2)) isa Oscar.fmpz
   @test fibonacci(ZZ(3)) isa Oscar.fmpz
   @test fibonacci(ZZ(0)) isa Oscar.fmpz
   @test fibonacci(ZZ(-1)) isa Oscar.fmpz

   @test moebius_mu(1) == 1
   @test moebius_mu(2) == -1
   @test moebius_mu(6) == 1

   @test moebius_mu(ZZ(1)) == 1
   @test moebius_mu(ZZ(2)) == -1
   @test moebius_mu(ZZ(6)) == 1

   @test moebius_mu(1) isa Int
   @test moebius_mu(2) isa Int

   @test moebius_mu(ZZ(1)) isa Int
   @test moebius_mu(ZZ(2)) isa Int

   @test_throws DomainError moebius_mu(0)
   @test_throws DomainError moebius_mu(ZZ(0))
   @test_throws DomainError moebius_mu(-1)
   @test_throws DomainError moebius_mu(ZZ(-1))

   @test jacobi_symbol(2, 3) == -1
   @test jacobi_symbol(1, 3) == 1
   @test jacobi_symbol(0, 3) == 0
   @test jacobi_symbol(-1, 3) == -1
   @test jacobi_symbol(14, 15) == -1
   @test jacobi_symbol(1, 15) == 1
   @test jacobi_symbol(5, 15) == 0
   @test jacobi_symbol(2, 45) == -1

   @test jacobi_symbol(ZZ(2), ZZ(3)) == -1
   @test jacobi_symbol(ZZ(1), ZZ(3)) == 1
   @test jacobi_symbol(ZZ(0), ZZ(3)) == 0
   @test jacobi_symbol(ZZ(-1), ZZ(3)) == -1
   @test jacobi_symbol(ZZ(14), ZZ(15)) == -1
   @test jacobi_symbol(ZZ(1), ZZ(15)) == 1
   @test jacobi_symbol(ZZ(5), ZZ(15)) == 0
   @test jacobi_symbol(ZZ(2), ZZ(45)) == -1

   @test jacobi_symbol(2, 3) isa Int
   @test jacobi_symbol(1, 3) isa Int
   @test jacobi_symbol(0, 3) isa Int
   @test jacobi_symbol(-1, 3) isa Int

   @test jacobi_symbol(ZZ(2), ZZ(3)) isa Int
   @test jacobi_symbol(ZZ(1), ZZ(3)) isa Int
   @test jacobi_symbol(ZZ(0), ZZ(3)) isa Int
   @test jacobi_symbol(ZZ(-1), ZZ(3)) isa Int

   @test_throws DomainError jacobi_symbol(2, 0)
   @test_throws DomainError jacobi_symbol(ZZ(2), ZZ(0))
   @test_throws DomainError jacobi_symbol(2, 2)
   @test_throws DomainError jacobi_symbol(ZZ(2), ZZ(2))
   @test_throws DomainError jacobi_symbol(2, -1)
   @test_throws DomainError jacobi_symbol(ZZ(2), ZZ(-1))

   @test divisor_sigma(6, 0) == 4
   @test divisor_sigma(6, 1) == 12
   @test divisor_sigma(6, 2) == 50

   @test divisor_sigma(ZZ(6), 0) == 4
   @test divisor_sigma(ZZ(6), 1) == 12
   @test divisor_sigma(ZZ(6), 2) == 50

   @test divisor_sigma(ZZ(6), ZZ(0)) == 4
   @test divisor_sigma(ZZ(6), ZZ(1)) == 12
   @test divisor_sigma(ZZ(6), ZZ(2)) == 50

   @test divisor_sigma(6, 0) isa Int
   @test divisor_sigma(6, 1) isa Int
   @test divisor_sigma(6, 2) isa Int

   @test divisor_sigma(ZZ(6), 0) isa Oscar.fmpz
   @test divisor_sigma(ZZ(6), 1) isa Oscar.fmpz
   @test divisor_sigma(ZZ(6), 2) isa Oscar.fmpz

   @test divisor_sigma(ZZ(6), ZZ(0)) isa Oscar.fmpz
   @test divisor_sigma(ZZ(6), ZZ(1)) isa Oscar.fmpz
   @test divisor_sigma(ZZ(6), ZZ(2)) isa Oscar.fmpz

   @test_throws DomainError divisor_sigma(0, 1)
   @test_throws DomainError divisor_sigma(ZZ(0), 1)
   @test_throws DomainError divisor_sigma(-1, 1)
   @test_throws DomainError divisor_sigma(ZZ(-1), 1)
   @test_throws DomainError divisor_sigma(6, -1)
   @test_throws DomainError divisor_sigma(ZZ(6), -1)
   @test_throws DomainError divisor_sigma(ZZ(6), ZZ(-1))

   @test euler_phi(1) == 1
   @test euler_phi(2) == 1

   @test euler_phi(ZZ(1)) == 1
   @test euler_phi(ZZ(2)) == 1

   @test euler_phi(1) isa Int
   @test euler_phi(2) isa Int

   @test euler_phi(ZZ(1)) isa Oscar.fmpz
   @test euler_phi(ZZ(2)) isa Oscar.fmpz

   @test_throws DomainError euler_phi(0)
   @test_throws DomainError euler_phi(ZZ(0))
   @test_throws DomainError euler_phi(-1)
   @test_throws DomainError euler_phi(ZZ(-1))
end
