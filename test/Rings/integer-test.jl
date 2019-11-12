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
   @test isunit(one(ZZ))
   @test isunit(-one(ZZ))
   @test !isunit(zero(ZZ))
   @test sign(ZZ(2)) == 1
   @test sign(ZZ(0)) == 0
   @test sign(-ZZ(2)) == -1
   @test sign(ZZ(2)) isa Int
   @test sign(ZZ(0)) isa Int
   @test sign(ZZ(-2)) isa Int

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

   @test_throws ErrorException divexact(ZZ(2), ZZ(3))
   @test_throws ErrorException divexact(ZZ(2), 3)
   @test_throws ErrorException divexact(2, ZZ(3))

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
end

