@testset "Rings.QQ.constructors" begin
   @test QQ === FractionField(ZZ)

   @test QQ(1, 2) isa Oscar.fmpq
   @test QQ(1, ZZ(2)) isa Oscar.fmpq
   @test QQ(ZZ(1), 2) isa Oscar.fmpq
   @test QQ(ZZ(1), ZZ(2)) isa Oscar.fmpq
   @test QQ(2) isa Oscar.fmpq
   @test QQ(ZZ(2)) isa Oscar.fmpq
   @test QQ() isa Oscar.fmpq

   @test 1//ZZ(2) isa Oscar.fmpq
   @test ZZ(1)//2 isa Oscar.fmpq
   @test ZZ(1)//ZZ(2) isa Oscar.fmpq

   @test zero(QQ) isa Oscar.fmpq
   @test one(QQ) isa Oscar.fmpq

   @test QQ(1, -2) == -1//2
   @test QQ(-1, 2) == -1//2
   @test QQ(-1, -2) == 1//2
   @test QQ(ZZ(1), -2) == -1//2
   @test QQ(ZZ(-1), 2) == -1//2
   @test QQ(ZZ(-1), -2) == 1//2
   @test QQ(1, ZZ(-2)) == -1//2
   @test QQ(-1, ZZ(2)) == -1//2
   @test QQ(-1, ZZ(-2)) == 1//2
   @test QQ(ZZ(1), ZZ(-2)) == -1//2
   @test QQ(ZZ(-1), ZZ(2)) == -1//2
   @test QQ(ZZ(-1), ZZ(-2)) == 1//2

   @test QQ(0, 2) == 0
   @test QQ(0, -1) == 0

   @test QQ(3) == 3
   @test QQ(0) == 0

   @test QQ() == 0

   @test ZZ(1)//(-2) == -1//2
   @test ZZ(-1)//2 == -1//2
   @test ZZ(-1)//-2 == 1//2
   @test 1//ZZ(-2) == -1//2
   @test -1//ZZ(2) == -1//2
   @test -1//ZZ(-2) == 1//2
   @test ZZ(1)//ZZ(-2) == -1//2
   @test ZZ(-1)//ZZ(2) == -1//2
   @test ZZ(-1)//ZZ(-2) == 1//2

   @test 0//ZZ(2) == 0
   @test 0//ZZ(-1) == 0

   @test QQ(1, typemax(Int)) == 1//typemax(Int)
   @test QQ(typemax(Int), 1) == typemax(Int)
   @test QQ(1, typemin(Int)) == 1//BigInt(typemin(Int)) # must use BigInt
   @test QQ(typemin(Int), 1) == typemin(Int)
   @test QQ(typemin(Int), -1) == -BigInt(typemin(Int)) # must use BigInt
   @test QQ(typemax(Int), -1) == -BigInt(typemax(Int)) # must use BigInt
   @test QQ(-1, typemin(Int)) == 1//-BigInt(typemin(Int)) # must use BigInt
   @test QQ(-1, typemax(Int)) == -1//typemax(Int)

   @test_throws DivideError ZZ(2)//0
   @test_throws DivideError 2//ZZ(0)
   @test_throws DivideError ZZ(2)//ZZ(0)
   @test_throws DivideError ZZ(0)//0
   @test_throws DivideError 0//ZZ(0)
   @test_throws DivideError ZZ(0)//ZZ(0)

   @test_throws DivideError QQ(2, 0)
   @test_throws DivideError QQ(ZZ(2), 0)
   @test_throws DivideError QQ(2, ZZ(0))
   @test_throws DivideError QQ(ZZ(2), ZZ(0))
   @test_throws DivideError QQ(0, 0)
   @test_throws DivideError QQ(ZZ(0), 0)
   @test_throws DivideError QQ(0, ZZ(0))
   @test_throws DivideError QQ(ZZ(0), ZZ(0))
end

@testset "Rings.QQ.properties" begin
   @test numerator(QQ(2, 1)) == 2
   @test denominator(QQ(2, 1)) == 1

   @test numerator(QQ(6, 3)) == 2
   @test denominator(QQ(6, 3)) == 1

   @test numerator(QQ(-2, 3)) == -2
   @test denominator(QQ(-2, 3)) == 3

   @test numerator(QQ(2, -3)) == -2
   @test denominator(QQ(2, -3)) == 3

   @test numerator(QQ(-2, -3)) == 2
   @test denominator(QQ(-2, -3)) == 3

   @test numerator(QQ()) == 0
   @test denominator(QQ()) == 1
end
