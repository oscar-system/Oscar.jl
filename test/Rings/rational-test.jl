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

   @test QQ(2//3) == 2//3
   @test QQ(-2//3) == -2//3
   @test QQ(0//3) == 0

   @test QQ(BigInt(2)//3) == 2//3

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

   @test QQ(2//typemin(Int)) == 1//(div(typemin(Int), 2))

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
   @test iszero(QQ())
   @test !iszero(QQ(-1, 3))
   @test isone(one(QQ))
   @test !isone(-one(QQ))
   @test isunit(QQ(2, 3))
   @test isunit(QQ(-1, 3))
   @test !isunit(zero(QQ))

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

   @test numerator(QQ(-2, 3)) isa Oscar.fmpz
   @test denominator(QQ(-2, 3)) isa Oscar.fmpz

   @test numerator(QQ()) isa Oscar.fmpz
   @test denominator(QQ()) isa Oscar.fmpz

   @test numerator(QQ(1)) isa Oscar.fmpz
   @test denominator(QQ(1)) isa Oscar.fmpz

   @test sign(QQ(-2, 3)) == -1
   @test sign(QQ()) == 0
   @test sign(QQ(2, 3)) == 1
   @test sign(QQ(1)) == 1
   @test sign(QQ(-1)) == -1
   @test sign(QQ(-2)) == -1
   @test sign(QQ(0)) == 0

   @test sign(QQ(-2, 3)) isa Oscar.fmpq
   @test sign(QQ()) isa Oscar.fmpq
   @test sign(QQ(1)) isa Oscar.fmpq
   @test sign(QQ(-1)) isa Oscar.fmpq

   @test abs(QQ(-2, 3)) == 2//3
   @test abs(QQ()) == 0
   @test abs(QQ(2)) == 2

   @test abs(QQ(-2, 3)) isa Oscar.fmpq
   @test abs(QQ()) isa Oscar.fmpq
   @test abs(QQ(2)) isa Oscar.fmpq
end
