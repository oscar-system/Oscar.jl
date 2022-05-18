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
   @test is_unit(QQ(2, 3))
   @test is_unit(QQ(-1, 3))
   @test !is_unit(zero(QQ))

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

   @test height(QQ(-2, 3)) isa Oscar.fmpz
   @test height(QQ()) isa Oscar.fmpz
   @test height(QQ(2)) isa Oscar.fmpz

   @test height(QQ(-3, 2)) == 3
   @test height(QQ(-2, 3)) == 3
   @test height(QQ(7, 5)) == 7
   @test height(QQ(5, 7)) == 7
   @test height(QQ()) == 1
   @test height(QQ(-1)) == 1
   @test height(QQ(2)) == 2

   @test floor(QQ(-2, 3)) isa Oscar.fmpq
   @test floor(QQ(7, 5)) isa Oscar.fmpq
   @test floor(QQ(2)) isa Oscar.fmpq
   @test floor(QQ()) isa Oscar.fmpq

   @test ceil(QQ(-2, 3)) isa Oscar.fmpq
   @test ceil(QQ(7, 5)) isa Oscar.fmpq
   @test ceil(QQ(2)) isa Oscar.fmpq
   @test ceil(QQ()) isa Oscar.fmpq

   @test floor(QQ(-1, 2)) == -1
   @test floor(QQ(1, 2)) == 0
   @test floor(QQ(1)) == 1
   @test floor(QQ(0)) == 0
   @test floor(QQ(-1)) == -1

   @test floor(QQ(-1, 2)) == -1
   @test floor(QQ(1, 2)) == 0
   @test floor(QQ(1)) == 1
   @test floor(QQ(0)) == 0
   @test floor(QQ(-1)) == -1
end

@testset "Rings.QQ.arithmetic" begin
   @test QQ(2, 3) + 1 isa Oscar.fmpq
   @test QQ(2, 3) - 1 isa Oscar.fmpq
   @test QQ(2, 3)*2 isa Oscar.fmpq
   @test QQ(1, 2) + 12345678901234567890 isa Oscar.fmpq
   @test QQ(1, 2) + 1234567890123456789012345678901234567890 isa Oscar.fmpq
   @test QQ(1, 2) - 12345678901234567890 isa Oscar.fmpq
   @test QQ(1, 2) - 1234567890123456789012345678901234567890 isa Oscar.fmpq
   @test QQ(1, 2)*12345678901234567890 isa Oscar.fmpq
   @test QQ(1, 2)*1234567890123456789012345678901234567890 isa Oscar.fmpq

   @test 1 + QQ(2, 3) isa Oscar.fmpq
   @test 1 - QQ(2, 3) isa Oscar.fmpq
   @test 2*QQ(2, 3) isa Oscar.fmpq
   @test 12345678901234567890 + QQ(1, 2) isa Oscar.fmpq
   @test 1234567890123456789012345678901234567890 + QQ(1, 2) isa Oscar.fmpq
   @test 12345678901234567890 - QQ(1, 2) isa Oscar.fmpq
   @test 1234567890123456789012345678901234567890 - QQ(1, 2) isa Oscar.fmpq
   @test 12345678901234567890*QQ(1, 2) isa Oscar.fmpq
   @test 1234567890123456789012345678901234567890*QQ(1, 2) isa Oscar.fmpq

   @test QQ(2, 3) + 2//5 isa Oscar.fmpq
   @test QQ(2, 3) - 2//5 isa Oscar.fmpq
   @test QQ(2, 3)*(2//5) isa Oscar.fmpq
   @test QQ(2, 3) + 12345678901234567890//5 isa Oscar.fmpq
   @test QQ(2, 3) - 12345678901234567890//5 isa Oscar.fmpq
   @test QQ(2, 3)*(12345678901234567890//5) isa Oscar.fmpq
   @test QQ(2, 3) + BigInt(2)//5 isa Oscar.fmpq
   @test QQ(2, 3) - BigInt(2)//5 isa Oscar.fmpq
   @test QQ(2, 3)*(BigInt(2)//5) isa Oscar.fmpq

   @test 2//5 + QQ(2, 3) isa Oscar.fmpq
   @test 2//5 - QQ(2, 3) isa Oscar.fmpq
   @test (2//5)*QQ(2, 3) isa Oscar.fmpq
   @test 12345678901234567890//5 + QQ(2, 3) isa Oscar.fmpq
   @test 12345678901234567890//5 - QQ(2, 3) isa Oscar.fmpq
   @test (12345678901234567890//5)*QQ(2, 3) isa Oscar.fmpq
   @test BigInt(2)//5 + QQ(2, 3) isa Oscar.fmpq
   @test BigInt(2)//5 - QQ(2, 3) isa Oscar.fmpq
   @test (BigInt(2)//5)*QQ(2, 3) isa Oscar.fmpq

   a = QQ(2, 3)

   @test a !== a + 0
   @test a !== a + BigInt(0)
   @test a !== a + ZZ(0)
   @test a !== a + QQ(0)
   @test a !== a + 0//1
   @test a !== a + BigInt(0)//1
   @test a !== 0 + a
   @test a !== BigInt(0) + a
   @test a !== ZZ(0) + a
   @test a !== QQ(0) + a
   @test a !== 0//1 + a
   @test a !== BigInt(0)//1 + a

   @test a !== a - 0
   @test a !== a - BigInt(0)
   @test a !== a - ZZ(0)
   @test a !== a - QQ(0)
   @test a !== a - 0//1
   @test a !== a - BigInt(0)//1
   

   @test a !== a*1
   @test a !== a*BigInt(1)
   @test a !== a*ZZ(1)
   @test a !== a*QQ(1)
   @test a !== a*(1//1)
   @test a !== a*(BigInt(1)//1)
   @test a !== 1*a
   @test a !== BigInt(1)*a
   @test a !== ZZ(1)*a
   @test a !== QQ(1)*a
   @test a !== (1//1)*a
   @test a !== (BigInt(1)//1)*a
end

@testset "Rings.QQ.comparison" begin
   @test QQ(2, 3) == QQ(2, 3)

   @test QQ(2, 3) == 2//3
   @test QQ(2, 3) == BigInt(2)//3

   @test QQ(2) == 2
   @test QQ(2) == BigInt(2)
   @test QQ(2) == ZZ(2)

   @test 2//3 == QQ(2, 3)
   @test BigInt(2)//3 == QQ(2, 3)

   @test 2 == QQ(2)
   @test BigInt(2) == QQ(2)
   @test ZZ(2) == QQ(2)

   @test QQ(2, 3) > 0
   @test QQ(2, 3) > BigInt(0)
   @test QQ(2, 3) > 0//1
   @test QQ(2, 3) > BigInt(0)//1
   @test QQ(2, 3) > ZZ(0)
   @test QQ(2, 3) > QQ(0)

   @test QQ(2, 3) < 1
   @test QQ(2, 3) < BigInt(1)
   @test QQ(2, 3) < 1//1
   @test QQ(2, 3) < BigInt(1)//1
   @test QQ(2, 3) < ZZ(1)
   @test QQ(2, 3) < QQ(1)
end

@testset "Rings.QQ.divexact" begin
   @test divexact(QQ(2, 3), QQ(3, 2)) isa Oscar.fmpq
   @test divexact(QQ(2, 3), 3//2) isa Oscar.fmpq
   @test divexact(QQ(2, 3), BigInt(3)//2) isa Oscar.fmpq
   @test divexact(2//3, QQ(2, 3)) isa Oscar.fmpq
   @test divexact(BigInt(2)//3, QQ(2, 3)) isa Oscar.fmpq

   @test divexact(QQ(2, 3), 2) isa Oscar.fmpq
   @test divexact(QQ(2, 3), BigInt(2)) isa Oscar.fmpq
   @test divexact(QQ(2, 3), ZZ(2)) isa Oscar.fmpq
   @test divexact(2, QQ(2, 3)) isa Oscar.fmpq
   @test divexact(BigInt(2), QQ(2, 3)) isa Oscar.fmpq
   @test divexact(ZZ(2), QQ(2, 3)) isa Oscar.fmpq

   @test divexact(QQ(2, 3), QQ(2, 3)) == 1
   @test divexact(QQ(2, 3), 5) == 2//15
   @test divexact(QQ(2, 3), ZZ(5)) == 2//15
   @test divexact(5, QQ(2, 3)) == 15//2
   @test divexact(ZZ(5), QQ(2, 3)) == 15//2

   @test_throws DivideError divexact(QQ(2, 3), 0)
   @test_throws DivideError divexact(QQ(2, 3), BigInt(0))
   @test_throws DivideError divexact(QQ(2, 3), ZZ(0))
   @test_throws DivideError divexact(QQ(2, 3), 0//1)
   @test_throws DivideError divexact(QQ(2, 3), BigInt(0)//1)
   @test_throws DivideError divexact(QQ(2, 3), QQ(0))
end

@testset "Rings.QQ.powering" begin
   @test QQ(2, 3)^0 isa Oscar.fmpq
   @test QQ(2, 3)^1 isa Oscar.fmpq
   @test QQ(2, 3)^2 isa Oscar.fmpq
   @test QQ(2, 3)^-1 isa Oscar.fmpq
   @test QQ(2, 3)^-1 isa Oscar.fmpq

   @test QQ(2, 3)^0 == 1
   @test QQ(2, 3)^1 == 2//3
   @test QQ(2, 3)^2 == 4//9
   @test QQ(2, 3)^-1 == 3//2
   @test QQ(2, 3)^-2 == 9//4

   @test QQ(1)^0 == 1
   @test QQ(1)^1 == 1
   @test QQ(-1)^0 == 1
   @test QQ(-1)^1 == -1
   @test QQ(-1)^2 == 1
   @test QQ(-1)^-1 == -1
   @test QQ(-1)^-2 == 1
   @test QQ(2)^3 == 8
   @test QQ(2)^-3 == 1//8

   @test QQ(0)^1 == 0
   @test QQ(0)^0 == 1
      
   @test_throws DivideError QQ(0)^-1
   @test_throws DivideError QQ(0)^-2
end
