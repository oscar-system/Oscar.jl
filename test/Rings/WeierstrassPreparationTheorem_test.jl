@testset "Weierstrass Preparation Theorem" begin 

	@testset "order of generality test" begin 
		R, (x,y,z)=PolynomialRing(QQ,["x","y","z"])
		f=x+y
		@test (true,1)==order_of_generality(f,2)
		@test (false,-1)==order_of_generality(f,3)
		f=x^3*y^4+3*z^2*y^3-x^2*y-x^4+x^6+z^5
		@test (true,4)==order_of_generality(f,1)
		@test (false,-1)==order_of_generality(f,2)
		@test (true,5)==order_of_generality(f,3)
	end
	
	@testset "Weierstrass Division Theorem test" begin
		R,(x,y,z)=PolynomialRing(QQ,["x","y","z"])
		f=x+y+z
		@test (1,-x-y,2)==weierstrass_division_theorem(z,f,2,3)
		f=x^2*y - 2*x^2*z^2 + 3*x*y*z + z^5 + z^3
		@test (3*x*y-z^2+1,0,3)==weierstrass_division_theorem(z^3,f,2,3)
		f=x+y
		@test (1,-x,2)==weierstrass_division_theorem(y,f,1,2)
		f=x^2*y^4*z^3 - x^2*z^3 + y^5 + y^3*z^2 - y^3
		@test (-y^4 - y^2*z^2 - y^2, 0, 2)==weierstrass_division_theorem(y^5,f,4,2)
	end
	
	@testset "Weierstrass Preparation Theorem test" begin
		R,(x,y,z)=PolynomialRing(QQ,["x","y","z"])
		f=x+y+z
		@test (1,x+y+z,2)==weierstrass_preparation_theorem(f,3,3)
		@test (1,x+y+z,2)==weierstrass_preparation_theorem(f,3,1)
		f=x^2*y^3*z^4-y^2*z^3+z^4+y^3
		@test (-x^2*y^3+1,y^3-y^2*z^3+z^4,2)==weierstrass_preparation_theorem(f,6,3)
		@test (-z^4*x^2+1,z^4-z^3*y^2+y^3,2)==weierstrass_preparation_theorem(f,6,2)
		f=x^3*y^4*z^5+x^3-y^4+z^2-x*z^3+y*z-z^3+x^3*z^5-y^3*z^2
		@test (-x^3 + 4*x*y^2 - 4*x*y*z - x*y + 2*x*z^2 + x*z - 4*y^3 + 5*y^2*z + 2*y^2 - 3*y*z^2 - 2*y*z - y + z^3 + z^2 + z + 1, x^3 - y^2*z + y*z + z^2, 5)==weierstrass_preparation_theorem(f,3,3)
	end
	
	@testset "var general test" begin
		R,(x,y,z)=PolynomialRing(QQ,["x","y","z"])
		f=x*y*z
		@test order_of_generality(var_general(f,1),1)[1]
		@test order_of_generality(var_general(f,3),3)[1]
		f=3*x*y-z^2+x^4*y^5*z^2-3*x*y*z
		@test order_of_generality(var_general(f,1),1)[1]
		@test order_of_generality(var_general(f,3),3)[1]
		
	end
	
	@testset "invert unit test" begin
		R,(x,y,z)=PolynomialRing(QQ,["x","y","z"])
		@test (true,1)==invert_unit(R(1),1)
		f=x+y+z
		@test (false,0)==invert_unit(f,3)
		f=x+y+z+1
		@test (true, x^2 + 2*x*y + 2*x*z - x + y^2 + 2*y*z - y + z^2 - z + 1)==invert_unit(f,2)
		f=x^2+3
		@test (true, 1//27*x^4-1//9*x^2+1//3)==invert_unit(f,5) 
	end
	
	@testset "jet test" begin
		R,(x,y,z)=PolynomialRing(QQ,["x","y","z"])
		f=x+y+z
		@test f==jet(f,1)
		f=x^2*y^3*z-x*z+y^7
		@test -x*z==jet(f, 4)
		@test R(0)==jet(f,0)
	end
	
	@testset "jet with weighted vector test" begin
		R,(x,y,z)=PolynomialRing(QQ,["x","y","z"])
		f=x+y+z
		vec=[2,4,0]
		@test z==jet(f, 1, vec)
		@test x+z==jet(f, 3, vec)
		f=x^3*y^2+y^2*z^3+x*y*z
		vec=[2,3,4]
		@test x^3*y^2+x*y*z==jet(f, 12, vec)
		
	end
end
