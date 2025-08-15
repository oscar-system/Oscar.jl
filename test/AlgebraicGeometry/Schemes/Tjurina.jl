@testset "global Tjurina number" begin
  R, (x, y) = QQ[:x, :y];
  @test 2 == tjurina_number(y^2 - x^3)
  @test 4 == tjurina_number(change_coefficient_ring(GF(2), y^2 - x^3))
  @test 3 == tjurina_number((y^2 - x^2)*(x-1))
  @test PosInf() == tjurina_number((x-1)^2)
  @test 0 == tjurina_number(AffineScheme(quo(R, ideal(R, one(R)))[1]))
  @test 10 == tjurina_number(AffineScheme(quo(R, ideal(R, x^5 + x^2*y^2 + y^5))[1]))
  @test 3 == tjurina_number(AffineScheme(quo(R, ideal(R, x^2 - y^4))[1]))
  @test PosInf() == tjurina_number(AffineScheme(quo(R, ideal(R, zero(R)))[1]))
end

@testset "local Tjurina number" begin
  R, (x, y) = QQ[:x, :y];  
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0]))
  @test 1 == tjurina_number(L((y^2 - x^2)*(x-1)))
  @test 0 == tjurina_number(L((x-3)^9-(y+5)^7))
  @test 7 == tjurina_number(L((x^6 + x*y^2)*(x+23)^12*(y-6)^12))
  @test PosInf() == tjurina_number(L(x^2))
  @test 12 == tjurina_number(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^4 + y^5))[1]), [0, 0]))
  @test 0 == tjurina_number(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, one(R)))[1]), [0, 0]))
  @test 4 == tjurina_number(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, (x-1)^3 + (x-1)*(y-2)^2))[1]), [1, 2]))
  @test PosInf() == tjurina_number(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^2))[1]), [0, 0]))
end

@testset "is_finitely_determined" begin
  R, (x,y) = QQ[:x, :y]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0]))
  L1,_  = localization(R, complement_of_point_ideal(R, [1, 1]))
  @test_throws ErrorException("Equivalence typ must be ':right' or ':contact'.") is_finitely_determined(L(0), :leftright)
  @test !is_finitely_determined(L(0))
  @test !is_finitely_determined(L(0), :right)  
  @test is_finitely_determined(L(x^2+y))
  @test is_finitely_determined(L(x^2+y), :right)
  @test !is_finitely_determined(L(x^2))
  @test !is_finitely_determined(L(x^2), :right)  
  @test is_finitely_determined(L(x^2+y^2))
  @test is_finitely_determined(L(x^2+y^2), :right)   
  @test is_finitely_determined(L(1))
  @test !is_finitely_determined(L(1), :right)    
  @test is_finitely_determined(L(x^2+1)) 
  @test !is_finitely_determined(L(x^2+1), :right) 
  @test is_finitely_determined(L(x^2+y+1))
  @test is_finitely_determined(L(x^2+y+1), :right)
  @test is_finitely_determined(L(x^2+y^2+1))
  @test is_finitely_determined(L(x^2+y^2+1), :right)
  @test is_finitely_determined(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+y^2))[1]), [0, 0]))  
  @test is_finitely_determined(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+y^2))[1]), [0, 0]), :right)  
  @test is_finitely_determined(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+x*y^3))[1]), [0, 1]))  
  @test is_finitely_determined(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+x*y^3))[1]), [0, 1]), :right)  
  @test is_finitely_determined(L1((x-1)^2+(y-1)^2))
  @test is_finitely_determined(L1((x-1)^2+(y-1)^2), :right)  
  @test is_finitely_determined(L1((x-1)^2+(y-1)^2+1))
  @test is_finitely_determined(L1((x-1)^2+(y-1)^2+1), :right) 
end

@testset "is_finitely_determined positive characteristic" begin
  R, (x,y) = GF(2)[:x, :y]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0]))
  @test !is_finitely_determined(L(0))
  @test !is_finitely_determined(L(0), :right)
  @test is_finitely_determined(L(1))
  @test !is_finitely_determined(L(1), :right)
  @test is_finitely_determined(L(x))
  @test is_finitely_determined(L(x), :right) 
  @test is_finitely_determined(L(x^3+y^2))
  @test !is_finitely_determined(L(x^3+y^2), :right)  
  @test !is_finitely_determined(L(x^2+y^2))
  @test !is_finitely_determined(L(x^2+y^2), :right)
end

@testset "determinacy_bound" begin
  R, (x,y) = QQ[:x, :y]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0]))
  @test_throws ErrorException("Equivalence typ must be ':right' or ':contact'.") determinacy_bound(L(0), :leftright)
  @test PosInf() == determinacy_bound(L(0)) 
  @test PosInf() == determinacy_bound(L(0), :right)  
  @test 1 == determinacy_bound(L(x^2+y))
  @test 1 == determinacy_bound(L(x^2+y), :right)
  @test PosInf() == determinacy_bound(L(x^2))
  @test PosInf() == determinacy_bound(L(x^2), :right)  
  @test 2 == determinacy_bound(L(x^2+y^2))
  @test 2 == determinacy_bound(L(x^2+y^2), :right)   
  @test 0 == determinacy_bound(L(1))
  @test PosInf() == determinacy_bound(L(1), :right)   
  @test 0 == determinacy_bound(L(x^2+y+1))
  @test 1 == determinacy_bound(L(x^2+y+1), :right)
  @test 0 == determinacy_bound(L(x^2+y^2+1))
  @test 2 == determinacy_bound(L(x^2+y^2+1), :right)
  @test 3 == determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+y^2))[1]), [0, 0]))  
  @test 3 == determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+y^2))[1]), [0, 0]), :right)  
  @test 8 == determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+x*y^3))[1]), [0, 0]))  
  @test 8 == determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+x*y^3))[1]), [0, 0]), :right) 
  @test 1 == determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+x*y^3))[1]), [0, 1]))  
  @test 1 == determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+x*y^3))[1]), [0, 1]), :right)  
  @test 17 == determinacy_bound(L(x^5+y^5))
  @test 17 == determinacy_bound(L(x^5+y^5), :right)
  @test 11 == determinacy_bound(L(x^5+x^2*y^2+y^5))
  @test 12 == determinacy_bound(L(x^5+x^2*y^2+y^5), :right)
end

@testset "determinacy_bound positive characteristic" begin
  R, (x,y) = GF(2)[:x, :y]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0]))
  @test PosInf() == determinacy_bound(L(0))
  @test PosInf() == determinacy_bound(L(0), :right)
  @test 0 == determinacy_bound(L(1))
  @test PosInf() == determinacy_bound(L(1), :right)
  @test 1 == determinacy_bound(L(x))
  @test 1 == determinacy_bound(L(x), :right) 
  @test 8 == determinacy_bound(L(x^3+y^2))
  @test PosInf() == determinacy_bound(L(x^3+y^2), :right)  
  @test PosInf() == determinacy_bound(L(x^2+y^2))
  @test PosInf() == determinacy_bound(L(x^2+y^2), :right)  
  @test 13 == determinacy_bound(L(x^3+x*y^3))
  @test 13 == determinacy_bound(L(x^3+x*y^3), :right)
  @test 22 == determinacy_bound(L(x^5+x^2*y^2+y^5))
  @test 30 == determinacy_bound(L(x^5+x^2*y^2+y^5), :right)
end

@testset "sharper_determinacy_bound" begin
  R, (x,y) = QQ[:x, :y]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0]))
  @test_throws ErrorException("Equivalence typ must be ':right' or ':contact'.") determinacy_bound(L(0), :leftright)
  @test PosInf() == sharper_determinacy_bound(L(0)) 
  @test PosInf() == sharper_determinacy_bound(L(0), :right)  
  @test 1 == sharper_determinacy_bound(L(x^2+y))
  @test 1 == sharper_determinacy_bound(L(x^2+y), :right)
  @test PosInf() == sharper_determinacy_bound(L(x^2))
  @test PosInf() == sharper_determinacy_bound(L(x^2), :right)  
  @test 3 == sharper_determinacy_bound(L(x^2+y^2))
  @test 3 == sharper_determinacy_bound(L(x^2+y^2), :right)
  @test 0 == sharper_determinacy_bound(L(1))
  @test PosInf() == sharper_determinacy_bound(L(1), :right)   
  @test 0 == sharper_determinacy_bound(L(x^2+y+1))
  @test 1 == sharper_determinacy_bound(L(x^2+y+1), :right)
  @test 0 == sharper_determinacy_bound(L(x^2+y^2+1))
  @test 3 == sharper_determinacy_bound(L(x^2+y^2+1), :right)
  @test 4 == sharper_determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+y^2))[1]), [0, 0]))
  @test 4 == sharper_determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+y^2))[1]), [0, 0]), :right)
  @test 6 == sharper_determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+x*y^3))[1]), [0, 0]))
  @test 6 == sharper_determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+x*y^3))[1]), [0, 0]), :right)
  @test 1 == sharper_determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+x*y^3))[1]), [0, 1]))
  @test 1 == sharper_determinacy_bound(HypersurfaceGerm(AffineScheme(quo(R, ideal(R, x^3+x*y^3))[1]), [0, 1]), :right)
  @test 7 == sharper_determinacy_bound(L(x^5+y^5))
  @test 7 == sharper_determinacy_bound(L(x^5+y^5), :right)
  @test 6 == sharper_determinacy_bound(L(x^5+x^2*y^2+y^5))
  @test 6 == sharper_determinacy_bound(L(x^5+x^2*y^2+y^5), :right)
end

@testset "sharper_determinacy_bound positive characteristic" begin
  R, (x,y) = GF(2)[:x, :y]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0]))
  @test PosInf() == sharper_determinacy_bound(L(0))
  @test PosInf() == sharper_determinacy_bound(L(0), :right)
  @test 0 == sharper_determinacy_bound(L(1))
  @test PosInf() == sharper_determinacy_bound(L(1), :right)
  @test 1 == sharper_determinacy_bound(L(x))
  @test 1 == sharper_determinacy_bound(L(x), :right) 
  @test 4 == sharper_determinacy_bound(L(x^3+y^2))
  @test PosInf() == sharper_determinacy_bound(L(x^3+y^2), :right)  
  @test PosInf() == sharper_determinacy_bound(L(x^2+y^2))
  @test PosInf() == sharper_determinacy_bound(L(x^2+y^2), :right)  
  @test 7 == sharper_determinacy_bound(L(x^3+x*y^3))
  @test 7 == sharper_determinacy_bound(L(x^3+x*y^3), :right)
  @test 6 == sharper_determinacy_bound(L(x^5+x^2*y^2+y^5))
  @test 8 == sharper_determinacy_bound(L(x^5+x^2*y^2+y^5), :right)
end

@testset "contact equivalence" begin
  R, (x,y) = QQ[:x, :y]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0]))
  ## order
  @test !is_contact_equivalent(L(x^5+y^5), L(x^2+y^2))
  ## unit
  @test is_contact_equivalent(L(x^3+12), L(y-12))
  ## smooth
  @test is_contact_equivalent(L(x^3-y), L(x-x*y^2))
  ## tjurina_number
  @test !is_contact_equivalent(L(x^23+y^2), L(x^6-y^2))  
  ## multiplication with unit
  @test is_contact_equivalent(L(x^3-y^2), L(x^3-y^2))
  @test is_contact_equivalent(L(x^3-x^2), L(x^2))
  @test is_contact_equivalent(L(x^6+x^2*y^3-y^4), L((x^6+x^2*y^3-y^4)*(y^3-x*y+x-1)))
  @test is_contact_equivalent(L(x^3+x*y^2+y^4), L((x^3+x*y^2+y^4)//(x^7-x*y^3+y-1)))
  ## different determinancy bound
  @test !is_contact_equivalent(L(x^5+x*y^2), L(x^3+y^4))
  @test !is_contact_equivalent(L(x^6+x*y^2), L(x^3+x*y^3))
  @test !is_contact_equivalent(L(x^7+x*y^2), L(x^3+y^5))
  ## same k-Jet
  @test is_contact_equivalent(L(x^3-y^2), L(x^5+x^3-y^2))
  @test is_contact_equivalent(L(x^5-y^4), L(x^23+x^5-x^6*y^6+y^12-y^4))
  ## Mather-Yau
  @test is_contact_equivalent(L(x^5+y^4), L(x^5*y^4-3*x^5+y^5-y^4))
  @test is_contact_equivalent(L(x^5+y^4), L(x^5-y^4))
  @test is_contact_equivalent(L(x*y), L(x^2-y^2*(x-3)))
  @test is_contact_equivalent(L(x^5+y^2), L(x^2-(x-y)^5))
  
  ## pos. char.
  R, (x,y) = GF(5)[:x, :y]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0]))
  @test !is_contact_equivalent(L(x^5+y^2), L(x^6+y^2))  
  ## not isolated
  @test_throws ErrorException("Unable to determine if is contact equivalent. (Singularities are not isolated)") is_contact_equivalent(L(x^5+x^3), L(-y^3))
end

@testset "Mather-Yau Isomorphismtest" begin
  R, (x,y) = QQ[:x, :y]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0]))  
  @test is_contact_equivalent(L(x^5+y^4), L(x^4+y^5))
  @test is_contact_equivalent(L(x^7+y^4), L((x+x*y)^4-y^7))  
  @test is_contact_equivalent(L(x^5-x^2*y^2+y^5), L(x^5+2*x^2*y^2+y^5))
  @test is_contact_equivalent(L(x^3+y^4), L(x^4+(x-y)^3))
  ## pos. char.
  R, (x,y) = GF(7)[:x, :y]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0])) 
  @test is_contact_equivalent(L(x^2+y^2),L(x*y)) 
  @test is_contact_equivalent(L(x^2+y^3),L(x^3+y^2))  
  @test is_contact_equivalent(L(x^4+y^3),L(x^3-5*y^4))
end

@testset "Mather-Yau with tjurina_number infinity" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0, 0]))  
  @test is_contact_equivalent(zero(L), L(0//(x+1)))
  @test is_contact_equivalent(L(x^2+y^2), L(y^2+z^2)) 
  @test !is_contact_equivalent(L(z^2), L(x^2+y^2))
  @test !is_contact_equivalent(L(x^2+y^2), L(x^2+y^2-z^2))
end

@testset "ADE-Singularities" begin  
  R, (x,y,z) = QQ[:x, :y, :z]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0, 0]))
  @test !is_contact_equivalent(L(x^5+y^2+z^2), L(x^3+x*y^2+z^2))
  @test !is_contact_equivalent(L(x^7+y^2+z^2), L(x^5+x*y^2))
  @test !is_contact_equivalent(L(x^7+y^2+z^2), L(x^3+y^4))
  @test !is_contact_equivalent(L(x^8+y^2+z^2), L(x^6+x*y^2))
  @test !is_contact_equivalent(L(x^8+y^2+z^2), L(x^3+x*y^3))
  @test !is_contact_equivalent(L(x^9+y^2+z^2), L(x^7+x*y^2))
  @test !is_contact_equivalent(L(x^9+y^2+z^2), L(x^3+y^5))
end

@testset "ADE-Singularities pos Char." begin  
  R, (x,y,z) = GF(7)[:x, :y, :z]
  L,_  = localization(R, complement_of_point_ideal(R, [0, 0, 0]))
  @test !is_contact_equivalent(L(x^5+y^2+z^2), L(x^3+x*y^2+z^2))
  @test !is_contact_equivalent(L(x^7+y^2+z^2), L(x^5+x*y^2))
  @test !is_contact_equivalent(L(x^7+y^2+z^2), L(x^3+y^4))
  @test !is_contact_equivalent(L(x^8+y^2+z^2), L(x^6+x*y^2))
  @test !is_contact_equivalent(L(x^8+y^2+z^2), L(x^3+x*y^3))
  @test !is_contact_equivalent(L(x^9+y^2+z^2), L(x^7+x*y^2))
  @test !is_contact_equivalent(L(x^9+y^2+z^2), L(x^3+y^5))
end



@testset "Tjurina number complete intersection germ" begin  
  A = affine_space(QQ, 3)
  R = coordinate_ring(A);
  (x,y,z) = gens(R);
  X = CompleteIntersectionGerm(spec(quo(R,ideal(R,[x^2+y^2+z^2, x^2+2*y^2+3*z^2]))[1]), [0,0,0])
  @test tjurina_number(X) == 5
  S = spec(quo(R,ideal(R,[x^5+y^6+z^7+x*y*z]))[1])
  @test tjurina_number(HypersurfaceGerm(S, [0,0,0])) == tjurina_number(CompleteIntersectionGerm(S, [0,0,0]))
  Y = CompleteIntersectionGerm(spec(quo(R,ideal(R,[x^2+y^2, x^2+2*y^2]))[1]), [0,0,0])
  @test tjurina_number(Y) == PosInf()  
end

