###
# Linear example /QQ(s)
###
import Random
K,s = RationalFunctionField(QQ,"s");
Kx,(x1,x2,x3,x4) = PolynomialRing(K,4);
I = ideal([x1-s*x2+(s+1)*x3,3*x2-s^2*x3+(s^2+1)*x4]);
val = ValuationMap(K,s,min); Random.seed!(3847598273423); TropI = tropical_variety(I,val) # works
val = ValuationMap(K,s,max); Random.seed!(3847598273423); TropI = tropical_variety(I,val) # does not work



###
# Linear example /QQ_p
###
import Random
Kx,(x1,x2,x3,x4) = PolynomialRing(QQ,4);
p = 2;
I = ideal([x1-p*x2+(p+1)*x3,3*x2-p^2*x3+(p^2+1)*x4]);
val = ValuationMap(QQ,p,min); Random.seed!(3847598273423); TropI = tropical_variety(I,val)
val = ValuationMap(QQ,p,max); Random.seed!(3847598273423); TropI = tropical_variety(I,val)


###
# Mixed degree example, space sextic /QQ_p (t-adic valuation is too hard)
# todo: send non-terminating GB computation to Christian
###
import Random
K = QQ;
p = K(3);
Kx,(x,y,z) = PolynomialRing(K,3);
val = ValuationMap(K,p,max);
Random.seed!(1337133713371337);
I = ideal(
  [rand(1:2)*p*1+rand(1:2)*x+rand(1:2)*y+rand(1:2)*z+rand(1:2)*x*y+rand(1:2)*x*z+rand(1:2)*y*z+rand(1:2)*p*x^2+rand(1:2)*p*y^2+rand(1:2)*p*z^2,
   rand(1:2)*p^-3*1+rand(1:2)*p^-3*x^3+rand(1:2)*p^-3*y^3+rand(1:2)*p^-3*z^3+rand(1:2)*p^0*x*y+rand(1:2)*p^0*x*z+rand(1:2)*p^0*y*z+rand(1:2)*p^0*x*y*z+rand(1:2)*p^-1*x+rand(1:2)*p^-1*y+rand(1:2)*p^-1*z+rand(1:2)*p^-1*x^2+rand(1:2)*p^-1*y^2+rand(1:2)*p^-1*z^2+rand(1:2)*p^-1*x^2*y+rand(1:2)*p^-1*x^2*z+rand(1:2)*p^-1*x*y^2+rand(1:2)*p^-1*x*z^2+rand(1:2)*p^-1*y^2*z+rand(1:2)*p^-1*y*z^2]);
TropI, wG = tropical_variety(I,val)





###
# Grass(2,5) /QQ(s)
###
import Random
K,s = RationalFunctionField(QQ,"s");
Kx,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10) = PolynomialRing(K,10);
p01 = x1
p02 = x2
p03 = x3
p04 = x4
p12 = x5
p13 = x6
p14 = x7
p23 = x8
p24 = x9
p34 = x10
I = ideal([p03*p12-p02*p13+p01*p23,
           p04*p12-p02*p14+p01*p24,
           p04*p13-p03*p14+p01*p34,
           p04*p23-p03*p24+p02*p34,
           p03*p12-p02*p13+p01*p23,
           p04*p12-p02*p14+p01*p24,
           p04*p13-p03*p14+p01*p34,
           p14*p23-p13*p24+p12*p34,
           p03*p12-p02*p13+p01*p23,
           p04*p12-p02*p14+p01*p24,
           p04*p23-p03*p24+p02*p34,
           p14*p23-p13*p24+p12*p34,
           p03*p12-p02*p13+p01*p23,
           p04*p13-p03*p14+p01*p34,
           p04*p23-p03*p24+p02*p34,
           p14*p23-p13*p24+p12*p34,
           p04*p12-p02*p14+p01*p24,
           p04*p13-p03*p14+p01*p34,
           p04*p23-p03*p24+p02*p34,
           p14*p23-p13*p24+p12*p34])
val = ValuationMap(K,s,min); Random.seed!(3847598273423); TropI,wG = tropical_variety(I,val)
val = ValuationMap(K,s,max); Random.seed!(3847598273423); TropI,wG = tropical_variety(I,val)


###
# Grass(2,5) /FF_p(s)
###
import Random
K,s = RationalFunctionField(GF(32003),"s");
Kx,(p01,p02,p03,p04,p12,p13,p14,p23,p24,p34) = PolynomialRing(K,10);
I = ideal([p03*p12-p02*p13+p01*p23,
           p04*p12-p02*p14+p01*p24,
           p04*p13-p03*p14+p01*p34,
           p04*p23-p03*p24+p02*p34,
           p03*p12-p02*p13+p01*p23,
           p04*p12-p02*p14+p01*p24,
           p04*p13-p03*p14+p01*p34,
           p14*p23-p13*p24+p12*p34,
           p03*p12-p02*p13+p01*p23,
           p04*p12-p02*p14+p01*p24,
           p04*p23-p03*p24+p02*p34,
           p14*p23-p13*p24+p12*p34,
           p03*p12-p02*p13+p01*p23,
           p04*p13-p03*p14+p01*p34,
           p04*p23-p03*p24+p02*p34,
           p14*p23-p13*p24+p12*p34,
           p04*p12-p02*p14+p01*p24,
           p04*p13-p03*p14+p01*p34,
           p04*p23-p03*p24+p02*p34,
           p14*p23-p13*p24+p12*p34])
val = ValuationMap(K,s,min); Random.seed!(3847598273423); TropI = tropical_variety(I,val)
val = ValuationMap(K,s,max); Random.seed!(3847598273423); TropI = tropical_variety(I,val)


###
# Grass(2,6) /QQ(s)
###
import Random
K,s = RationalFunctionField(QQ,"s");
Kx,(p01,p02,p03,p04,p05,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45) = PolynomialRing(K,15);
val = ValuationMap(K,s);
I = ideal([p03*p12-p02*p13+p01*p23,  p04*p12-p02*p14+p01*p24,
	         p05*p12-p02*p15+p01*p25,  p04*p13-p03*p14+p01*p34,
	         p05*p13-p03*p15+p01*p35,  p05*p14-p04*p15+p01*p45,
	         p04*p23-p03*p24+p02*p34,  p05*p23-p03*p25+p02*p35,
	         p05*p24-p04*p25+p02*p45,  p05*p34-p04*p35+p03*p45,
	         p03*p12-p02*p13+p01*p23,  p04*p12-p02*p14+p01*p24,
	         p05*p12-p02*p15+p01*p25,  p04*p13-p03*p14+p01*p34,
	         p05*p13-p03*p15+p01*p35,  p05*p14-p04*p15+p01*p45,
	         p14*p23-p13*p24+p12*p34,  p15*p23-p13*p25+p12*p35,
	         p15*p24-p14*p25+p12*p45,  p15*p34-p14*p35+p13*p45,
	         p03*p12-p02*p13+p01*p23,  p04*p12-p02*p14+p01*p24,
	         p05*p12-p02*p15+p01*p25,  p04*p23-p03*p24+p02*p34,
	         p05*p23-p03*p25+p02*p35,  p05*p24-p04*p25+p02*p45,
	         p14*p23-p13*p24+p12*p34,  p15*p23-p13*p25+p12*p35,
	         p15*p24-p14*p25+p12*p45,  p25*p34-p24*p35+p23*p45,
	         p03*p12-p02*p13+p01*p23,  p04*p13-p03*p14+p01*p34,
	         p05*p13-p03*p15+p01*p35,  p04*p23-p03*p24+p02*p34,
	         p05*p23-p03*p25+p02*p35,  p05*p34-p04*p35+p03*p45,
	         p14*p23-p13*p24+p12*p34,  p15*p23-p13*p25+p12*p35,
	         p15*p34-p14*p35+p13*p45,  p25*p34-p24*p35+p23*p45,
	         p04*p12-p02*p14+p01*p24,  p04*p13-p03*p14+p01*p34,
	         p05*p14-p04*p15+p01*p45,  p04*p23-p03*p24+p02*p34,
	         p05*p24-p04*p25+p02*p45,  p05*p34-p04*p35+p03*p45,
	         p14*p23-p13*p24+p12*p34,  p15*p24-p14*p25+p12*p45,
	         p15*p34-p14*p35+p13*p45,  p25*p34-p24*p35+p23*p45,
	         p05*p12-p02*p15+p01*p25,  p05*p13-p03*p15+p01*p35,
	         p05*p14-p04*p15+p01*p45,  p05*p23-p03*p25+p02*p35,
	         p05*p24-p04*p25+p02*p45,  p05*p34-p04*p35+p03*p45,
	         p15*p23-p13*p25+p12*p35,  p15*p24-p14*p25+p12*p45,
	         p15*p34-p14*p35+p13*p45,  p25*p34-p24*p35+p23*p45])
Random.seed!(133713371337); TropI,wG = tropical_variety(I,val)



###
# Grass(2,6) /FF_p(s)
###
import Random
K,s = RationalFunctionField(GF(32003),"s");
Kx,(p01,p02,p03,p04,p05,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45) = PolynomialRing(K,15);
val = ValuationMap(K,s);
I = ideal([p03*p12-p02*p13+p01*p23,  p04*p12-p02*p14+p01*p24,
	         p05*p12-p02*p15+p01*p25,  p04*p13-p03*p14+p01*p34,
	         p05*p13-p03*p15+p01*p35,  p05*p14-p04*p15+p01*p45,
	         p04*p23-p03*p24+p02*p34,  p05*p23-p03*p25+p02*p35,
	         p05*p24-p04*p25+p02*p45,  p05*p34-p04*p35+p03*p45,
	         p03*p12-p02*p13+p01*p23,  p04*p12-p02*p14+p01*p24,
	         p05*p12-p02*p15+p01*p25,  p04*p13-p03*p14+p01*p34,
	         p05*p13-p03*p15+p01*p35,  p05*p14-p04*p15+p01*p45,
	         p14*p23-p13*p24+p12*p34,  p15*p23-p13*p25+p12*p35,
	         p15*p24-p14*p25+p12*p45,  p15*p34-p14*p35+p13*p45,
	         p03*p12-p02*p13+p01*p23,  p04*p12-p02*p14+p01*p24,
	         p05*p12-p02*p15+p01*p25,  p04*p23-p03*p24+p02*p34,
	         p05*p23-p03*p25+p02*p35,  p05*p24-p04*p25+p02*p45,
	         p14*p23-p13*p24+p12*p34,  p15*p23-p13*p25+p12*p35,
	         p15*p24-p14*p25+p12*p45,  p25*p34-p24*p35+p23*p45,
	         p03*p12-p02*p13+p01*p23,  p04*p13-p03*p14+p01*p34,
	         p05*p13-p03*p15+p01*p35,  p04*p23-p03*p24+p02*p34,
	         p05*p23-p03*p25+p02*p35,  p05*p34-p04*p35+p03*p45,
	         p14*p23-p13*p24+p12*p34,  p15*p23-p13*p25+p12*p35,
	         p15*p34-p14*p35+p13*p45,  p25*p34-p24*p35+p23*p45,
	         p04*p12-p02*p14+p01*p24,  p04*p13-p03*p14+p01*p34,
	         p05*p14-p04*p15+p01*p45,  p04*p23-p03*p24+p02*p34,
	         p05*p24-p04*p25+p02*p45,  p05*p34-p04*p35+p03*p45,
	         p14*p23-p13*p24+p12*p34,  p15*p24-p14*p25+p12*p45,
	         p15*p34-p14*p35+p13*p45,  p25*p34-p24*p35+p23*p45,
	         p05*p12-p02*p15+p01*p25,  p05*p13-p03*p15+p01*p35,
	         p05*p14-p04*p15+p01*p45,  p05*p23-p03*p25+p02*p35,
	         p05*p24-p04*p25+p02*p45,  p05*p34-p04*p35+p03*p45,
	         p15*p23-p13*p25+p12*p35,  p15*p24-p14*p25+p12*p45,
	         p15*p34-p14*p35+p13*p45,  p25*p34-p24*p35+p23*p45])
Random.seed!(133713371337); TropI,wG = tropical_variety(I,val)



###
# Grass(3,6) /FF_p(s)
###
import Random
K,s = RationalFunctionField(GF(32003),"s");
Kx,(p012,p013,p014,p015,p023,p024,p025,p034,p035,p045,p123,p124,p125,p134,p135,p145,p234,p235,p245,p345) = PolynomialRing(K,20);
val = ValuationMap(K,s);
I = ideal([p014*p023-p013*p024+p012*p034,  p015*p023-p013*p025+p012*p035,
	         p015*p024-p014*p025+p012*p045,  p015*p034-p014*p035+p013*p045,
	         p014*p123-p013*p124+p012*p134,  p015*p123-p013*p125+p012*p135,
	         p015*p124-p014*p125+p012*p145,  p015*p134-p014*p135+p013*p145,
	         p015*p234-p014*p235+p013*p245-p012*p345,  p014*p023-p013*p024+p012*p034,
	         p015*p023-p013*p025+p012*p035,  p015*p024-p014*p025+p012*p045,
	         p025*p034-p024*p035+p023*p045,  p024*p123-p023*p124+p012*p234,
	         p025*p123-p023*p125+p012*p235,  p025*p124-p024*p125+p012*p245,
	         p025*p134-p024*p135+p023*p145+p012*p345,  p025*p234-p024*p235+p023*p245,
	         p014*p023-p013*p024+p012*p034,  p015*p023-p013*p025+p012*p035,
	         p015*p034-p014*p035+p013*p045,  p025*p034-p024*p035+p023*p045,
	         p034*p123-p023*p134+p013*p234,  p035*p123-p023*p135+p013*p235,
	         p035*p124-p034*p125-p023*p145+p013*p245,  p035*p134-p034*p135+p013*p345,
	         p035*p234-p034*p235+p023*p345,  p014*p023-p013*p024+p012*p034,
	         p015*p024-p014*p025+p012*p045,  p015*p034-p014*p035+p013*p045,
	         p025*p034-p024*p035+p023*p045,  p034*p124-p024*p134+p014*p234,
	         p045*p123+p034*p125-p024*p135+p014*p235,  p045*p124-p024*p145+p014*p245,
	         p045*p134-p034*p145+p014*p345,  p045*p234-p034*p245+p024*p345,
	         p015*p023-p013*p025+p012*p035,  p015*p024-p014*p025+p012*p045,
	         p015*p034-p014*p035+p013*p045,  p025*p034-p024*p035+p023*p045,
	         p045*p123-p035*p124+p025*p134-p015*p234,  p035*p125-p025*p135+p015*p235,
	         p045*p125-p025*p145+p015*p245,  p045*p135-p035*p145+p015*p345,
	         p045*p235-p035*p245+p025*p345,  p014*p123-p013*p124+p012*p134,
	         p015*p123-p013*p125+p012*p135,  p015*p124-p014*p125+p012*p145,
	         p024*p123-p023*p124+p012*p234,  p025*p123-p023*p125+p012*p235,
	         p025*p124-p024*p125+p012*p245,  p045*p123-p035*p124+p034*p125-p012*p345,
	         p125*p134-p124*p135+p123*p145,  p125*p234-p124*p235+p123*p245,
	         p014*p123-p013*p124+p012*p134,  p015*p123-p013*p125+p012*p135,
	         p015*p134-p014*p135+p013*p145,  p034*p123-p023*p134+p013*p234,
	         p035*p123-p023*p135+p013*p235,  p045*p123+p025*p134-p024*p135+p013*p245,
	         p035*p134-p034*p135+p013*p345,  p125*p134-p124*p135+p123*p145,
	         p135*p234-p134*p235+p123*p345,  p014*p123-p013*p124+p012*p134,
	         p015*p124-p014*p125+p012*p145,  p015*p134-p014*p135+p013*p145,
	         p034*p124-p024*p134+p014*p234,  p035*p124-p025*p134-p023*p145+p014*p235,
	         p045*p124-p024*p145+p014*p245,  p045*p134-p034*p145+p014*p345,
	         p125*p134-p124*p135+p123*p145,  p145*p234-p134*p245+p124*p345,
	         p015*p123-p013*p125+p012*p135,  p015*p124-p014*p125+p012*p145,
	         p015*p134-p014*p135+p013*p145,  p034*p125-p024*p135+p023*p145+p015*p234,
	         p035*p125-p025*p135+p015*p235,  p045*p125-p025*p145+p015*p245,
	         p045*p135-p035*p145+p015*p345,  p125*p134-p124*p135+p123*p145,
	         p145*p235-p135*p245+p125*p345,  p024*p123-p023*p124+p012*p234,
	         p025*p123-p023*p125+p012*p235,  p034*p123-p023*p134+p013*p234,
	         p035*p123-p023*p135+p013*p235,  p045*p123-p023*p145-p015*p234+p014*p235,
	         p025*p234-p024*p235+p023*p245,  p035*p234-p034*p235+p023*p345,
	         p125*p234-p124*p235+p123*p245,  p135*p234-p134*p235+p123*p345,
	         p024*p123-p023*p124+p012*p234,  p025*p124-p024*p125+p012*p245,
	         p034*p124-p024*p134+p014*p234,  p035*p124-p024*p135+p015*p234+p013*p245,
	         p045*p124-p024*p145+p014*p245,  p025*p234-p024*p235+p023*p245,
	         p045*p234-p034*p245+p024*p345,  p125*p234-p124*p235+p123*p245,
	         p145*p234-p134*p245+p124*p345,  p025*p123-p023*p125+p012*p235,
	         p025*p124-p024*p125+p012*p245,  p034*p125-p025*p134+p014*p235-p013*p245,
	         p035*p125-p025*p135+p015*p235,  p045*p125-p025*p145+p015*p245,
	         p025*p234-p024*p235+p023*p245,  p045*p235-p035*p245+p025*p345,
	         p125*p234-p124*p235+p123*p245,  p145*p235-p135*p245+p125*p345,
	         p034*p123-p023*p134+p013*p234,  p034*p124-p024*p134+p014*p234,
	         p034*p125-p025*p134+p015*p234-p012*p345,  p035*p134-p034*p135+p013*p345,
	         p045*p134-p034*p145+p014*p345,  p035*p234-p034*p235+p023*p345,
	         p045*p234-p034*p245+p024*p345,  p135*p234-p134*p235+p123*p345,
	         p145*p234-p134*p245+p124*p345,  p035*p123-p023*p135+p013*p235,
	         p035*p124-p024*p135+p014*p235+p012*p345,  p035*p125-p025*p135+p015*p235,
	         p035*p134-p034*p135+p013*p345,  p045*p135-p035*p145+p015*p345,
	         p035*p234-p034*p235+p023*p345,  p045*p235-p035*p245+p025*p345,
	         p135*p234-p134*p235+p123*p345,  p145*p235-p135*p245+p125*p345,
	         p045*p123-p023*p145+p013*p245-p012*p345,  p045*p124-p024*p145+p014*p245,
	         p045*p125-p025*p145+p015*p245,  p045*p134-p034*p145+p014*p345,
	         p045*p135-p035*p145+p015*p345,  p045*p234-p034*p245+p024*p345,
	         p045*p235-p035*p245+p025*p345,  p145*p234-p134*p245+p124*p345,
	         p145*p235-p135*p245+p125*p345])
Random.seed!(133713371337); TropI,wG = tropical_variety(I,val)
