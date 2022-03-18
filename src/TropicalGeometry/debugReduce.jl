R,(t,x,y,z) = Singular.PolynomialRing(QQ,["t","x","y","z","w"],ordering=Singular.ordering_a([-1,1,1,1])*Singular.ordering_lp())
I = Singular.Ideal(R,[x,y,z])
