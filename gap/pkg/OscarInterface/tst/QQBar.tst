#@local F, x, z, y, r, M
gap> START_TEST( "QQBar.tst" );

#
gap> F:= QQBarField;
QQBarField
gap> x:= One( F );
<{a1: 1.00000}>
gap> Print( x, "\n" );
<{a1: 1.00000}>
gap> String( x );
"<{a1: 1.00000}>"
gap> z:= Zero( F );
<{a1: 0}>
gap> Zero( x ) = z;
true
gap> z < x;
true
gap> x < z;
false
gap> x < x;
false
gap> z < 1;
true
gap> 0 < x;
true
gap> x = QQBarFieldElement( 1 );
true
gap> y:= x + x;;
gap> 2 * x = y;
true
gap> x * 2 = y;
true
gap> x - x = z;
true
gap> - x = z - x;
true
gap> x * x = x;
true
gap> r:= Sqrt( y );
<{a2: 1.41421}>
gap> r^2 = y;
true
gap> r^-1 = 1 / r;
true
gap> r^-1 = x / r;
true
gap> x = x / 1;
true
gap> AbsoluteValue( -r ) = r;
true
gap> r = QQBarFieldElement( Sqrt( 2 ) );
true

#
gap> Sqrt(-x) = QQBarFieldElement(E(4));
true
gap> Sqrt(-x) = x*E(4);
true

#
gap> M:= [ [ z, r ], [ r, z ] ];;
gap> Determinant( M ) = -2;
true
gap> MinimalPolynomial( M );
x_1^2+(<{a1: -2.00000}>)
gap> Eigenvalues( F, M );
[ <{a2: 1.41421}>, <{a2: -1.41421}> ]
gap> Display( M );
2x2 matrix over QQBarField:
               . <{a2: 1.41421}>
 <{a2: 1.41421}>               .
gap> Display( M : short );
2x2 matrix over QQBarField:
               . <{a2: 1.41421}>
 <{a2: 1.41421}>               .

#
gap> r:= Sqrt( QQBarFieldElement( 2 ) );;
gap> IsRealQQBarFieldElement( r );
true
gap> ComplexConjugate( r ) = r;
true
gap> r:= Sqrt( QQBarFieldElement( -2 ) );;
gap> IsRealQQBarFieldElement( r );
false
gap> ComplexConjugate( r ) = -r;
true

#
gap> Random( F ) in F;
true

#
gap> STOP_TEST( "QQBar.tst" );
