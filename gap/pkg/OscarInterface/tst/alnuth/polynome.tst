# this needs RadiRoot to be loaded
gap> START_TEST("Factorisation of polynomials using PARI/GP");  

# some splitting fields
gap> pol := UnivariatePolynomial( Rationals, [ -1, 0, 0, 0, 1, 0, 1 ] );;
gap> SplittingField( pol );
<algebraic extension over the Rationals of degree 24>
gap> pol := UnivariatePolynomial( Rationals, [ 4, -1, 0, 6, -1, -2, 2, 1 ] );;
gap> SplittingField( pol );
<algebraic extension over the Rationals of degree 42>
gap> pol := UnivariatePolynomial( Rationals, [ 4, 0, 4, 0, -3, 0, -1, 0, 1 ] );;
gap> SplittingField( pol );
<algebraic extension over the Rationals of degree 24>
gap> pol := UnivariatePolynomial(Rationals, 
>                                [ -3, 1, -1, 2, 1, 3, -2, -1, 0, 1 ]);;
gap> SplittingField( pol );
<algebraic extension over the Rationals of degree 54>
gap> pol := UnivariatePolynomial(Rationals, 
>                                [ 47, 0, 103, 0, 41, 0, 7, 0, -2, 0, 1 ]);;
gap> SplittingField( pol );
<algebraic extension over the Rationals of degree 10>

# example from email to Bill Allombert (10/05/11)
gap> pol := UnivariatePolynomial( Rationals, [ 1, 2, 2, 2, 2, 1, 1 ] );
x_1^6+x_1^5+2*x_1^4+2*x_1^3+2*x_1^2+2*x_1+1
gap> f := UnivariatePolynomial( Rationals, [ 1, 1, 1 ] );
x_1^2+x_1+1
gap> K := SplittingField( pol );
<algebraic extension over the Rationals of degree 48>
gap> FactorsPolynomialAlgExt( K, f );
[ x_1^2+x_1+!1 ]

# example from Bill Allombert (10/05/11)
gap> pol := UnivariatePolynomial( Rationals, [ 20736, 0, 4147200, 0, 60632064, 
> 0, 347286528, 0, 1078555392, 0, 2069549568, 0, 2613917952, 0, 2241209088, 0,
> 1318874976, 0, 530669952, 0, 143684928, 0, 25510464, 0, 2872752, 0, 195936, 0,
> 7536, 0, 144, 0, 1 ] );;
gap> L := FieldByPolynomial( pol );
<algebraic extension over the Rationals of degree 32>
gap> facs := FactorsPolynomialAlgExt( L, pol );;
gap> Collected(List(facs, Degree));
[ [ 1, 32 ] ]
gap> Product(facs) = AlgExtEmbeddedPol(L, pol);
true

#
gap> STOP_TEST( "polynome.tst", 10000000);   
