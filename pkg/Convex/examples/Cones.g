#! @Chapter Cones

#! @Section Examples

LoadPackage( "Convex" );

#! The following demonstrates cone operations performed with Polymake:

#! @Example
P:= Cone( [ [ 2, 7 ], [ 0, 12 ], [ -2, 5 ] ] );
#! <A cone in |R^2>
d:= DefiningInequalities( P );
#! [ [ -7, 2 ], [ 5, 2 ] ]
Q:= ConeByInequalities( d );
#! <A cone in |R^2>
DefiningInequalities( Q );
#! [ [ -7, 2 ], [ 5, 2 ] ]
RayGenerators( P );
#! [ [ -2, 5 ], [ 2, 7 ] ]
RayGenerators( Q );
#! [ [ -2, 5 ], [ 2, 7 ] ]
IsPointed( P );
#! true
IsPointed( Q );
#! true
P=Q;
#! true
HilbertBasis( P );
#! [ [ -2, 5 ], [ -1, 3 ], [ 0, 1 ], [ 1, 4 ], [ 2, 7 ] ]
HilbertBasis( Q );
#! [ [ -2, 5 ], [ -1, 3 ], [ 0, 1 ], [ 1, 4 ], [ 2, 7 ] ]
P_dual:= DualCone( P );
#! <A cone in |R^2>
RayGenerators( P_dual );
#! [ [ -7, 2 ], [ 5, 2 ] ]
Dimension( P );
#! 2
List( Facets( P ), RayGenerators );
#! [ [ [ 2, 7 ] ], [ [ -2, 5 ] ] ]
IsRegularCone( P );
#! false
IsRay( P );
#! false
R:= Cone( [ [ 4, 5 ], [ -2, 1 ] ] );
#! <A cone in |R^2>
T:= IntersectionOfCones( P, R );
#! <A cone in |R^2>
RayGenerators( T );
#! [ [ -2, 5 ], [ 2, 7 ] ]
W:= Cone( [ [-3,-4 ] ] );
#! <A ray in |R^2>
I:= IntersectionOfCones( P, W );
#! <A cone in |R^2>
RayGenerators( I );
#! [  ]
Contains( P, I );
#! true
Contains( W, I );
#! true
Contains( P, R );
#! false
Contains( R, P );
#! true
LinealitySpaceGenerators( P );
#! [  ]
P:= Cone( [ [ 1, 1, -3 ], [ -1, -1, 3 ], [ 1, 2, 1 ], [ 2, 1, 2 ] ] );
#! < A cone in |R^3>
IsPointed( P );
#! false
Dimension( P );
#! 3
IsRegularCone( P );
#! false
P;
#! < A cone in |R^3 of dimension 3 with 4 ray generators>
RayGenerators( P );
#! [ [ -1, -1, 3 ], [ 1, 1, -3 ], [ 1, 2, 1 ], [ 2, 1, 2 ] ]
d:= DefiningInequalities( P );
#! [ [ -5, 8, 1 ], [ 7, -4, 1 ] ]
facets:= Facets( P );
#! [ <A cone in |R^3>, <A cone in |R^3> ]
DualCone( P );
#! < A cone in |R^3>
RayGenerators( DualCone( P ) );
#! [ [ 0, -1, 2 ], [ 0, 2, -1 ] ]
LinealitySpaceGenerators( DualCone( P ) );
#! [ [ -1, 0, 1 ] ]
#! @EndExample
