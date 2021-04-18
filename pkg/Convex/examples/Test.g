#! @Chapter Functionality

#! @Section Examples

LoadPackage( "Convex" );

#! Eventually, we will perform non-trivial tests and exemplify how this package works. This is a first (trivial) test when setting up the package.

#! @Example
1+1;
#! 2
#! @EndExample

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
#! [ [ [ -2, 5 ] ], [ [ 2, 7 ] ] ]
faces := FacesOfCone( P );
#! [ <A cone in |R^2>, <A cone in |R^2>, <A ray in |R^2>, 
#!  <A ray in |R^2> ]
RelativeInteriorRay( P );
#! [ -2, 41 ]
IsRelativeInteriorRay( [ -2, 41 ], P );
#! true
IsRelativeInteriorRay( [ 2, 7 ], P );
#! false
LinealitySpaceGenerators( P );
#! [  ]
IsRegularCone( P );
#! false
IsRay( P );
#! false
proj_x1:= FourierProjection( P, 2 );
#! <A cone in |R^1>
RayGenerators( proj_x1 );
#! [ [ -1 ], [ 1 ] ]
DefiningInequalities( proj_x1 );
#! [ [ 0 ] ]
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
cdd_cone:= ExternalCddCone( P );
#! < Polyhedron given by its V-representation >
Display( cdd_cone );
#! V-representation 
#! begin 
#! 3 X 3  rational
#!                
#!    0   2   7 
#!    0   0  12 
#!    0  -2   5 
#! end
Cdd_Dimension( cdd_cone );
#! 2
H:= Cdd_H_Rep( cdd_cone );
#! < Polyhedron given by its H-representation >
Display( H );
#! H-representation 
#! begin 
#!    2 X 3  rational
#!                
#!    0   5   2 
#!    0  -7   2 
#! end
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
faces := FacesOfCone( P );
#! [ <A cone in |R^3>, <A cone in |R^3>, <A cone in |R^3>, 
#!  <A cone in |R^3>, <A cone in |R^3> ]
FVector( P );
#! [ 1, 2, 1 ]
List( faces, Dimension );
#! [ 0, 3, 2, 1, 2 ]
L_using_4ti2 := [ [ [ 0, 0, 0 ] ], [ [ -2, -1, 10 ], 
[ 0, 0, 1 ], [ 2, 1, 2 ] ],  [ [ 1, 1, -3 ] ] ];;
L_using_Normaliz := [ [ [ 0, 0, 0 ] ], [ [ -1, 0, 7 ], 
[ 0, 0, 1 ], [ 1, 0, 5 ] ], [ [ 1, 1, -3 ] ] ];;
L := LatticePointsGenerators( P );;
L = L_using_4ti2 or L = L_using_Normaliz;
#! true
DualCone( P );
#! < A cone in |R^3>
RayGenerators( DualCone( P ) );
#! [ [ -5, 8, 1 ], [ 7, -4, 1 ] ]
Q_x1x3:= FourierProjection(P, 2 );
#! <A cone in |R^2>
RayGenerators( Q_x1x3 );
#! [ [ -1, 3 ], [ 1, -3 ], [ 1, 1 ] ]
#! @EndExample
