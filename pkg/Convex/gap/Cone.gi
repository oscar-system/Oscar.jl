#############################################################################
##
##  Cone.gi             Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
## Chapter: Cones
##
#############################################################################

InstallMethod( RayGenerators,
               [ IsCone ],
               
    function( cone )
    Error( "Test new method" );
    return false;
    
    #return Cdd_GeneratingRays( ExternalCddCone( cone ) );
    
end );
