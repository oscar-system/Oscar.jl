#############################################################################
##
##  Functions.gi        Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
## Chapter: Functionality
##
#############################################################################

## Section: Availability and loading of Polymake

##
InstallMethod( FacetInequalities,
               " for external polytopes",
               [ IsExternalPolytopeRep ],
  function( polytope )
    local P, s;
    
    if PolymakeAvailable() then
        
        # parse the vertices of polytope into format for Polymake
        
        # issue commands in Julia
        JuliaEvalString( "F = Julia.Polymake.polytope.Polytope( POINTS = [1 -1 -1; 1 1 -1; 1 -1 1; 1 1 1; 1 0 0] ).FACETS" );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        
        # cast the resulting string into a list of integers
        Error( Concatenation( "Test", s ) );
        return false;
        
    fi;
    
    # otherwise try next method
    TryNextMethod();
    
end );
