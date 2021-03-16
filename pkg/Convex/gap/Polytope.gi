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

InstallMethod( VerticesOfPolytope,
               "for polytopes",
               [ IsPolytope ],
               
  function( polytope )
    
    #return Polymake_GeneratingVertices( ExternalPolymakePolytope( polyt ) );
    return polytope!.input_points;
    
end );


##
InstallMethod( FacetInequalities,
               " for external polytopes",
               [ IsExternalPolytopeRep ],
  function( polytope )
    local v, P, s, string_list, l;
    
    if PolymakeAvailable() then
        
        # parse the vertices of polytope into format for Polymake
        v := VerticesOfPolytope( polytope );
        
        # add a 1 as first argument to all these vertices as polymake requires homogeneous input
        v := List( [ 1 .. Length( v ) ], i -> Concatenation( [ 1 ], v[ i ] ) );
        s := String( v );
        Error( Concatenation( "Test: ", s ) );
        
        # issue commands in Julia
        JuliaEvalString( "F = Julia.Polymake.polytope.Polytope( POINTS = [1 -1 -1; 1 1 -1; 1 -1 1; 1 1 1; 1 0 0] ).FACETS" );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        
        # cast the resulting string into a list of integers
        string_list := SplitString( s, '\n' );
        string_list := List( [ 2 .. Length( string_list ) ], i -> Concatenation( "[", ReplacedString( string_list[ i ], " ", "," ), "]" ) );
        l := EvalString( Concatenation( "[", JoinStringsWithSeparator( string_list, "," ), "]" ) );
        
        Error( Concatenation( "Test", String( l ) ) );
        
        # return this list of inequalities
        return l;
        
    fi;
    
    # otherwise try next method
    TryNextMethod();
    
end );


####################################
##
## Attributes
##
####################################

##
InstallMethod( ExternalPolymakePolytope,
               "for polytopes",
               [ IsPolytope ],
   function( polyt )
   local old_pointlist, new_pointlist, ineqs, i,j;
   
   if IsBound( polyt!.input_points ) and IsBound( polyt!.input_ineqs ) then
        
        Error( "points and inequalities at the same time are not supported\n" );
        
   fi;
    
   if IsBound( polyt!.input_points ) then 
        
        old_pointlist := polyt!.input_points;
        
        new_pointlist:= [ ];
        
        for i in old_pointlist do 
            
            j:= ShallowCopy( i );
            
            Add( j, 1, 1 );
            
            Add( new_pointlist, j );
            
        od;
        
        #return Polymake_PolyhedronByGenerators( new_pointlist );
        return false;
        
    elif  IsBound( polyt!.input_ineqs ) then
        
        ineqs := ShallowCopy( polyt!.input_ineqs );
        
        #return Polymake_PolyhedronByInequalities( ineqs );
        return false;
        
    else
        
        Error( "something went wrong\n" );
        
   fi;
   
end );
