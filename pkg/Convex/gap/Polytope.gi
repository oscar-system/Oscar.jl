#############################################################################
##
##  Functions.gi        Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
## Chapter: Polytopes
##
#############################################################################


####################################
##
##  Attributes of polyopes
##
####################################

##
InstallMethod( ExternalPolymakePolytope,
               "for polyopes",
               [ IsPolytope ],
    function( poly )
    local old_pointlist, new_pointlist, ineqs, i,j;
    
    if IsBound( poly!.input_points ) and IsBound( poly!.input_ineqs ) then
        Error( "points and inequalities at the same time are not supported\n" );
    fi;
    
    if IsBound( poly!.input_points ) then
        
        old_pointlist := poly!.input_points;
        new_pointlist:= [ ];
        for i in old_pointlist do 
            j:= ShallowCopy( i );
            Add( j, 1, 1 );
            Add( new_pointlist, j );
        od;
        return Polymake_PolytopeByGenerators( new_pointlist );
        
    elif  IsBound( poly!.input_ineqs ) then
        
        ineqs := ShallowCopy( poly!.input_ineqs );
        return Polymake_PolytopeFromInequalities( ineqs );
        
    else
        
        Error( "something went wrong\n" );
        
    fi;
    
end );


InstallMethod( VerticesOfPolytope,
               "for polyopes",
               [ IsPolytope ],
  function( poly )
    
    if PolymakeAvailable() then
        return Polymake_V_Rep( ExternalPolymakePolytope( poly ) )!.generating_vertices;
    fi;
    TryNextMethod();
    
end );

InstallMethod( DefiningInequalities,
               "for polyopes",
               [ IsPolytope ],
  function( poly )
    
    if PolymakeAvailable() then
        return Polymake_H_Rep( ExternalPolymakePolytope( poly ) )!.inequalities;
    fi;
    TryNextMethod();
    
end );


##
InstallMethod( FacetInequalities,
               " for external polyopes",
               [ IsExternalPolytopeRep ],
  function( poly )
    
    if PolymakeAvailable() then
        return Polymake_H_Rep( ExternalPolymakePolytope( poly ) )!.inequalities;
    fi;
    TryNextMethod();
    
end );


InstallMethod( IsBounded,
               " for external polyopes.",
               [ IsPolytope ],
    function( poly )
    
    if PolymakeAvailable() then
        return Polymake_IsBounded( ExternalPolymakePolytope( poly ) );
    fi;
    
    TryNextMethod();
    
end );
