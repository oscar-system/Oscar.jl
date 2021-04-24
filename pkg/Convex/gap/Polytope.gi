#############################################################################
##
##  Functions.gi        Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
##  Chapter: Polytopes
##
#############################################################################


####################################
##
##  Attributes of polyopes
##
####################################

    
InstallMethod( ExternalPolymakePolytope,
               "for polyopes",
               [ IsPolytope ],
    function( poly )
    
    if IsBound( poly!.input_points ) and Length( poly!.input_points ) = 1 and IsZero( poly!.input_points ) then
        return Polymake_PolytopeByGenerators( poly!.input_points[ 1 ] );
    fi;
    
    if IsBound( poly!.input_points ) then
        return Polymake_PolytopeByGenerators( poly!.input_points );
    fi;
    
    if IsBound( poly!.input_equalities ) then
        return Polymake_PolytopeFromInequalities( poly!.input_inequalities, poly!.input_equalities );
    fi;
    
    # otherwise our fallback is poly from inequalities
    return Polymake_PolytopeFromInequalities( poly!.input_inequalities );
    
end );


InstallMethod( VerticesOfPolytope,
               "for polyopes",
               [ IsPolytope ],
  function( poly )
    
    if PolymakeAvailable() then
        return Polymake_V_Rep( ExternalPolymakePolytope( poly ) )!.vertices;
    fi;
    TryNextMethod();
    
end );


InstallMethod( DefiningInequalities,
               "for polyopes",
               [ IsPolytope ],
  function( poly )
    local ineqs, eqs;
    
    if PolymakeAvailable() then
        ineqs := Polymake_Inequalities( ExternalPolymakePolytope( poly ) );
        eqs := Polymake_Equalities( ExternalPolymakePolytope( poly ) );
        return Set( Concatenation( eqs, (-1) * eqs, ineqs ) );
    fi;
    TryNextMethod();
    
end );


InstallMethod( FacetInequalities,
               " for external polyopes",
               [ IsExternalPolytopeRep ],
  function( poly )
    
    if PolymakeAvailable() then
        return Polymake_Inequalities( ExternalPolymakePolytope( poly ) );
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


InstallMethod( LatticePoints,
               [ IsPolytope ],
    function( poly )
    
    if PolymakeAvailable() then
        return Polymake_LatticePoints( ExternalPolymakePolytope( poly ) );
    fi;
    
    TryNextMethod();
    
end );


InstallMethod( Dimension,
               [ IsPolytope ],
    function( poly )
    
    if PolymakeAvailable() then
        return Polymake_Dimension( ExternalPolymakePolytope( poly ) );
    fi;
    
    TryNextMethod();
    
end );
