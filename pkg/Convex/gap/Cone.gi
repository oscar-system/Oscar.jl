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


InstallMethod( ExternalPolymakeCone,
               [ IsCone ],
   function( cone )
    local list, new_list, number_of_equalities, linearity, i, u ;
    
    new_list:= [ ];
    if IsBound( cone!.input_rays ) and Length( cone!.input_rays )= 1 and IsZero( cone!.input_rays ) then
        
        new_list:= [ Concatenation( [ 1 ], cone!.input_rays[ 1 ] ) ];
        return Polymake_ConeByGenerators( new_list );
        
    fi;
    
    if IsBound( cone!.input_rays ) then 
        
        list := cone!.input_rays;
        
        for i in [1..Length( list ) ] do 
            u:= ShallowCopy( list[ i ] );
            Add( u, 0, 1 );
            Add( new_list, u );
        od;
        
        return Polymake_ConeByGenerators( new_list );
        
    fi;
    
    
    if IsBound( cone!.input_equalities ) then
        
        list := StructuralCopy( cone!.input_equalities );
        number_of_equalities:= Length( list );
        linearity := [1..number_of_equalities];
        
        Append( list, StructuralCopy( cone!.input_inequalities ) );
        
        for i in [1..Length( list ) ] do 
            u:= ShallowCopy( list[ i ] );
            Add( u, 0, 1 );
            Add( new_list, u );
        od;
        
        return Polymake_ConeFromInequalities( new_list, linearity );
        
    else
        
        list:= StructuralCopy( cone!.input_inequalities );
        
        for i in [1..Length( list ) ] do 
            u:= ShallowCopy( list[ i ] );
            Add( u, 0, 1 );
            Add( new_list, u );
        od;
        
        return Polymake_ConeFromInequalities( new_list );
        
    fi;
    
end );


InstallMethod( RayGenerators,
               [ IsCone ],
    function( cone )
    
    if PolymakeAvailable() then
        return Polymake_GeneratingRays( ExternalPolymakeCone( cone ) );
    fi;
    
    TryNextMethod();
    
end );

InstallMethod( IsPointed,
               [ IsCone ],
    function( cone )
    
    if PolymakeAvailable() then
        return Polymake_IsPointed( ExternalPolymakeCone( cone ) );
    fi;
    
    TryNextMethod();
    
end );

InstallMethod( Dimension,
               [ IsCone ],
    function( cone )
    
    if PolymakeAvailable() then
        return Polymake_Dimension( ExternalPolymakeCone( cone ) );
    fi;
    
    TryNextMethod();
    
end );

InstallMethod( RaysInFacets,
               [ IsCone ],
    function( cone )
    
    if PolymakeAvailable() then
        return Polymake_RaysInFacets( ExternalPolymakeCone( cone ) );
    fi;
    
    TryNextMethod();
    
end );


InstallMethod( DefiningInequalities,
               [ IsCone ],
  function( cone )
    local ineqs, eqs;
    
    if PolymakeAvailable() then
        ineqs := Polymake_Inequalities( ExternalPolymakeCone( cone ) );
        eqs := Polymake_Equalities( ExternalPolymakeCone( cone ) );
        return Set( Concatenation( equ, (-1) * equ, ineq ) );
    fi;
    
    TryNextMethod();
    
end );
