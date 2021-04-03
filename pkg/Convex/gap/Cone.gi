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
        
        return Polymake_ConeByInequalities( new_list, linearity );
        
    else
        
        list:= StructuralCopy( cone!.input_inequalities );
        
        for i in [1..Length( list ) ] do 
            u:= ShallowCopy( list[ i ] );
            Add( u, 0, 1 );
            Add( new_list, u );
        od;
        
        return Polymake_ConeByInequalities( new_list );
        
    fi;
    
end );


InstallMethod( RayGenerators,
               [ IsCone ],
    function( cone )
    
    if PolymakeAvailable() then
        return Set( Polymake_V_Rep( ExternalPolymakeCone( cone ) )!.generating_rays );
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
    local input_rays, string_list, command_string, s;
    
    # compute the ray generators in Polymake
    if PolymakeAvailable() then
        return Polymake_Dimension( ExternalPolymakeCone( cone ) );
    fi;
    
    TryNextMethod();
    
end );

InstallMethod( RaysInFacets,
               " for cones",
               [ IsCone ],
    
    function( cone )
    local ineqs, input_rays, ray_list, i, ray_list_for_facet, j, product, k;
    
    # compute the ray generators in Polymake
    if PolymakeAvailable() then
        
        # compute the inequalities which each define a facet of the cone in question
        ineqs := DefiningInequalities( cone );
        input_rays := RayGenerators( cone );
        
        # now compute the incident matrix
        ray_list := [];
        for i in [ 1 .. Length( ineqs ) ] do
            ray_list_for_facet := [];
            for j in [ 1 .. Length( input_rays ) ] do
                product := List( [ 1 .. Length( ineqs[ i ] ) ], k -> ineqs[ i ][ k ] * input_rays[ j ][ k ] );
                if ( product > 0 ) then
                    Append( ray_list_for_facet, [ 1 ] );
                else
                    Append( ray_list_for_facet, [ 0 ] );
                fi;
            od;
            Append( ray_list, [ ray_list_for_facet ] );
        od;
        
        # return the result
        return ray_list;
        
    fi;
    
    # otherwise try next method
    TryNextMethod();
    
end );


InstallMethod( DefiningInequalities,
               [ IsCone ],
  function( cone )
    
    # compute the ray generators in Polymake
    if PolymakeAvailable() then
        return Polymake_Inequalities( External_PolymakeCone( cone ) );
    fi;
    
    TryNextMethod();
    
end );
