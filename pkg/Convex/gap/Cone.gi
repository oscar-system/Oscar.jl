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
        
        new_list:= [ Concatenation( cone!.input_rays[ 1 ] ) ];
        return Polymake_ConeByGenerators( new_list );
        
    fi;
    
    if IsBound( cone!.input_rays ) then 
        return Polymake_ConeByGenerators( new_list );
    fi;
    
    if IsBound( cone!.input_equalities ) then
        return Polymake_ConeFromInequalities( cone!.input_inequalities, cone!.input_equalities );
    else
        return Polymake_ConeFromInequalities( cone!.input_inequalities );
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
        return Set( Concatenation( eqs, (-1) * eqs, ineqs ) );
    fi;
    
    TryNextMethod();
    
end );


InstallMethod( IntersectionOfCones,
               "for homalg cones",
               [ IsCone, IsCone ],
  function( cone1, cone2 )
    local ext_cone, ineqs, eqs, cone;
    
    if not Rank( ContainingGrid( cone1 ) ) = Rank( ContainingGrid( cone2 ) ) then
        Error( "cones are not from the same grid" );
    fi;
    
    # compute the intersection cone
    ext_cone := Polymake_Intersection( ExternalPolymakeCone( cone1), ExternalPolymakeCone( cone2 ) );
    ineqs := ext_cone!.inequalities;
    eqs := ext_cone!.equalities;
    cone := ConeByEqualitiesAndInequalities( ineqs, eqs );
    SetContainingGrid( cone, ContainingGrid( cone1 ) );
    
    # and return it
    return cone;
    
end );
