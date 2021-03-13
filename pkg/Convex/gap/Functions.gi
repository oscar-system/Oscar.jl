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
InstallMethod( PolymakeAvailability, [  ],
  function( )
    local available;
    available := false;
    
    # Check if Polymake and the gap-julia interface are available
    if ( TestPackageAvailability( "JuliaInterface", ">= 0.5.2" ) and IsPackageMarkedForLoading( "JuliaInterface", ">= 0.5.2" ) ) then
        if JuliaImportPackage("Polymake") then
            available := true;
        fi;
    fi;
    
    return available;
    
end );

##
InstallMethod( LoadPolymake, [  ],
  function( )
    local result;
    result := false;
    
    if PolymakeAvailability() then
        ImportJuliaModuleIntoGAP( "Polymake" );
        result := true;
    fi;
    
    return result;
    
end );
