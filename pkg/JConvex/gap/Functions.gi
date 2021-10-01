#############################################################################
##
##  Functions.gi        JConvex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake
##
## Chapter: Functionality
##
#############################################################################

#############################################################################
## Availability of required software
#############################################################################

##
InstallMethod( PolymakeAvailable, [  ],
  function( )
    local available;
    available := false;

    # Check if Polymake and the gap-julia interface are available
    if ( TestPackageAvailability( "JuliaInterface", ">= 0.5.2" ) <> fail ) then
        if IsPackageMarkedForLoading( "JuliaInterface", ">= 0.5.2" ) then
            if JuliaImportPackage("Polymake") then
                available := true;
                ImportJuliaModuleIntoGAP( "Polymake" );
            fi;
        fi;
    fi;

    return available;

end );

##
InstallMethod( CddInterfaceAvailable, [  ],
  function( )
    return IsPackageMarkedForLoading( "CddInterface", ">= 2020.06.24" );
end );

##
InstallMethod( NormalizInterfaceAvailable, [  ],
  function( )
    return IsPackageMarkedForLoading( "NormalizInterface", ">= 1.2.0" );
end );
