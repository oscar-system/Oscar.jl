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

## Section: Availability of Polymake

##
InstallMethod( PolymakeAvailability, [  ],
  function( )
    local available;
    
    # initialize available
    available := true;
    
    # Check if TopcomInterface is available
    if TestPackageAvailability( "JuliaInterface", ">= 0.5.2" ) = fail then
      Print( "The Gap package JuliaInterface is not available\n" );
      available := false;
    fi;
    
    # and marked to be loaded as suggested package
    if not IsPackageMarkedForLoading( "JuliaInterface", ">= 0.5.2" ) then
      Print( "The Gap package JuliaInterface has not been marked to be loaded by Convex if available.\n" );
      available := false;
    fi;
    
    # finally test if polymake is available in Julia
    #ImportJuliaModuleInGAP("polymake")
    
    # return the result of the operation
    return available;
    
end );
