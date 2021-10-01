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
    # Check if Polymake and the gap-julia interface are available
    if IsPackageMarkedForLoading( "JuliaInterface", ">= 0.5.2" ) then
        if IsBoundGlobal("_Polymake_jl") and IsJuliaObject(ValueGlobal("_Polymake_jl")) then
            return true;
        fi;
        if JuliaEvalString("try import Polymake ; return true\ncatch\n return false end") then
            BindGlobal("_Polymake_jl", JuliaPointer(Julia.Polymake));
            return true;
        fi;
    fi;
    return false;
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
