#############################################################################
##
##  Functions.gd        JConvex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake
##
#! @Chapter Functionality
##
#############################################################################

#############################################################################
#! @Section Availability of required software
#############################################################################

#! @Arguments
#! @Returns a boolean
#! @Description
#! Checks if the polymake functionality is available in Julia.
DeclareOperation( "PolymakeAvailable", [ ] );

#! @Arguments
#! @Returns a boolean
#! @Description
#! Checks if the CddInterface is available in Julia.
DeclareOperation( "CddInterfaceAvailable", [ ] );

#! @Arguments
#! @Returns a boolean
#! @Description
#! Checks if the NormalizInterface is available in Julia.
DeclareOperation( "NormalizInterfaceAvailable", [ ] );
