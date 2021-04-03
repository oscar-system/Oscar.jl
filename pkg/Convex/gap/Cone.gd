#############################################################################
##
##  Cone.gd             Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
#! @Chapter Cones
##
#############################################################################

###############################
##
#! @Section Attributes of Cones
##
###############################

#! @Arguments C
#! @Returns a polymake cone
#! @Description
#! Converts the cone to a Polymake cone.
DeclareAttribute( "ExternalPolymakeCone",  IsCone  );

#! @Arguments C
#! @Returns a list of all rays contained in the facets of the cone
#! @Description
#! Converts the cone to a Polymake cone.
DeclareAttribute( "RaysInFacets", IsCone );

#! @Arguments C
#! @Returns a list of the defining inequalities of the given cone.
#! @Description
#! Converts the cone to a Polymake cone.
DeclareAttribute( "DefiningInequalities", IsCone );

