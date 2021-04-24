#############################################################################
##
##  Cone.gd             Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake
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
#! Return the rays in the facets of the cone.
#DeclareAttribute( "RaysInFacets", IsCone );

#! @Arguments C
#! @Returns a list of all rays contained in the facets of the cone
#! @Description
#! Returns the rays in all faces of the cone.
#DeclareAttribute( "RaysInFaces", IsCone );

#! @Arguments C
#! @Returns a list of the defining inequalities of the given cone.
#! @Description
#! Return the defining inequalities of the cone.
#DeclareAttribute( "DefiningInequalities", IsCone );

#! @Arguments C1, C2
#! @Returns a cone
#! @Description
#! Returns the intersection of two cones.
#DeclareOperation( "IntersectionOfCones",
#                  [ IsCone, IsCone ] );
