#############################################################################
##
##  Cone.gd             JConvex package
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
