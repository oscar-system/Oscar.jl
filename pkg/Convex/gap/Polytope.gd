#############################################################################
##
##  Polytope.gd         Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake
##
#! @Chapter Polytopes
##
#############################################################################


####################################
##
#! @Section Attributes of polytopes
##
####################################

#! @Arguments P
#! @Returns a PolymakePolyhedron
#! @Description
#! Converts the polytope to a PolymakePolytope. The operations of Polymake can then be applied.
DeclareAttribute( "ExternalPolymakePolytope", IsPolytope );
