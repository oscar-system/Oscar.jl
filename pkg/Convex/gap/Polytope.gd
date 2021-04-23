#############################################################################
##
##  Polytope.gd         Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
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

#! @Arguments P
#! @Returns a list of lists
#! @Description
#! The operation returns the vertices of the polytope
DeclareAttribute( "VerticesOfPolytope", IsPolytope );

#! @Arguments P
#! @Returns a list of lists
#! @Description
#! The operation returns the list of the inequalities of the facets.
#! Each defining inequality that is not defining-equality of the 
#! polytope is a facet inequality.
DeclareAttribute( "FacetInequalities", IsPolytope );
