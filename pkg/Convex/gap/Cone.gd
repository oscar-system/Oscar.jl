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

##############################
##
##  Attributes
##
##############################

#! @Section Attributes of Cones

# DeclareAttribute( "RayGenerators",
#                    IsCone );


#! @Arguments C
#! @Returns a cdd object
#! @Description
#! Converts the cone to a polymake polyhedron. The operations of Polymake can then be applied.
DeclareAttribute( "ExternalPolymakeCone",  IsCone  );

DeclareAttribute( "RaysInFacets",
                  IsCone );

#! @Arguments C
#! @Returns a list
#! @Description
#! Returns the list of the defining inequalities of the cone <C>C</C>.
DeclareAttribute( "DefiningInequalities", 
                   IsCone );
