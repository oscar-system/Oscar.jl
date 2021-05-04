#############################################################################
##
##  ExternalPolymakeCone.gd      JConvex package
##                               Martin Bies
##
##  Copyright 2021               University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake
##
#! @Chapter External PolymakeCones
##
#############################################################################


##############################################################################################
##
#! @Section GAP category of PolymakeCones
##
##############################################################################################

#! @Description
#! The GAP category for cones residing in Polymake
#! @Returns true or false
#! @Arguments object
DeclareCategory( "IsPolymakeCone",
                 IsObject );


##############################################################################################
##
#! @Section Constructors for PolymakeCones
##
##############################################################################################

#! @Arguments genes[, linearities_list ]
#! @Returns a PolymakeCone
#! @Description
#! The function takes a list in which every entry represents a vertex in the ambient vector space.
#! In case we want some vertices to be free (the vertex and its negative belong to the PolymakeCone) we should refer
#! in a second list to their indices.
DeclareGlobalFunction( "Polymake_ConeByGenerators" );

#! @Arguments ineq [, linearities_list ]
#! @Returns a PolymakeCone
#! @Description
#! The function takes a list in which every entry represents an inequality (or equality).
#! In case we want some entries to represent equalities we should refer in a second list to their indices.
DeclareGlobalFunction( "Polymake_ConeFromInequalities" );

##############################################################################################
##
#! @Section Attributes of PolymakeCones
##
##############################################################################################

#! @Arguments P
#! @Returns a PolymakeCone
#! @Description
#! The function takes a PolymakeCone and returns its canonical V-rep.
DeclareAttribute( "Polymake_CanonicalConeByGenerators",  IsPolymakeCone  );

#! @Arguments P
#! @Returns a PolymakeCone
#! @Description
#! The function takes a PolymakeCone and returns its canonical H-rep.
DeclareAttribute( "Polymake_CanonicalConeFromInequalities",  IsPolymakeCone  );


#! @Arguments P
#! @Returns a PolymakeCone
#! @Description
#! The function takes a PolymakeCone and returns its reduced $V$-representation.
DeclareAttribute( "Polymake_V_Rep",  IsPolymakeCone  );

#! @Arguments P
#! @Returns a PolymakeCone
#! @Description 
#! The function takes a PolymakeCone and returns its reduced $H$-representation. 
DeclareAttribute( "Polymake_H_Rep",  IsPolymakeCone  );
#! @InsertChunk Example4

#! @Arguments P
#! @Returns The dimension of the ambient space of the PolymakeCone(i.e., the space that contains $P$).
DeclareAttribute( "Polymake_AmbientSpaceDimension", IsPolymakeCone );

#! @Arguments P
#! @Returns The dimension of the PolymakeCone, where the dimension, $\mathrm{dim}(P)$, of a PolymakeCone $P$
#! is the maximum number of affinely independent points in $P$ minus 1.
DeclareAttribute( "Polymake_Dimension", IsPolymakeCone );

#! @Arguments P
#! @Returns list
#! @Description
#! The output is the set of rays of the PolymakeCone.
DeclareAttribute( "Polymake_Rays", IsPolymakeCone );

#! @Arguments P
#! @Returns list
#! @Description
#! The output is the lineality of the PolymakeCone.
DeclareAttribute( "Polymake_Lineality", IsPolymakeCone );

#! @Arguments P
#! @Returns a list
#! @Description
#! The output is the set of equalities of the PolymakeCone.
DeclareAttribute( "Polymake_Equalities", IsPolymakeCone );

#! @Arguments P
#! @Description
#! The output is the set of inequalities of the PolymakeCone.
DeclareAttribute( "Polymake_Inequalities", IsPolymakeCone );

#! @Arguments P
#! @Description
#! The output is the incident matrix of the ray generators in the facets.
DeclareAttribute( "Polymake_RaysInFacets", IsPolymakeCone );

#! @Arguments P
#! @Description
#! The output is the incident matrix of the rays in the faces.
DeclareAttribute( "Polymake_RaysInFaces", IsPolymakeCone );


##############################################################################################
##
#! @Section Properties of PolymakeCones
##
##############################################################################################

#! @Arguments P
#! @Returns true or false
#! @Description
#! The output is <C>true</C> if the PolymakeCone is empty and <C>false</C> otherwise
DeclareProperty( "Polymake_IsEmpty", IsPolymakeCone );

#! @Arguments P
#! @Returns true or false
#! @Description
#! The output is <C>true</C> if the PolymakeCone is pointed and <C>false</C> otherwise
DeclareProperty( "Polymake_IsPointed", IsPolymakeCone );
#! @InsertChunk demo


##############################################################################################
##
#! @Section Operations with cones
##
##############################################################################################

#! @Arguments C
#! @Returns a string
#! @Description
#! Returns the Polymake cone that defines the intersection of the two cones.
DeclareOperation( "Polymake_Intersection", [ IsPolymakeCone, IsPolymakeCone ] );


##############################################################################################
##
#! @Section Command strings
##
##############################################################################################

#! @Arguments C
#! @Returns a string
#! @Description
#! Returns the string which, when executed in Julia, constructs the cone in question as H-representation.
DeclareOperation( "Polymake_H_Rep_command_string", [ IsPolymakeCone ] );

#! @Arguments C
#! @Returns a string
#! @Description
#! Returns the string which, when executed in Julia, constructs the cone in question as V-representation.
DeclareOperation( "Polymake_V_Rep_command_string", [ IsPolymakeCone ] );
