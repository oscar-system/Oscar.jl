#############################################################################
##
##  ExternalPolymakeCone.gd      Convex package
##                               Martin Bies
##
##  Copyright 2021               University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
#! @Chapter Cones in Polymake
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
#! in a second list to their indices .
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
#! @Returns The reduced generating vertices of the PolymakeCone.
DeclareAttribute( "Polymake_GeneratingVertices", IsPolymakeCone );

#! @Arguments P
#! @Returns list
#! @Description
#! The output is the reduced generating rays of the PolymakeCone.
DeclareAttribute( "Polymake_GeneratingRays", IsPolymakeCone );

#! @Arguments P
#! @Returns a list
#! @Description
#! The output is the reduced equalities of the PolymakeCone.
DeclareAttribute( "Polymake_Equalities", IsPolymakeCone );

#! @Arguments P
#! @Description
#! The output is the reduced inequalities of the PolymakeCone.
DeclareAttribute( "Polymake_Inequalities", IsPolymakeCone );


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
#! @Section Command strings
##
##############################################################################################

#! @Arguments C
#! @Returns a string
#! @Description
#! Returns the string which, when executed in Julia, constructs the polytope in question as H-representation.
DeclareProperty( "Polymake_H_Rep_command_string", IsPolymakeCone );

#! @Arguments C
#! @Returns a string
#! @Description
#! Returns the string which, when executed in Julia, constructs the polytope in question as V-representation.
DeclareProperty( "Polymake_V_Rep_command_string", IsPolymakeCone );
