#############################################################################
##
##  Cone.gi             Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
## Chapter: Cones
##
#############################################################################

InstallMethod( RayGenerators,
               [ IsCone ],
               
    function( cone )
    local input_rays;
    
    # find the input rays
    input_rays := [ Concatenation( [ 1 ], cone!.input_rays[ 1 ] ) ];
    
    # create polyhedron by generators in polymake
    # finds the ry generators and return them
    
    Error( "Test new method" );
    return false;
    
    #return Cdd_GeneratingRays( ExternalCddCone( cone ) );
    
end );


InstallMethod( ExternalPolymakeCone,
               [ IsCone ],
               
   function( cone )
   
   local list, new_list, number_of_equalities, linearity, i, u ;
   
   new_list:= [ ];
   if IsBound( cone!.input_rays ) and Length( cone!.input_rays )= 1 and IsZero( cone!.input_rays ) then
   
      new_list:= [ Concatenation( [ 1 ], cone!.input_rays[ 1 ] ) ];
      
      return false;
      #return Polymake_PolyhedronByGenerators( new_list );
      
   fi;
   
   if IsBound( cone!.input_rays ) then 
   
      list := cone!.input_rays;
      
      for i in [1..Length( list ) ] do 
          
          u:= ShallowCopy( list[ i ] );
          
          Add( u, 0, 1 );
          
          Add( new_list, u );
      
      od;
      
      return false;
      #return Polymake_PolyhedronByGenerators( new_list );
   
   fi;
   
   
   if IsBound( cone!.input_equalities ) then
   
      list := StructuralCopy( cone!.input_equalities );
      
      number_of_equalities:= Length( list );
      
      linearity := [1..number_of_equalities];
      
      Append( list, StructuralCopy( cone!.input_inequalities ) );
      
      for i in [1..Length( list ) ] do 
      
          u:= ShallowCopy( list[ i ] );
          
          Add( u, 0, 1 );
          
          Add( new_list, u );
      
      od;
      
      return false;
      #return Polymake_PolyhedronByInequalities( new_list, linearity );
   
   else 
   
      list:= StructuralCopy( cone!.input_inequalities );
      
      for i in [1..Length( list ) ] do 
          
          u:= ShallowCopy( list[ i ] );
          
          Add( u, 0, 1 );
          
          Add( new_list, u );
          
      od;
      
      #return Polymake_PolyhedronByInequalities( new_list );
      return false;
   
   fi;
   
end );
