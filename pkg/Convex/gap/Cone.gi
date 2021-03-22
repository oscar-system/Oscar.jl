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
    local input_rays, string_list, command_string, s, l;
    
    # compute the ray generators in Polymake
    if PolymakeAvailable() then
        
        # Parse the rays into format recognized by Polymake
        input_rays := cone!.input_rays;
        string_list := List( [ 1 .. Length( input_rays ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( input_rays[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ", JoinStringsWithSeparator( string_list, "; " ), " ] ).RAYS" );
        
        # issue command in Julia and fetch result as string
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        
        # cast the result into a list of integers
        string_list := SplitString( s, '\n' );
        string_list := List( [ 2 .. Length( string_list ) ], i -> Concatenation( "[", ReplacedString( string_list[ i ], " ", "," ), "]" ) );
        l := EvalString( Concatenation( "[", JoinStringsWithSeparator( string_list, "," ), "]" ) );
        
        #return Cdd_GeneratingRays( ExternalCddCone( cone ) );
        #Error( "Test new method" );
        return l;
        
    fi;
    
    # otherwise try next method
    TryNextMethod();
    
end );

InstallMethod( IsPointed,
               [ IsCone ],
               
    function( cone )
    local input_rays, string_list, command_string, s, l;
    
    # compute the ray generators in Polymake
    if PolymakeAvailable() then
        
        # Parse the rays into format recognized by Polymake
        input_rays := cone!.input_rays;
        string_list := List( [ 1 .. Length( input_rays ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( input_rays[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ", JoinStringsWithSeparator( string_list, "; " ), " ] ).POINTED" );
        
        # issue command in Julia and fetch result as string
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        
        # return the result
        return EvalString( s );
        
    fi;
    
    # otherwise try next method
    TryNextMethod();
    
end );

InstallMethod( Dimension,
               [ IsCone ],
    function( cone )
    local input_rays, string_list, command_string, s, l;
    
    # compute the ray generators in Polymake
    if PolymakeAvailable() then
        
        # Parse the rays into format recognized by Polymake
        input_rays := cone!.input_rays;
        string_list := List( [ 1 .. Length( input_rays ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( input_rays[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ", JoinStringsWithSeparator( string_list, "; " ), " ] ).CONE_DIM" );
        
        # issue command in Julia and fetch result as string
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        
        # return the result
        return EvalString( s );
        
    fi;
    
    # otherwise try next method
    TryNextMethod();
    
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
