#############################################################################
##
##  ExternalPolymakeCone.gd      Convex package
##                               Martin Bies
##
##  Copyright 2021               University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
##  Chapter Cones in Polymake
##
#############################################################################


##############################################################################################
##
##  Section GAP category of PolymakeCones
##
##############################################################################################

DeclareRepresentation( "IsPolymakeConeRep", IsPolymakeCone and IsAttributeStoringRep, [ ] );

BindGlobal( "TheFamilyOfPolymakeCones", NewFamily( "TheFamilyOfPolymakeCones" ) );

BindGlobal( "TheTypeOfPolymakeCone", NewType( TheFamilyOfPolymakeCones, IsPolymakeConeRep ) );


##############################################################################################
##
##  Constructors for PolymakeCones
##
##############################################################################################

InstallGlobalFunction( Polymake_ConeByGenerators,
  function( arg )
    local cone, i, matrix, temp, dim;
    
    if Length( arg )= 0 then
        
        Error( "Wronge input" );
        
    elif Length( arg ) = 1 and IsList( arg[1] ) then
        
        return Polymake_ConeByGenerators( arg[ 1 ], [ ] );
        
    elif Length( arg ) = 2 and IsList( arg[ 1 ] ) and IsList( arg[ 2 ] ) then
        
        if IsEmpty( arg[ 1 ] ) or ForAny( arg[ 1 ], IsEmpty ) then
            
            cone := rec( generating_vertices := [ ],
                        generating_rays := [ ],
                        matrix:= arg[ 1 ],
                        number_type := "rational",
                        rep_type := "V-rep" );
            
            ObjectifyWithAttributes( cone, TheTypeOfPolymakeCone );
            
            return cone;
            
        fi;
        
        if not IsMatrix( arg[ 1 ] ) then
            
            Error( "Wronge input: The first argument should be a Gap matrix!" );
            
        fi;
        
        if not ForAll( arg[ 1 ], row -> row[ 1 ] in [ 0, 1 ] ) then
            
            Error( "Wronge input: Please see the documentation!" );
            
        fi;
        
        if not ForAll( arg[ 2 ], i -> i in [ 1 .. NrRows( arg[ 1 ] ) ] ) then
            
            Error( "Wronge input for lineality" );
        
        fi;
        
        dim := Length( arg[ 1 ][ 1 ] ) - 1;
        matrix := Filtered( arg[ 1 ], row -> not IsZero( row ) );
        
        if IsEmpty( matrix ) then
            
            Error( "Wronge input: Please make sure the input has sensable direction vectors!" );
            
        fi;
        
        temp := GeneratingVerticesAndGeneratingRays( arg[ 1 ], arg[ 2 ] );
        
        cone := rec( generating_vertices := temp[ 1 ],
                    generating_rays := temp[ 2 ],
                    matrix :=arg[ 1 ],
                    lineality := arg[ 2 ],
                    number_type := "rational",
                    rep_type := "V-rep" );
                    
        ObjectifyWithAttributes( cone, TheTypeOfPolymakeCone );
        
        return cone;
        
    fi;
    
end );

InstallGlobalFunction( Polymake_ConeFromInequalities,
  function( arg )
    local cone, i, temp, matrix, dim;
    
    if Length( arg ) = 0 then
        
        Error( "Wronge input: Please provide some input!" );
        
    elif Length( arg ) = 1 and IsList( arg[ 1 ] ) then
        
        return Polymake_ConeFromInequalities( arg[ 1 ], [ ] );
        
    elif Length( arg ) = 2 and IsList( arg[ 1 ] ) and IsList( arg[ 2 ] ) then
        
        if IsEmpty( arg[ 1 ] ) or ForAny( arg[ 1 ], IsEmpty ) then 
            
            Error( "Wronge input: Please remove the empty lists from the input!" );
            
        fi;
        
        if not IsMatrix( arg[ 1 ] ) then
            
            Error( "Wronge input: The first argument should be a Gap matrix!" );
            
        fi;
        
        if not ForAll( arg[ 2 ], i -> i in [ 1 .. NrRows( arg[ 1 ] ) ] ) then
            
            Error( "Wronge input for lineality" );
            
        fi;
        
        dim := Length( arg[ 1 ][ 1 ] ) - 1;
        
        matrix := Filtered( arg[ 1 ], row -> not IsZero( row ) );
        
        if IsEmpty( matrix ) then
            
            matrix := [ Concatenation( [ 1 ], ListWithIdenticalEntries( dim, 0 ) ) ];
            
        fi;
        
        temp := InequalitiesAndEqualities( matrix, arg[ 2 ] );
        
        cone := rec( matrix := matrix,
                    inequalities := temp[ 1 ],
                    equalities := temp[ 2 ],
                    lineality := arg[ 2 ],
                    number_type := "rational",
                    rep_type := "H-rep" );
        
        ObjectifyWithAttributes( cone, TheTypeOfPolymakeCone );
        
        return cone;
        
    fi;
    
end );


##############################################################################################
##
##  Attributes of PolymakeCones
##
##############################################################################################

InstallMethod( Polymake_V_Rep,
               [ IsPolymakeCone ],
  function( cone )
    local command_string, s, res_string, rays, lineality, lineality_gens;
    
    if cone!.rep_type = "V-rep" then
        
        return cone;
        
    #elif Polymake_IsEmpty( cone ) then
        
    #    return Polymake_ConeByGenerators( [ ] );
        
    else
        
        # compute rays
        command_string := Concatenation( Polymake_H_Rep_command_string( cone ), ".RAYS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        rays := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # extract lineality
        command_string := Concatenation( Polymake_H_Rep_command_string( cone ), ".LINEALITY_SPACE" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        lineality_gens := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # convert lineality -- for NConvex we add the lineality generators to the rays and add a second list, which label the position of the lineality-generators
        lineality := [ 1 .. Length( lineality_gens ) ];
        rays := Concatenation( lineality_gens, rays );
        
        # return the V-representation
        return Polymake_ConeByGenerators( rays, lineality );
        
    fi;
    
end );


InstallMethod( Polymake_H_Rep,
               [ IsPolymakeCone ],
  function( cone )
    local command_string, s, res_string, ineqs, eqs_gens, eqs;
    #dir, file, output;
    
    if cone!.rep_type = "H-rep" then
        
        return cone;
        
    #elif Polymake_IsEmpty( cone ) then
        
    #    return Polymake_ConeFromInequalities( [ ] );
        
    else
        
        if cone!.rep_type = "V-rep" and cone!.matrix = [] then
            return Polymake_ConeFromInequalities( [ [ 0, 1 ], [ -1, -1 ] ] );
        fi;
        
        # compute facets
        command_string := Concatenation( Polymake_V_Rep_command_string( cone ), ".FACETS" );
        
        #dir := Directory( "/home/i" );
        #file := Filename( dir, "test.txt" );
        #output := OutputTextFile( file, true );;
        #AppendTo( output, "Next test:\n\n" );
        #AppendTo( output, command_string );
        #AppendTo( output, "\n\n" );
        #CloseStream(output);
        
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        ineqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # compute linear span
        command_string := Concatenation( Polymake_V_Rep_command_string( cone ), ".LINEAR_SPAN" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        eqs_gens := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # convert eqs -- for NConvex we add the eqs to the ineqs and add a second list, which label the position of the equations
        eqs := [ 1 .. Length( eqs_gens ) ];
        ineqs := Concatenation( eqs_gens, ineqs );
        
        # return cone by inequalities
        return Polymake_ConeFromInequalities( ineqs, eqs );
        
    fi;
    
end );


InstallMethod( Polymake_AmbientSpaceDimension,
              "finding the dimension of the ambient space of the cone",
              [ IsPolymakeCone ],
  function( cone )
    
    return Length( Polymake_H_Rep( cone )!.matrix[1] )-1;
    
end );


InstallMethod( Polymake_Dimension,
              " returns the dimension of the cone",
            [ IsPolymakeCone ],
  function( cone )
    local help_cone, command_string, s;
    
    if Polymake_IsEmpty( cone ) then 
        return -1;
    fi;
    
    # compute v-representation
    help_cone := Polymake_V_Rep( cone );
    
    # produce command string
    command_string := Concatenation( Polymake_V_Rep_command_string( help_cone ), ".CONE_DIM" );
    
    # issue command in Julia and fetch result
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
    
    return EvalString( s );
    
end );


InstallMethod( Polymake_GeneratingVertices,
              " return the list of generating vertices",
              [ IsPolymakeCone ],
  function( cone )
    
    return Set( Polymake_V_Rep( cone )!.generating_vertices );
    
end );


InstallMethod( Polymake_GeneratingRays,
              " return the list of generating vertices",
              [ IsPolymakeCone ],
  function( cone )
    
    return Set( Polymake_V_Rep( cone )!.generating_rays );
    
end );


InstallMethod( Polymake_Equalities,
              " return the list of equalities of a cone",
              [ IsPolymakeCone ],
  function( cone )
    
    return Set( ( Polymake_H_Rep( cone ) )!.equalities );
    
end );


InstallMethod( Polymake_Inequalities,
              " return the list of inequalities of a cone",
              [ IsPolymakeCone ],
  function( cone )
    
    return Set( ( Polymake_H_Rep( cone ) )!.inequalities );
    
end );


InstallMethod( Polymake_RaysInFacets,
              " returns the incident matrix of the rays in the facets",
            [ IsPolymakeCone ],
  function( cone )
    local help_cone, command_string, s, res_string, number_rays, ray_list, i, dummy, j, helper;
    #dir, file, output;
    
    # compute V-representation
    help_cone := Polymake_V_Rep( cone );
    
    # produce command string
    command_string := Concatenation( Polymake_V_Rep_command_string( help_cone ), ".RAYS_IN_FACETS" );
    
    # issue command in Julia and fetch result
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
    
    # process the string
    res_string := SplitString( s, '\n' );
    res_string := List( [ 2 .. Length( res_string ) ], i -> ReplacedString( ReplacedString( ReplacedString( res_string[ i ], " ", "," ), "{", "[" ), "}", "]" ) );
    res_string := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
    number_rays := Length( help_cone!.generating_rays );
    ray_list := [];
    for i in [ 1 .. Length( res_string ) ] do
        dummy := [];
        for j in [ 1 .. Length( res_string[ i ] ) ] do
            helper := List( [ 1 .. number_rays ], i -> 0 );
            helper[ res_string[ i ][ j ] + 1 ] := 1;
            Append( dummy, [ helper ] );
        od;
        Append( ray_list, dummy );
    od;
    
    #dir := Directory( "/home/i" );
    #file := Filename( dir, "test.txt" );
    #output := OutputTextFile( file, true );;
    #AppendTo( output, "\n" );
    #AppendTo( output, "InRaysInFacets\n" );
    #AppendTo( output, command_string );
    #AppendTo( output, "\n" );
    #AppendTo( output, s );
    #    AppendTo( output, "\n" );
    #AppendTo( output, res_string );
    #AppendTo( output, "\n" );
    #AppendTo( output, String( number_rays ) );
    #AppendTo( output, "\n" );
    #AppendTo( output, String( ray_list ) );
    #AppendTo( output, "\n\n" );
    #CloseStream(output);
    
    return ray_list;
    
end );


##############################################################################################
##
##  Properties of PolymakeCones
##
##############################################################################################

InstallMethod( Polymake_IsEmpty,
               "finding if the cone empty is or not",
               [ IsPolymakeCone ],
  function( cone )
    
    return Length( Polymake_V_Rep( cone )!.matrix ) = 0;
    
end );


InstallMethod( Polymake_IsPointed,
               "finding if the cone is pointed or not",
               [ IsPolymakeCone ],
  function( cone )
    local help_cone, rays, lin, command_string, s;
    
    # compute V-representation
    help_cone := Polymake_V_Rep( cone );
    
    # parse the rays into format recognized by Polymake
    command_string := Concatenation( Polymake_V_Rep_command_string( help_cone ), ".POINTED" );
    
    # issue command in Julia and fetch result as string
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.ConeByGAP4PackageConvex ) );
    return EvalString( s );
    
end );



##############################################################################################
##
##  Command strings
##
##############################################################################################

InstallMethod( Polymake_H_Rep_command_string,
               "construct command string for H-Representation of Cone in Julia",
               [ IsPolymakeCone ],
  function( cone )
    local ineqs, eqs, command_string;
    
        # check if the given cone is a V-rep
        if not ( cone!.rep_type = "H-rep" ) then
            return fail;
        fi;
        
        # prepare string with inequalities
        ineqs := cone!.inequalities;
        ineqs := List( [ 1 .. Length( ineqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( ineqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "ConeByGAP4PackageConvex", " = Julia.Polymake.polytope.Cone( INEQUALITIES = [ ", JoinStringsWithSeparator( ineqs, "; " ), " ] " );
        
        # check if we also need equalities
        eqs := cone!.equalities;
        if ( Length( eqs ) > 0 ) then
            eqs := List( [ 1 .. Length( eqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( eqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
            command_string := Concatenation( command_string, " EQUALITIES = [ ", JoinStringsWithSeparator( eqs, "; " ), " ] " );
        fi;
        
        # append closing bracket
        command_string := Concatenation( command_string, ")" );
        
        # add return
        return command_string;
        
end );

InstallMethod( Polymake_V_Rep_command_string,
               "construct command string for V-Representation of Cone in Julia",
               [ IsPolymakeCone ],
  function( cone )
    local rays, lin, command_string;
        
        # check if the given cone is a V-rep
        if not ( cone!.rep_type = "V-rep" ) then
            return "fail";
        fi;
        
        # prepare string with rays
        rays := cone!.generating_rays;
        rays := List( [ 1 .. Length( rays ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( rays[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "ConeByGAP4PackageConvex", " = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ", JoinStringsWithSeparator( rays, "; " ), "] " );
        
        # see if we need lineality
        lin := cone!.lineality;
        if ( Length( lin ) > 0 ) then
            lin := List( [ 1 .. Length( lin ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( lin[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
            command_string := Concatenation( command_string, " INPUT_LINEALITY = [ ", JoinStringsWithSeparator( lin, "; " ), " ] " );
        fi;
        
        # append closing bracket
        command_string := Concatenation( command_string, ")" );
        
        # add return
        return command_string;
    
end );
