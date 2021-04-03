#############################################################################
##
##  ExternalPolymakePolytope.gd  Convex package
##                               Martin Bies
##
##  Copyright 2021               University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
##  Chapter Polytopes in Polymake
##
#############################################################################


##############################################################################################
##
##  Section GAP category of PolymakePolytopes
##
##############################################################################################

DeclareRepresentation( "IsPolymakePolytopeRep", IsPolymakePolytope and IsAttributeStoringRep, [ ] );

BindGlobal( "TheFamilyOfPolymakePolytopes", NewFamily( "TheFamilyOfPolymakePolytopes" ) );

BindGlobal( "TheTypeOfPolymakePolytope", NewType( TheFamilyOfPolymakePolytopes, IsPolymakePolytopeRep ) );


##############################################################################################
##
##  Constructors for PolymakePolytopes
##
##############################################################################################


InstallGlobalFunction( Polymake_PolytopeByGenerators,
  function( arg )
    local poly, i, matrix, temp, dim;
    
    if Length( arg )= 0 then
        
        Error( "Wronge input" );
        
    elif Length( arg ) = 1 and IsList( arg[1] ) then
        
        return Polymake_PolytopeByGenerators( arg[ 1 ], [ ] );
        
    elif Length( arg ) = 2 and IsList( arg[ 1 ] ) and IsList( arg[ 2 ] ) then
        
        if IsEmpty( arg[ 1 ] ) or ForAny( arg[ 1 ], IsEmpty ) then
            
            poly := rec( generating_vertices := [ ],
                        generating_rays := [ ],
                        matrix:= arg[ 1 ],
                        number_type := "rational",
                        rep_type := "V-rep" );
            
            ObjectifyWithAttributes( poly, TheTypeOfPolymakePolytope );
            
            return poly;
            
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
        
        poly := rec( generating_vertices := temp[ 1 ],
                    generating_rays := temp[ 2 ],
                    matrix :=arg[ 1 ],
                    lineality := arg[ 2 ],
                    number_type := "rational",
                    rep_type := "V-rep" );
                    
        ObjectifyWithAttributes( poly, TheTypeOfPolymakePolytope );
        
        return poly;
        
    fi;
    
end );

InstallGlobalFunction( Polymake_PolytopeFromInequalities,
  function( arg )
    local poly, i, temp, matrix, dim;
    
    if Length( arg ) = 0 then
        
        Error( "Wronge input: Please provide some input!" );
        
    elif Length( arg ) = 1 and IsList( arg[ 1 ] ) then
        
        return Polymake_PolytopeFromInequalities( arg[ 1 ], [ ] );
        
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
        
        poly := rec( matrix := matrix,
                    inequalities := temp[ 1 ],
                    equalities := temp[ 2 ],
                    lineality := arg[ 2 ],
                    number_type := "rational",
                    rep_type := "H-rep" );
        
        ObjectifyWithAttributes( poly, TheTypeOfPolymakePolytope );
        
        return poly;
        
    fi;
    
end );


##############################################################################################
##
##  Attributes of PolymakeCones
##
##############################################################################################

InstallMethod( Polymake_V_Rep,
               [ IsPolymakePolytope ],
  function( poly )
    local ineqs, eqs, command_string, s, rays, vertices;
    
    if poly!.rep_type = "V-rep" then
        
        return poly;
        
    else
        
        # compute rays
        command_string := Concatenation( Polymake_H_Rep_command_string( poly ), ".VERTICES" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        rays := EvalString( s );
        
        # compute vertices
        command_string := Concatenation( Polymake_H_Rep_command_string( poly ), ".LINEALITY_SPACE" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        vertices := EvalString( s );
        
        # return the V-representation
        return Polymake_PolytopeByGenerators( rays, vertices );
        
    fi;
    
end );


InstallMethod( Polymake_H_Rep,
               [ IsPolymakePolytope ],
  function( poly )
    local command_string, s, res_string, ineqs, eqs, dir, file, output;
    
    if poly!.rep_type = "H-rep" then
        
        return poly;
        
    else
        
        if poly!.rep_type = "V-rep" and poly!.matrix = [] then
            return Polymake_PolytopeFromInequalities( [ [ 0, 1 ], [ -1, -1 ] ] );
        fi;
        
        # compute facets
        command_string := Concatenation( Polymake_V_Rep_command_string( poly ), ".FACETS" );
        
        dir := Directory( "/home/i" );
        file := Filename( dir, "test.txt" );
        output := OutputTextFile( file, true );;
        AppendTo( output, command_string );
        CloseStream(output);
        
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        ineqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # compute linear span
        command_string := Concatenation( Polymake_V_Rep_command_string( poly ), ".AFFINE_HULL" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        eqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # return poly by inequalities
        return Polymake_PolytopeFromInequalities( ineqs, eqs );
        
    fi;
    
end );


InstallMethod( Polymake_AmbientSpaceDimension,
              "finding the dimension of the ambient space of the poly",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Length( Polymake_H_Rep( poly )!.matrix[1] )-1;
    
end );


InstallMethod( Polymake_Dimension,
              " returns the dimension of the poly",
            [ IsPolymakePolytope ],
  function( poly )
    local help_poly, rays, lin, command_string, s;
    
    if Polymake_IsEmpty( poly ) then 
        return -1;
    fi;
    
    # compute v-representation
    help_poly := Polymake_V_Rep( poly );
    
    # produce command string
    command_string := Concatenation( Polymake_V_Rep_command_string( help_poly ), ".POLYTOPE_DIM" );
    
    # issue command in Julia and fetch result
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
    return EvalString( s );
    
end );


InstallMethod( Polymake_GeneratingVertices,
              " return the list of generating vertices",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Set( Polymake_V_Rep( poly )!.generating_vertices );
    
end );


InstallMethod( Polymake_GeneratingRays,
              " return the list of generating vertices",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Set( Polymake_V_Rep( poly )!.generating_rays );
    
end );


InstallMethod( Polymake_Equalities,
              " return the list of equalities of a poly",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Set( ( Polymake_H_Rep( poly ) )!.equalities );
    
end );


InstallMethod( Polymake_Inequalities,
              " return the list of inequalities of a poly",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Set( ( Polymake_H_Rep( poly ) )!.inequalities );
    
end );


##############################################################################################
##
##  Properties of PolymakeCones
##
##############################################################################################

InstallMethod( Polymake_IsEmpty,
               "finding if the poly empty is or not",
               [ IsPolymakePolytope ],
  function( poly )
    
    return Length( Polymake_V_Rep( poly )!.matrix ) = 0;
    
end );


InstallMethod( Polymake_IsPointed,
               "finding if the poly is pointed or not",
               [ IsPolymakePolytope ],
  function( poly )
    local help_poly, rays, lin, command_string, s;
    
    # compute V-representation
    help_poly := Polymake_V_Rep( poly );
    
    # parse the rays into format recognized by Polymake
    command_string := Concatenation( Polymake_V_Rep_command_string( help_poly ), ".POINTED" );
    
    # issue command in Julia and fetch result as string
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
    return EvalString( s );
    
end );


InstallMethod( Polymake_IsBounded,
              " returns if the polytope is bounded or not",
              [ IsPolymakePolytope ],
  function( poly )
    local help_poly, command_string, s;
    
    help_poly := Polymake_H_Rep( poly );
    command_string := Concatenation( Polymake_H_Rep_command_string( help_poly ), ".BOUNDED" );
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
    return EvalString( s );
    
end );


##############################################################################################
##
##  Command strings
##
##############################################################################################

InstallMethod( Polymake_H_Rep_command_string,
               "construct command string for H-Representation of Cone in Julia",
               [ IsPolymakePolytope ],
  function( poly )
    local ineqs, eqs, command_string;
    
        # check if the given poly is a V-rep
        if not ( poly!.rep_type = "H-rep" ) then
            return fail;
        fi;
        
        # prepare string with inequalities
        ineqs := poly!.inequalities;
        ineqs := List( [ 1 .. Length( ineqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( ineqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "PolytopeByGAP4PackageConvex", " = Julia.Polymake.polytope.Polytope( INEQUALITIES = [ ", JoinStringsWithSeparator( ineqs, "; " ), " ] " );
        
        # check if we also need equalities
        eqs := poly!.equalities;
        if ( Length( eqs ) > 0 ) then
            eqs := List( [ 1 .. Length( eqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( eqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
            command_string := Concatenation( command_string, " EQUATIONS = [ ", JoinStringsWithSeparator( eqs, "; " ), " ] " );
        fi;
        
        # append closing bracket
        command_string := Concatenation( command_string, ")" );
        
        # add return
        return command_string;
        
end );

InstallMethod( Polymake_V_Rep_command_string,
               "construct command string for V-Representation of Cone in Julia",
               [ IsPolymakePolytope ],
  function( poly )
    local vertices, lin, command_string, dir, file, output;
        
        # check if the given poly is a V-rep
        if not ( poly!.rep_type = "V-rep" ) then
            return "fail";
        fi;
        
        dir := Directory( "/home/i" );
        file := Filename( dir, "test.txt" );
        output := OutputTextFile( file, true );;
        AppendTo( output, String( poly!.generating_vertices ) );
        AppendTo( output, "\n" );
        AppendTo( output, String( poly!.generating_rays ) );
        AppendTo( output, "\n" );
        CloseStream(output);
        
        # prepare string with vertices -- as polymake considers them affine, we have to add a 1 at the beginning
        vertices := poly!.matrix;
        #vertices := Filtered( vertices, v -> not IsZero( v ) );
        vertices := List( [ 1 .. Length( vertices ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( vertices[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "PolytopeByGAP4PackageConvex", " = Julia.Polymake.polytope.Polytope( POINTS = [ ", JoinStringsWithSeparator( vertices, "; " ), "] " );
        
        # see if we need lineality
        lin := poly!.lineality;
        if ( Length( lin ) > 0 ) then
            lin := List( [ 1 .. Length( lin ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( lin[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
            command_string := Concatenation( command_string, " INPUT_LINEALITY = [ ", JoinStringsWithSeparator( lin, "; " ), " ] " );
        fi;
        
        # append closing bracket
        command_string := Concatenation( command_string, ")" );
        
        # add return
        return command_string;
        
end );
