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

DeclareRepresentation( "IsPolymakeCone", IsPolymakeCone and IsAttributeStoringRep, [ ] );

BindGlobal( "TheFamilyOfPolymakeCones", NewFamily( "TheFamilyOfPolymakeCones" ) );

BindGlobal( "TheTypeOfPolymakeCone", NewType( TheFamilyOfPolymakeCones, IsPolymakeConeRep ) );


##############################################################################################
##
##  Constructors for PolymakeCones
##
##############################################################################################

InstallMethod( Polymake_ConeByGenerators,
               " a list of vertices",
               [ IsList ],
  function( arg )
    local poly, i, matrix, temp, dim;
    
    if Length( arg )= 0 then
        
        Error( "Wronge input" );
        
    elif Length( arg ) = 1 and IsList( arg[1] ) then
        
        return Polymake_ConeByGenerators( arg[ 1 ], [ ] );
        
    elif Length( arg ) = 2 and IsList( arg[ 1 ] ) and IsList( arg[ 2 ] ) then
        
        if IsEmpty( arg[ 1 ] ) or ForAny( arg[ 1 ], IsEmpty ) then
            
            poly := rec( generating_vertices := [ ],
                        generating_rays := [ ],
                        matrix:= arg[ 1 ],
                        number_type := "rational",
                        rep_type := "V-rep" );
            
            ObjectifyWithAttributes( poly, TheTypeOfPolymakeCone );
            
            return poly;
            
        fi;
        
        if not IsMatrix( arg[ 1 ] ) then
            
            Error( "Wronge input: The first argument should be a Gap matrix!" );
            
        fi;
        
        if not ForAll( arg[ 1 ], row -> row[ 1 ] in [ 0, 1 ] ) then
            
            Error( "Wronge input: Please see the documentation!" );
            
        fi;
        
        if not ForAll( arg[ 2 ], i -> i in [ 1 .. NrRows( arg[ 1 ] ) ] ) then
            
            Error( "Wronge input for linearity" );
        
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
                    linearity := arg[ 2 ],
                    number_type := "rational",
                    rep_type := "V-rep" );
                    
        ObjectifyWithAttributes( poly, TheTypeOfPolymakeCone );
        
        return poly;
        
    fi;
    
end );

InstallMethod( Polymake_ConeByInequalities,
               " a list of vertices",
               [ IsList ],
  function( arg )
    local poly, i, temp, matrix, dim;
    
    if Length( arg ) = 0 then
        
        Error( "Wronge input: Please provide some input!" );
        
    elif Length( arg ) = 1 and IsList( arg[ 1 ] ) then
        
        return Cdd_PolyhedronByInequalities( arg[ 1 ], [ ] );
        
    elif Length( arg ) = 2 and IsList( arg[ 1 ] ) and IsList( arg[ 2 ] ) then
        
        if IsEmpty( arg[ 1 ] ) or ForAny( arg[ 1 ], IsEmpty ) then 
            
            Error( "Wronge input: Please remove the empty lists from the input!" );
            
        fi;
        
        if not IsMatrix( arg[ 1 ] ) then
            
            Error( "Wronge input: The first argument should be a Gap matrix!" );
            
        fi;
        
        if not ForAll( arg[ 2 ], i -> i in [ 1 .. NrRows( arg[ 1 ] ) ] ) then
            
            Error( "Wronge input for linearity" );
            
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
                    linearity := arg[ 2 ],
                    number_type := "rational",
                    rep_type := "H-rep" );
        
        ObjectifyWithAttributes( poly, TheTypeOfPolymakeCone );
        
        return poly;
        
    fi;
    
end );


##############################################################################################
##
##  Attributes of PolymakeCones
##
##############################################################################################

InstallMethod( Polymake_V_Rep,
               [ IsPolymakeCone ],
  function( poly )
    local ineqs, eqs, command_string, s, rays, vertices;
    
    if poly!.rep_type = "V-rep" then
        
        return poly;
        
    else
        
        # Prepare ineqs and eqs to be parsed
        ineqs := poly!.inequalities;
        ineqs_string_list := List( [ 1 .. Length( ineqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( ineqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        eqs := poly!.equalities;
        eqs_string_list := List( [ 1 .. Length( eqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( eqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        
        # compute the rays
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INEQUALITIES = [ ",
                                         JoinStringsWithSeparator( ineqs_string_list, "; " ),
                                         " ], EQUALITIES = [ ",
                                         JoinStringsWithSeparator( equ_string_list, "; " ),
                                         " ] ).RAYS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        rays := EvalString( s );
        
        # compute the vertices
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INEQUALITIES = [ ",
                                         JoinStringsWithSeparator( ineqs_string_list, "; " ),
                                         " ], EQUALITIES = [ ",
                                         JoinStringsWithSeparator( equ_string_list, "; " ),
                                         " ] ).VERTICES" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        vertices := EvalString( s );
        
        # return the V-representation
        return Polymake_ConeByGenerators( rays, vertices );
        
    fi;
    
end );


InstallMethod( Polymake_H_Rep,
               [ IsPolymakeCone ],
  function( poly )
    local rays, string_list, command_string, s, res_string, ineqs, eqs;
    
    if poly!.rep_type = "H-rep" then
        
        return poly;
        
    else
        
        if poly!.rep_type = "V-rep" and poly!.matrix = [] then
            return Polymake_ConeByInequalities( [ [ 0, 1 ], [ -1, -1 ] ] );
        fi;
        
        # prepare string with rays
        rays := cone!.generating_rays;
        string_list := List( [ 1 .. Length( rays ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( rays[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        
        # compute the inequalities
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ", JoinStringsWithSeparator( string_list, "; " ), " ] ).FACETS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        ineqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # compute the equalities
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ", JoinStringsWithSeparator( string_list, "; " ), " ] ).EQUALITIES" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        eqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # return cone by inequalities
        return Polymake_ConeByInequalities( ineqs, eqs );
        
    fi;
    
end );


InstallMethod( Polymake_AmbientSpaceDimension,
              "finding the dimension of the ambient space of the cone",
              [ IsPolymakeCone ],
  function( poly )
    
    return Length( Polymake_H_Rep( poly )!.matrix[1] )-1;
    
end );


InstallMethod( Polymake_Dimension,
              " returns the dimension of the cone",
            [ IsPolymakeCone ],
  function( poly )
    local help_poly, input_rays, string_list, command_string, s;
    
    if Polymake_IsEmpty( poly ) then 
        return -1;
    else
        
        # compute v-representation
        help_poly := Polymake_V_Rep( poly );
        
        # parse the rays into format recognized by Polymake
        input_rays := help_poly!.generating_rays;
        string_list := List( [ 1 .. Length( input_rays ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( input_rays[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ", JoinStringsWithSeparator( string_list, "; " ), " ] ).CONE_DIM" );
        
        # issue command in Julia and fetch result
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        return EvalString( s );
        
    fi;
    
end );


InstallMethod( Polymake_GeneratingVertices,
              " return the list of generating vertices",
              [ IsPolymakeCone ],
  function( poly )
    
    return Set( Polymake_V_rep( poly )!.generating_vertices );
    
end );


InstallMethod( Polymake_GeneratingRays,
              " return the list of generating vertices",
              [ IsPolymakeCone ],
  function( poly )
    
    return Set( Polymake_V_rep( poly )!.generating_rays );
    
end );


InstallMethod( Polymake_Equalities,
              " return the list of equalities of a cone",
              [ IsPolymakeCone ],
  function( poly )
    
    return Set( Polymake_H_Rep( poly ) )!.equalities );
    
end );


InstallMethod( Polymake_Inequalities,
              " return the list of inequalities of a cone",
              [ IsPolymakeCone ],
  function( poly )
    
    return Set( Polymake_H_Rep( poly ) )!.inequalities );
    
end );


##############################################################################################
##
##  Properties of PolymakeCones
##
##############################################################################################

InstallMethod( Polymake_IsEmpty,
               "finding if the cone empty is or not",
               [ IsPolymakeCone ],
  function( poly )
    
    return Length( Polymake_V_Rep( poly )!.matrix ) = 0;
    
end );


InstallMethod( Polymake_IsPointed,
               "finding if the cone is pointed or not",
               [ IsPolymakeCone ],
  function( poly )
    local help_poly, input_rays, string_list, command_string, s;
    
    if poly!.rep_type = "H-rep" then
        help_poly := Polymake_V_Rep( poly );
    else
        help_poly := poly;
    fi;
    
    # Parse the rays into format recognized by Polymake
    input_rays := cone!.generating_rays;
    string_list := List( [ 1 .. Length( input_rays ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( input_rays[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
    command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ", JoinStringsWithSeparator( string_list, "; " ), " ] ).POINTED" );
    
    # issue command in Julia and fetch result as string
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
    return EvalString( s );
    
end );
