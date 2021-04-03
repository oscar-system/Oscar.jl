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
                    
        ObjectifyWithAttributes( poly, TheTypeOfPolymakeCone );
        
        return poly;
        
    fi;
    
end );

InstallGlobalFunction( Polymake_ConeFromInequalities,
  function( arg )
    local poly, i, temp, matrix, dim;
    
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
        
        poly := rec( matrix := matrix,
                    inequalities := temp[ 1 ],
                    equalities := temp[ 2 ],
                    lineality := arg[ 2 ],
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
        ineqs := List( [ 1 .. Length( ineqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( ineqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        eqs := poly!.equalities;
        eqs := List( [ 1 .. Length( eqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( eqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        
        # compute the rays
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INEQUALITIES = [ ",
                                         JoinStringsWithSeparator( ineqs, "; " ),
                                         " ], EQUALITIES = [ ",
                                         JoinStringsWithSeparator( eqs, "; " ),
                                         " ] ).RAYS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        rays := EvalString( s );
        
        # compute the vertices
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INEQUALITIES = [ ",
                                         JoinStringsWithSeparator( ineqs, "; " ),
                                         " ], EQUALITIES = [ ",
                                         JoinStringsWithSeparator( eqs, "; " ),
                                         " ] ).LINEALITY_SPACE" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        vertices := EvalString( s );
        
        # return the V-representation
        return Polymake_ConeByGenerators( rays, vertices );
        
    fi;
    
end );


InstallMethod( Polymake_H_Rep,
               [ IsPolymakeCone ],
  function( cone )
    local rays, lin, command_string, s, res_string, ineqs, eqs;
    
    if cone!.rep_type = "H-rep" then
        
        return cone;
        
    else
        
        if cone!.rep_type = "V-rep" and cone!.matrix = [] then
            return Polymake_ConeFromInequalities( [ [ 0, 1 ], [ -1, -1 ] ] );
        fi;
        
        # prepare string with rays
        rays := cone!.generating_rays;
        rays := List( [ 1 .. Length( rays ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( rays[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        lin := cone!.lineality;
        lin := List( [ 1 .. Length( eqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( lin[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        
        # compute the inequalities
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ",
                                         JoinStringsWithSeparator( rays, "; " ),
                                         " ], INPUT_LINEALITY = [ ",
                                         JoinStringsWithSeparator( lin, "; " ),
                                         " ] ).FACETS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        ineqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # compute the equalities
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ",
                                         JoinStringsWithSeparator( rays, "; " ),
                                         " ], INPUT_LINEALITY = [ ",
                                         JoinStringsWithSeparator( lin, "; " ),
                                         " ] ).LINEAR_SPAN" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        eqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # return cone by inequalities
        return Polymake_ConeFromInequalities( ineqs, eqs );
        
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
  function( cone )
    local help_cone, rays, lin, command_string, s;
    
    if Polymake_IsEmpty( cone ) then 
        return -1;
    else
        
        # compute v-representation
        help_cone := Polymake_V_Rep( cone );
        
        # parse the rays into format recognized by Polymake
        rays := help_cone!.generating_rays;
        rays := List( [ 1 .. Length( rays ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( rays[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        lin := help_cone!.lineality;
        lin := List( [ 1 .. Length( lin ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( lin[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        
        # prepare command string
        command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ",
                                         JoinStringsWithSeparator( rays, "; " ),
                                         " ], INPUT_LINEALITY = [ ",
                                         JoinStringsWithSeparator( lin, "; " ),
                                         " ] ).CONE_DIM" );
        
        # issue command in Julia and fetch result
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
        return EvalString( s );
        
    fi;
    
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
    rays := help_cone!.generating_rays;
    rays := List( [ 1 .. Length( rays ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( rays[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
    lin := help_cone!.lineality;
    lin := List( [ 1 .. Length( lin ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( lin[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
    
    # prepare command string
    command_string := Concatenation( "F = Julia.Polymake.polytope.Cone( INPUT_RAYS = [ ",
                                    JoinStringsWithSeparator( rays, "; " ),
                                    " ], INPUT_LINEALITY = [ ",
                                    JoinStringsWithSeparator( lin, "; " ),
                                    " ] ).POINTED" );
    
    # issue command in Julia and fetch result as string
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.F ) );
    return EvalString( s );
    
end );
