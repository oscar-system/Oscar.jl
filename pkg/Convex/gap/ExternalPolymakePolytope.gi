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
    
    if Length( arg )= 0 or ForAll( arg, IsEmpty ) then
        
        Error( "Wronge input: Please provide some input!" );
        
    elif Length( arg ) = 1 and IsList( arg[1] ) then
        
        return Polymake_PolytopeByGenerators( arg[ 1 ], [ ] );
        
    elif Length( arg ) = 2 and IsList( arg[ 1 ] ) and IsList( arg[ 2 ] ) then
        
        if ( not IsEmpty( arg[ 1 ] ) ) and not ( IsMatrix( arg[ 1 ] ) ) then
            Error( "Wronge input: The first argument should be a Gap matrix!" );
        fi;
        
        if ( not IsEmpty( arg[ 2 ] ) ) and not ( IsMatrix( arg[ 2 ] ) ) then
            Error( "Wronge input: The second argument should be a Gap matrix!" );
        fi;
        
        poly := rec( vertices := arg[ 1 ],
                     lineality := arg[ 2 ],
                     number_type := "rational",
                     rep_type := "V-rep" );
        ObjectifyWithAttributes( poly, TheTypeOfPolymakePolytope );
        return Polymake_CanonicalPolytopeByGenerators( poly );
        
    fi;
    
end );


InstallGlobalFunction( Polymake_PolytopeFromInequalities,
  function( arg )
    local poly, i, temp, matrix, dim;
    
    if Length( arg ) = 0 or ForAll( arg, IsEmpty ) then
        
        Error( "Wronge input: Please provide some input!" );
        
    elif Length( arg ) = 1 and IsList( arg[ 1 ] ) then
        
        return Polymake_PolytopeFromInequalities( arg[ 1 ], [ ] );
        
    elif Length( arg ) = 2 and IsList( arg[ 1 ] ) and IsList( arg[ 2 ] ) then
        
        if ( not IsEmpty( arg[ 1 ] ) ) and not ( IsMatrix( arg[ 1 ] ) ) then
            Error( "Wronge input: The first argument should be a Gap matrix!" );
        fi;
        
        if ( not IsEmpty( arg[ 2 ] ) ) and not ( IsMatrix( arg[ 2 ] ) ) then
            Error( "Wronge input: The second argument should be a Gap matrix!" );
        fi;
        
        poly := rec( inequalities := arg[ 1 ],
                     equalities := arg[ 2 ],
                     number_type := "rational",
                     rep_type := "H-rep" );
        ObjectifyWithAttributes( poly, TheTypeOfPolymakePolytope );
        return Polymake_CanonicalPolytopeFromInequalities( poly );
        
    fi;
    
end );


##############################################################################################
##
##  Canonicalize polytopes
##
##############################################################################################

InstallMethod( Polymake_CanonicalPolytopeByGenerators,
               [ IsPolymakePolytope ],
  function( poly )
    local command_string, s, res_string, vertices, v_copy, scaled_vertices, i, scale, lineality, scaled_lineality, new_poly;
    
    if poly!.rep_type = "H-rep" then
        
        return fail;
        
    else
        
        # compute vertices
        command_string := Concatenation( Polymake_V_Rep_command_string( poly ), ".VERTICES" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        vertices := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational vertices - we turn them into integral vectors
        # also, Polymake requires x0 = 1 in affine coordinates - we remove this 1
        scaled_vertices := [];
        for i in [ 1 .. Length( vertices ) ] do
            scale := Lcm( List( vertices[ i ], r -> DenominatorRat( r ) ) );
            v_copy := ShallowCopy( vertices[ i ] );
            Remove( v_copy, 1 );
            Append( scaled_vertices, [ scale * v_copy ] );
        od;
        
        # extract lineality
        command_string := Concatenation( Polymake_V_Rep_command_string( poly ), ".LINEALITY_SPACE" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        lineality := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational lineality - we turn them into integral vectors
        scaled_lineality := [];
        for i in [ 1 .. Length( lineality ) ] do
            scale := Lcm( List( lineality[ i ], r -> DenominatorRat( r ) ) );
            v_copy := ShallowCopy( lineality[ i ] );
            Remove( v_copy, 1 );
            Append( scaled_lineality, [ scale * v_copy ] );
        od;
        
        # construct the new poly
        new_poly := rec( vertices := scaled_vertices,
                         lineality := scaled_lineality,
                         number_type := "rational",
                         rep_type := "V-rep" );
        ObjectifyWithAttributes( new_poly, TheTypeOfPolymakePolytope );
        return new_poly;
        
    fi;
    
end );

InstallMethod( Polymake_CanonicalPolytopeFromInequalities,
               [ IsPolymakePolytope ],
  function( poly )
    local command_string, s, res_string, ineqs, scaled_ineqs, i, scale, eqs, scaled_eqs, new_poly;
    
    if poly!.rep_type = "V-rep" then
        
        return fail;
        
    else
        
        # compute facets
        command_string := Concatenation( Polymake_H_Rep_command_string( poly ), ".FACETS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        ineqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational facets - we turn them into integral vectors
        scaled_ineqs := [];
        for i in [ 1 .. Length( ineqs ) ] do
            scale := Lcm( List( ineqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_ineqs, [ scale * ineqs[ i ] ] );
        od;
        
        # compute affine hull
        command_string := Concatenation( Polymake_H_Rep_command_string( poly ), ".AFFINE_HULL" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        eqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational affine hulls - we turn them into integral vectors
        scaled_eqs := [];
        for i in [ 1 .. Length( eqs ) ] do
            scale := Lcm( List( eqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_eqs, [ scale * eqs[ i ] ] );
        od;
        
        # construct the new poly
        new_poly := rec( inequalities := scaled_ineqs,
                         equalities := scaled_eqs,
                         number_type := "rational",
                         rep_type := "H-rep" );
        ObjectifyWithAttributes( new_poly, TheTypeOfPolymakePolytope );
        return new_poly;
        
    fi;
    
end );


##############################################################################################
##
##  Conversion of polytopes
##
##############################################################################################

InstallMethod( Polymake_V_Rep,
               [ IsPolymakePolytope ],
  function( poly )
    local command_string, s, res_string, vertices, v_copy, scaled_vertices, i, scale, lineality, scaled_lineality, new_poly;
    
    if poly!.rep_type = "V-rep" then
        
        return poly;
        
    else
        
        # compute vertices
        command_string := Concatenation( Polymake_H_Rep_command_string( poly ), ".VERTICES" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        vertices := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational vertices - we turn them into integral vectors
        scaled_vertices := [];
        for i in [ 1 .. Length( vertices ) ] do
            scale := Lcm( List( vertices[ i ], r -> DenominatorRat( r ) ) );
            v_copy := ShallowCopy( vertices[ i ] );
            Remove( v_copy, 1 );
            Append( scaled_vertices, [ scale * v_copy ] );
        od;
        
        # compute lineality
        command_string := Concatenation( Polymake_H_Rep_command_string( poly ), ".LINEALITY_SPACE" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        lineality := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational lineality - we turn them into integral vectors
        scaled_lineality := [];
        for i in [ 1 .. Length( lineality ) ] do
            scale := Lcm( List( lineality[ i ], r -> DenominatorRat( r ) ) );
            v_copy := ShallowCopy( lineality[ i ] );
            Remove( v_copy, 1 );
            Append( scaled_lineality, [ scale * v_copy ] );
        od;
        
        # construct the new poly
        new_poly := rec( vertices := scaled_vertices,
                         lineality := scaled_lineality,
                         number_type := "rational",
                         rep_type := "V-rep" );
        ObjectifyWithAttributes( new_poly, TheTypeOfPolymakePolytope );
        return new_poly;
        
    fi;
    
end );


InstallMethod( Polymake_H_Rep,
               [ IsPolymakePolytope ],
  function( poly )
    local command_string, s, res_string, ineqs, i, scale, scaled_ineqs, eqs, scaled_eqs, new_poly;
    
    if poly!.rep_type = "H-rep" then
        
        return poly;
        
    else
        
        if poly!.rep_type = "V-rep" and poly!.matrix = [] then
            return Polymake_PolytopeFromInequalities( [ [ 0, 1 ], [ -1, -1 ] ] );
        fi;
        
        # compute inequalities
        command_string := Concatenation( Polymake_V_Rep_command_string( poly ), ".FACETS" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        ineqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational facets - we turn them into integral vectors
        scaled_ineqs := [];
        for i in [ 1 .. Length( ineqs ) ] do
            scale := Lcm( List( ineqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_ineqs, [ scale * ineqs[ i ] ] );
        od;
        
        # compute equalities
        command_string := Concatenation( Polymake_V_Rep_command_string( poly ), ".AFFINE_HULL" );
        JuliaEvalString( command_string );
        s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
        res_string := SplitString( s, '\n' );
        res_string := List( [ 2 .. Length( res_string ) ], i -> Concatenation( "[", ReplacedString( res_string[ i ], " ", "," ), "]" ) );
        eqs := EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
        
        # sometimes, Polymake returns rational affine hulls - we turn them into integral vectors
        scaled_eqs := [];
        for i in [ 1 .. Length( eqs ) ] do
            scale := Lcm( List( eqs[ i ], r -> DenominatorRat( r ) ) );
            Append( scaled_eqs, [ scale * eqs[ i ] ] );
        od;
        
        # construct the new poly
        new_poly := rec( inequalities := scaled_ineqs,
                         equalities := scaled_eqs,
                         number_type := "rational",
                         rep_type := "H-rep" );
        ObjectifyWithAttributes( new_poly, TheTypeOfPolymakePolytope );
        return new_poly;
        
    fi;
    
end );


##############################################################################################
##
##  Attributes of PolymakeCones
##
##############################################################################################

InstallMethod( Polymake_AmbientSpaceDimension,
              "finding the dimension of the ambient space of the poly",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Length( Polymake_V_Rep( poly )!.vertices[0] );
    
end );


InstallMethod( Polymake_Dimension,
              " returns the dimension of the poly",
            [ IsPolymakePolytope ],
  function( poly )
    local command_string, s;
    
    if Polymake_IsEmpty( poly ) then
        return -1;
    fi;
    
    # produce command string
    command_string := Concatenation( Polymake_V_Rep_command_string( Polymake_V_Rep( poly ) ), ".CONE_DIM" );
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
    return EvalString( s ) - 1;
    
end );


InstallMethod( Polymake_Vertices,
              " return the list of generating vertices",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Set( Polymake_V_Rep( poly )!.vertices );
    
end );


InstallMethod( Polymake_Linealities,
              " return the list of linealities",
              [ IsPolymakePolytope ],
  function( poly )
    
    return Set( Polymake_V_Rep( poly )!.linealities );
    
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


InstallMethod( Polymake_LatticePoints,
              " return the list of the lattice points of poly",
              [ IsPolymakePolytope ],
  function( poly )
    local help_poly, command_string, s, res_string;
    
    # compute v-representation
    help_poly := Polymake_V_Rep( poly );
    
    # produce command string
    command_string := Concatenation( Polymake_V_Rep_command_string( help_poly ), ".LATTICE_POINTS_GENERATORS" );
    
    # issue command in Julia and fetch result
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
    res_string := SplitString( s, '\n' );
    res_string := List( [ 2 .. Length( res_string ) - 3 ], i -> Concatenation( "[", ReplacedString( ReplacedString( res_string[ i ], " ", "," ), "<", "" ), "]" ) );
    return EvalString( Concatenation( "[", JoinStringsWithSeparator( res_string, "," ), "]" ) );
    
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
    
    return Length( Polymake_V_Rep( poly )!.vertices ) = 0;
    
end );


InstallMethod( Polymake_IsPointed,
               "finding if the poly is pointed or not",
               [ IsPolymakePolytope ],
  function( poly )
    local command_string, s;
    
    command_string := Concatenation( Polymake_V_Rep_command_string( Polymake_V_Rep( poly ) ), ".POINTED" );
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
    return EvalString( s );
    
end );


InstallMethod( Polymake_IsBounded,
              " returns if the polytope is bounded or not",
              [ IsPolymakePolytope ],
  function( poly )
    local command_string, s;
    
    command_string := Concatenation( Polymake_H_Rep_command_string( Polymake_H_Rep( poly ) ), ".BOUNDED" );
    JuliaEvalString( command_string );
    s := JuliaToGAP( IsString, Julia.string( Julia.PolytopeByGAP4PackageConvex ) );
    return EvalString( s );
    
end );


##############################################################################################
##
##  Command strings
##
##############################################################################################

InstallMethod( Polymake_V_Rep_command_string,
               "construct command string for V-Representation of polytope in Julia",
               [ IsPolymakePolytope ],
  function( poly )
    local vertices, lin, new_vertices, new_lin, i, v, command_string;
        
        # check if the given poly is a V-rep
        if not ( poly!.rep_type = "V-rep" ) then
            return "fail";
        fi;
        
        # extract data
        vertices := poly!.vertices;
        lin := poly!.lineality;
        
        # add 1s at the beginning, since Polymake always considers polytopes as intersections with the hyperplane x0 = 1
        new_vertices := [];
        for i in [ 1 .. Length( vertices ) ] do
            v := ShallowCopy( vertices[ i ] );
            Add( v, 1, 1 );
            Add( new_vertices, v );
        od;
        new_lin := [];
        for i in [ 1 .. Length( lin ) ] do
            v := ShallowCopy( lin[ i ] );
            Add( v, 1, 1 );
            Add( new_lin, v );
        od;
        
        # check for degenerate case
        if Length( new_vertices ) = 0 then
            new_vertices := new_lin;
        fi;
        
        # prepare string with vertices -- as polymake considers them affine, we have to add a 1 at the beginning
        new_vertices := List( [ 1 .. Length( new_vertices ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( new_vertices[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "PolytopeByGAP4PackageConvex", " = Julia.Polymake.polytope.Polytope( POINTS = [ ", JoinStringsWithSeparator( new_vertices, "; " ), "] " );
        
        # see if we need lineality
        if ( Length( new_lin ) > 0 ) then
            new_lin := List( [ 1 .. Length( new_lin ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( new_lin[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
            command_string := Concatenation( command_string, ", INPUT_LINEALITY = [ ", JoinStringsWithSeparator( new_lin, "; " ), " ] " );
        fi;
        
        # return command string
        return Concatenation( command_string, ")" );
        
end );

InstallMethod( Polymake_H_Rep_command_string,
               "construct command string for H-Representation of polytope in Julia",
               [ IsPolymakePolytope ],
  function( poly )
    local ineqs, eqs, command_string;
    
        # check if the given poly is a V-rep
        if not ( poly!.rep_type = "H-rep" ) then
            return fail;
        fi;
        
        # extract data
        ineqs := poly!.inequalities;
        eqs := poly!.equalities;
        
        # check for degenerate case
        if Length( ineqs ) = 0 then
            ineqs := eqs;
        fi;
        
        # prepare string with inequalities
        ineqs := List( [ 1 .. Length( ineqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( ineqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
        command_string := Concatenation( "PolytopeByGAP4PackageConvex", " = Julia.Polymake.polytope.Polytope( INEQUALITIES = [ ", JoinStringsWithSeparator( ineqs, "; " ), " ] " );
        
        # check if we also need equalities
        if ( Length( eqs ) > 0 ) then
            eqs := List( [ 1 .. Length( eqs ) ], i -> ReplacedString( ReplacedString( ReplacedString( String( eqs[ i ] ), ",", "" ), "[ ", "" ), " ]", "" ) );
            command_string := Concatenation( command_string, ", EQUATIONS = [ ", JoinStringsWithSeparator( eqs, "; " ), " ] " );
        fi;
        
        # return command string
        return Concatenation( command_string, ")" );
        
end );
