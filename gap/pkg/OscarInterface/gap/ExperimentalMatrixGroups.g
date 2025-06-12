# This file implements a minimal GAP MatrixObj type which wraps around Julia
# matrices, to make them usable within GAP.

BindGlobal("JuliaMatrixFamily", NewFamily("JuliaMatrixFamily"));

DeclareRepresentation(
  "IsJuliaMatrixRep",
  IsComponentObjectRep and IsAttributeStoringRep and IsMatrixObj);
  
JuliaMatrixIdentificationType :=
    NewType( JuliaMatrixFamily, IsJuliaMatrixRep and IsMutable);

BindGlobal("MakeJuliaMatrixRep",function(m)
  if not IsJuliaObject(m) then
    Error("<m> must be a Julia object");
  fi;
  return Objectify(JuliaMatrixIdentificationType,
               rec(m := m));
end);
    
InstallOtherMethod( NumberRows, [IsJuliaMatrixRep], m -> Julia.nrows(m!.m));

InstallOtherMethod( NumberColumns, [IsJuliaMatrixRep], m -> Julia.ncols(m!.m));

InstallOtherMethod( BaseDomain, [IsJuliaMatrixRep], m -> Julia.base_ring(m!.m));

InstallOtherMethod( \[\,\],
    [ IsJuliaMatrixRep, IsPosInt and IsSmallIntRep,
                       IsPosInt and IsSmallIntRep ],
    function( m, i, j )
      return (m!.m)[i,j];
    end );
    
# TODO: Not working. JuliaMatrixRep is immutable
InstallOtherMethod( \[\,\]\:\=,
    [ IsJuliaMatrixRep and IsMutable, IsPosInt and IsSmallIntRep,
                       IsPosInt and IsSmallIntRep, IsObject ],
    function( m, i, j, val )
    local mat;
        mat := m!.m;
        mat[i,j] := val;
    end );
    
InstallOtherMethod( \*, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        return MakeJuliaMatrixRep(m1!.m * m2!.m);
    end );
    
InstallOtherMethod( \/, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        MakeJuliaMatrixRep(m1!.m * Inverse(m2!.m));
    end);
    
InstallMethod( InverseMutable, [IsJuliaMatrixRep], m -> MakeJuliaMatrixRep(InverseOp(m!.m)));
InstallMethod( InverseImmutable, [IsJuliaMatrixRep], m -> MakeImmutable(MakeJuliaMatrixRep(InverseOp(m!.m))));

InstallOtherMethod( AdditiveInverseMutable, [IsJuliaMatrixRep], m -> MakeJuliaMatrixRep(AdditiveInverseOp(m!.m)));
InstallOtherMethod( AdditiveInverseImmutable, [IsJuliaMatrixRep], m -> MakeImmutable(MakeJuliaMatrixRep(AdditiveInverseOp(m!.m))));

InstallOtherMethod(OneMutable, [IsJuliaMatrixRep], m -> MakeJuliaMatrixRep(OneOp(m!.m)));
InstallOtherMethod(OneImmutable, [IsJuliaMatrixRep], m -> MakeImmutable(MakeJuliaMatrixRep(OneOp(m!.m))));

InstallOtherMethod(ZeroMutable, [IsJuliaMatrixRep], m -> MakeJuliaMatrixRep(ZeroOp(m!.m)));
InstallOtherMethod(ZeroImmutable, [IsJuliaMatrixRep], m -> MakeImmutable(MakeJuliaMatrixRep(ZeroOp(m!.m))));

InstallOtherMethod(\+, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        return MakeJuliaMatrixRep(m1!.m + m2!.m);
    end);
    
InstallOtherMethod(\-, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        return MakeJuliaMatrixRep(m1!.m - m2!.m);
    end);

InstallOtherMethod(\^, [IsJuliaMatrixRep, IsInt and IsSmallIntRep], function( m1, k )
        if k < 0 then
            return MakeJuliaMatrixRep( (Inverse(m1)!.m)^(-k) );
        else
            return MakeJuliaMatrixRep( (m1!.m)^k );
        fi;
    end);
    
InstallOtherMethod(Display, [IsJuliaMatrixRep], function( m )
    local i,j;
        Print("[ ");
        for i in [1..NumberRows(m)] do
            Print("[ ");
            for j in [1..NumberColumns(m)] do
                Print(m[i,j]);
                if (j < NumberColumns(m)) then
                    Print(", ");
                fi;
            od;
            if (i < NumberRows(m)) then
                Print(" ] \n");
                Print("  ");
            else
                Print(" ] ] \n");
            fi;
        od;
    end);

InstallOtherMethod(Unpack, [IsJuliaMatrixRep], function(m)
    local v, i, j;
    v := EmptyPlist(NumberRows(m)*NumberColumns(m));
    for i in [1..NumberRows(m)] do
        for j in [1..NumberColumns(m)] do
            v[j+NumberColumns(m)*(i-1)] := m[i,j];
        od;
    od;
    return v;
end);


# Method: SetNiceMorphismForJuliaMatrixRepGroup
# Input: Matrix group G
# Output: -
# This function sets a nice morphism for G. The morphism is defined from the Julia matrix to the corresponding GAP matrix

DeclareOperation("SetNiceMorphismForJuliaMatrixRepGroup", [IsGroup]);

InstallMethod(SetNiceMorphismForJuliaMatrixRepGroup, [IsGroup], function(G)
    local hom, f, f_inv, ele, GAPGenerators, i, gens, GAPGroup, JuliaGAPMap;
    
        gens := GeneratorsOfGroup(G);
        Assert(0, ForAll(gens, IsJuliaMatrixRep));
        
        ele := gens[1];
        hom := Oscar_jl.iso_oscar_gap(Julia.base_ring(ele!.m));

        GAPGenerators := List(gens, g -> Oscar_jl.map_entries(hom, g!.m));
        
        GAPGroup := GroupByGenerators(GAPGenerators);
        
        f := m -> Oscar_jl.map_entries(hom,m!.m);
        f_inv := m -> MakeJuliaMatrixRep(Oscar_jl.preimage_matrix(hom,m));

        JuliaGAPMap := GroupHomomorphismByFunction(G,GAPGroup,f,f_inv);
        SetNiceMonomorphism(G,JuliaGAPMap);
        SetIsHandledByNiceMonomorphism(G, true);
    end);
    
InstallOtherMethod(\=, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        return m1!.m = m2!.m;
    end);
    
InstallOtherMethod(\<, [IsJuliaMatrixRep, IsJuliaMatrixRep], function( m1, m2 )
        Error("comparing IsJuliaMatrixRep not supported");
    end);

InstallOtherMethod(MinimalPolynomial, [IsJuliaMatrixRep], m -> Oscar_jl.minpoly(m!.m));
InstallOtherMethod(CharacteristicPolynomial, [IsJuliaMatrixRep], m -> Oscar_jl.charpoly(m!.m));
#InstallOtherMethod(IsSymmetric, [IsJuliaMatrixRep], m -> Oscar_jl.is_symmetric(m!.m));
InstallOtherMethod(IsMonomialMatrix, [IsJuliaMatrixRep], m -> Oscar_jl.is_monomial(m!.m));
InstallOtherMethod(RankMat, [IsJuliaMatrixRep], m -> Oscar_jl.rank(m!.m));
InstallOtherMethod(TraceMat, [IsJuliaMatrixRep], m -> Oscar_jl.tr(m!.m));

BindGlobal("TransformPolynomialFromJuliaToGAP", function(pol)
    local x, hom, res, i;

        hom := Oscar_jl.iso_oscar_gap(Julia.base_ring(pol.parent));
        x := Indeterminate(Julia.codomain(hom),"x");
        
        res := Zero(hom.header.codomain);
        for i in [0..(pol.length-1)] do
            res := res + x^i * hom(Julia.coeff(pol,i));
        od;
        
        return res;
    end);

# The following code is a workaround for the problem described at
# https://github.com/gap-system/gap/issues/5422
DeclareFilter( "IsPreImagesByAction" );

InstallMethod( PreImagesRepresentative,
    FamRangeEqFamElm,
    [ IsBijective and IsSPGeneralMapping and IsPreImagesByAction,
      IsMultiplicativeElementWithInverse ], 2000,
    function( map, elm )
    local R, gens_imgs, nice, intermed, comp;

    if not IsBound( map!.nice_inverse ) then
      R:= Range( map );
      Assert( 1, IsHandledByNiceMonomorphism( R ) );
      gens_imgs:= MappingGeneratorsImages( map );
      nice:= NiceMonomorphism( R );
      intermed:= List( gens_imgs[2], x -> ImagesRepresentative( nice, x ) );
      comp:= CompositionMapping(
                 GroupHomomorphismByImagesNC( NiceObject( R ), Source( map ),
                                              intermed, gens_imgs[1] ),
                 nice );
      map!.nice_inverse:= comp;
    fi;

    return ImagesRepresentative( map!.nice_inverse, elm );
    end );
