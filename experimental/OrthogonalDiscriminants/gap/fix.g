###############################################################################
##
##  This file contains bugfixes and improvements
##  (They are contained in the master branch of the GAP development
##  and will be contained in GAP from version 4.13 on.)
##


###############################################################################
##
##  Improve the behavior of 'Indicator' in characteristic 2.
##
InstallMethod( IndicatorOp,
    "for a Brauer character table and <n> = 2",
    [ IsBrauerTable, IsHomogeneousList, IsPosInt ], 10,
    function( modtbl, ibr, n )
    local ind, princ, i, real, ordtbl, ordindicator, dec, j;

    if UnderlyingCharacteristic( modtbl ) <> 2 or ibr <> Irr( modtbl ) then
      TryNextMethod();
    fi;

    ind:= [];
    princ:= BlocksInfo( modtbl )[1].modchars;

    for i in [ 1 .. Length( ibr ) ] do
      if ibr[i] <> ComplexConjugate( ibr[i] ) then
        # Non-real characters have indicator 0.
        ind[i]:= 0;
      elif not i in princ then
        # Real characters outside the principal block have indicator 1.
        ind[i]:= 1;
      elif Set( ibr[i] ) = [ 1 ] then
        # The trivial character is defined to have indicator 1.
        ind[i]:= 1;
      else
        # Set 'Unknown()' for all other characters.
        ind[i]:= Unknown();
      fi;
    od;

    # We use a criterion described in [GowWillems95, Lemma 1.2].
    # The 2-modular restriction of an orthogonal (indicator '+')
    # ordinary character preserves a quadratic form.
    # If the trivial character is not a constituent of this restriction
    # then any real constituent of this restriction that occurs with odd
    # multiplicity has indicator +.
    real:= PositionsProperty( ibr, x -> x = ComplexConjugate( x ) );
    if ForAny( ind, IsUnknown ) then
      ordtbl:= OrdinaryCharacterTable( modtbl );
      ordindicator:= Indicator( ordtbl, 2 );
      dec:= DecompositionMatrix( modtbl );
      for i in [ 1 .. Length( dec ) ] do
        if ordindicator[i] = 1 and dec[i][1] = 0 then
          for j in Filtered( real, x -> IsOddInt( dec[i, x] ) ) do
            if IsUnknown( ind[j] ) then
              ind[j]:= 1;
            else
              Assert( 1, ind[j] = 1 );
            fi;
          od;
        fi;
      od;
    fi;

    return ind;
    end );
