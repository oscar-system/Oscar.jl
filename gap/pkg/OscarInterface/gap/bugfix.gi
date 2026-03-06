##  `IsomorphismPermGroupForMatrixGroup` has been fixed in
##  https://github.com/gap-system/gap/pull/5339,
##  the function below makes this fix available to OSCAR.
##
##  Additionally, this function tries to improve the behaviour
##  by avoiding `NicomorphismOfGeneralMatrixGroup` in the situation
##  that a `NiceMonomorphism` was already stored
##  whose image is *not* a permutation group.
##  In this situation, we compute an `IsomorphismPermGroup` for the
##  image of the stored map, and return the composition of the two maps.
##
##  So the idea is that we trust GAP that the stored `NiceMonomorphism`
##  of `G` is really nice, in the sense that its image is better than `G`.
##  A situation where this holds (in particular in OSCAR) is
##  that a matrix group of cyclotomics knows a faithful reduction to
##  a matrix group over a finite field.
##  `NicomorphismOfGeneralMatrixGroup` is typically very expensive
##  in this situation.
BindGlobal( "IsomorphismPermGroupForMatrixGroup", function( G )
  local map;

  if HasNiceMonomorphism( G ) then
    map:= NiceMonomorphism( G );
    if not IsPermGroup( Range( map ) ) then
      # Trust GAP that this map is still useful.
      map:= CompositionMapping( IsomorphismPermGroup( Image( map ) ), map );
    fi;
  else
    # We cannot be sure about `IsHandledByNiceMonomorphism` for `G`.
    if not HasIsFinite( G ) then
      Info( InfoWarning,1,
            "IsomorphismPermGroup: The group is not known to be finite" );
    fi;
    map:= NicomorphismOfGeneralMatrixGroup( G, false, false );
    SetNiceMonomorphism( G, map );
  fi;
  # Now `G` stores a `NiceMonomorphism`.
  if IsPermGroup( Range( NiceMonomorphism( G ) ) ) then
    # Set an attribute if available.
    if IsAttribute( RestrictedNiceMonomorphism ) then
      map:= RestrictedNiceMonomorphism( G );
    else
      map:= RestrictedNiceMonomorphism( NiceMonomorphism( G ), G );
    fi;
  fi;
  if IsIdenticalObj( Source( map ), G ) and IsSurjective( map ) then
    return map;
  fi;
  map:= GeneralRestrictedMapping( map, G, ImagesSet( map, G ) );
  SetIsBijective( map, true );

  return map;
end );

InstallMethod( IsomorphismPermGroup,
    "matrix group",
    [ IsMatrixGroup ],
    IsomorphismPermGroupForMatrixGroup );

InstallMethod( IsomorphismPermGroup,
    "finite matrix group",
    [ IsMatrixGroup and IsFinite and IsHandledByNiceMonomorphism ],
    # We do not want the upranking via 'IsHandledByNiceMonomorphism',
    # analogous to the situation with the method for
    # 'IsGroup and IsFinite and IsHandledByNiceMonomorphism'in
    # 'lib/grpnice.gi'.
    [ [ IsMatrixGroup and IsFinite ], 1 ],
    IsomorphismPermGroupForMatrixGroup );
