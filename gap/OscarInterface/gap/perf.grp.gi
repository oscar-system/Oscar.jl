BindGlobal( "PerfGrpLoadOscar", function(sz)
  local p,pos,name,libname;
  if PERFRec=fail then
    ReadGrp("perf0.grp");
  fi;
  # get the index
  pos:=PositionSet(PERFRec.sizes,sz);
  if pos=fail then
    return fail;
  fi;
  if IsBound(PERFGRP[pos]) and PERFGRP[pos]<>fail and PERFGRP[pos][1]<>fail then
    return pos;
  fi;
  # get the file number
  p:=Position(PERFRec.newlyAdded,sz);
  if p=fail then
    p:=PositionSorted(PERFRec.covered,pos);
  else
    p:=12+p;
  fi;
  name:=Concatenation("perf",String(p),".grp");
  libname := SHALLOW_COPY_OBJ( "grp" );
  APPEND_LIST_INTR( libname, "/" );
  APPEND_LIST_INTR( libname, name );
  if not READ_GAP_ROOT( libname ) then
    if p=27 or p=33 then
        Oscar._load_hulpke_extraperfect();
    fi;
  fi;
  if not READ_GAP_ROOT( libname ) then
    Error("\n\n",
    "Unfortunately, the file ",name," is not available in the `grp` subdirectory ",
    "of your GAP installation.\n\n\n");
  fi;

  ReadGrp(name);
  return pos;
end );

ReplaceGapFunc("PerfGrpLoad", PerfGrpLoadOscar);
