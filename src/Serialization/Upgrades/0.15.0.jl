################################################################################
# Upgrade Summary
# This upgrade currentlz only adapts to type renamings.

push!(upgrade_scripts_set, UpgradeScript(
  v"0.15.0",
  function upgrade_0_15_0(s::UpgradeState, dict::Dict)

    upgraded_dict = dict

    renamings = Dict{String,String}([
      ("FlintPadicField", "PadicField"),
      ("padic", "PadicFieldElem"),
      ("arb", "ArbFieldElem"),
      ("acb", "AcbFieldElem"),
      ("AnticNumberField", "AbsSimpleNumField"),
      ("nf_elem", "AbsSimpleNumFieldElem"),
      ("NfAbsNS", "AbsNonSimpleNumField"),
      ("NfAbsNSElem", "AbsNonSimpleNumFieldElem"),
      ("Hecke.NfRel", "Hecke.RelSimpleNumField"),
      ("Hecke.NfRelElem", "Hecke.RelSimpleNumFieldElem"),
      ("Hecke.NfRelNS", "Hecke.RelNonSimpleNumField"),
      ("Hecke.NfRelNSElem", "Hecke.RelNonSimpleNumFieldElem"),
      ("Hecke.NumFieldEmbNfAbs", "Hecke.AbsSimpleNumFieldEmbedding"),
      ("Hecke.NumFieldEmbNfAbsNS", "Hecke.AbsNonSimpleNumFieldEmbedding"),
      ("Hecke.NumFieldEmbNfRel", "Hecke.RelSimpleNumFieldEmbedding"),
      ("Hecke.NumFieldEmbNfNS", "Hecke.RelNonSimpleNumFieldEmbedding"),
      ("GrpAbFinGen", "FinGenAbGroup"),
      ("GrpAbFinGenElem", "FinGenAbGroupElem"),
      ("Hecke.EmbeddedField", "EmbeddedNumField"),
      ("EmbeddedElem", "EmbeddedNumFieldElem"),
    ])

    function upgrade_type(d::String)
      return get(renamings, d, d)
    end
    function upgrade_type(d::Dict)
      upg_d = d
      upg_d[:name] = get(renamings, d[:name], d[:name])
      if d[:params] isa Dict && haskey(d[:params], :_type)
        upg_d[:params][:_type] = upgrade_type(d[:params][:_type])
      elseif d[:params] isa Vector
        upg_d[:params] = [upgrade_type(v) for v in d[:params]]
      end
      return upg_d
    end

    if haskey(dict, :_type)
      upgraded_dict[:_type] = upgrade_type(dict[:_type])
    end
        
    if haskey(dict, :data) && dict[:data] isa Dict
      upgraded_dict[:data] = upgrade_0_15_0(s, dict[:data])
    end

    if haskey(dict, :_refs)
      upgraded_refs = Dict()
      for (k, v) in dict[:_refs]
        upgraded_refs[k] = upgrade_0_15_0(s, v)
      end
      upgraded_dict[:_refs] = upgraded_refs
    end

    return upgraded_dict
  end
))
