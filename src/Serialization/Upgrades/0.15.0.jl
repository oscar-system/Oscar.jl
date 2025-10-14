################################################################################
# Upgrade Summary
# This upgrade currentlz only adapts to type renamings.

push!(upgrade_scripts_set, UpgradeScript(
  v"0.15.0",
  function upgrade_0_15_0(s::UpgradeState, dict::AbstractDict{Symbol, Any})
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

    upgraded_dict = rename_types(dict, renamings)
    
    if haskey(dict, :data) && dict[:data] isa AbstractDict
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
