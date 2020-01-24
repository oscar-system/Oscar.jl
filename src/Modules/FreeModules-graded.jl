mutable struct FreeModule_dec <: Module
  d::Array{GrpAbFinGenElem, 1}
  R::MPolyRing_dec
end

function FreeModule(R::MPolyRing_dec, n::Int; cached::Bool = false)
  return FreeModule_dec([R.D[0] for i=1:n], R)
end

function FreeModule(R::MPolyRing_dec, d::Array{GrpAbFinGenElem, 1})
  return FreeModule_dec(d, R)
end

function show(io::IO, F::FreeModule_dec)
  @show_name(io, F)
  @show_special(io, F)
  print(io, "Free module of rank $(length(F.d)) over ")
  print(IOContext(io, :compact =>true), F.R)
  if isgraded(F.R)
    print(io, ", graded by\n")
  else
    print(io, ", filtrated by\n")
  end
  for i=1:dim(F)
  end
end


dim(F::FreeModule_dec)  = length(F.d)
ngens(F::FreeModule_dec) = dim(F)


struct FreeModuleElem_dec <: ModuleElem
  r::SRow{MPolyRing_dec}
  parent::FreeModule_dec
end

function basis(F::FreeModule_dec)
  bas = FreeModuleElem_dec[]
  for i=1:dim(F)
    s = SRow(F.R)
    push!(s, [i, W(1)])
end

