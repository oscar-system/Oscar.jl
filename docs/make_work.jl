using Documenter, Oscar, DocumenterMarkdown, DocumenterCitations

import Documenter:
    Anchors,
    DocTests,
    Documents,
    Documenter,
    Utilities

import Documenter.Utilities: Selectors

using Documenter.Builder
import Documenter.Builder: SetupBuildDirectory, walk_navpages
import Documenter.Documents: Document, Page

#=
 copied from Documenter (0.26.1)

 Adds the loop to collect the AA/Nemo/Hecke doc files into the
 build tree of Oscar.

=#
function Selectors.runner(::Type{SetupBuildDirectory}, doc::Documents.Document)
    @info "SetupBuildDirectory: setting up build directory."

    # Frequently used fields.
    build  = doc.user.build
    source = doc.user.source
    workdir = doc.user.workdir

    # The .user.source directory must exist.
    isdir(source) || error("source directory '$(abspath(source))' is missing.")

    # We create the .user.build directory.
    # If .user.clean is set, we first clean the existing directory.
    doc.user.clean && isdir(build) && rm(build; recursive = true)
    isdir(build) || mkpath(build)

    # We'll walk over all the files in the .user.source directory.
    # The directory structure is copied over to .user.build. All files, with
    # the exception of markdown files (identified by the extension) are copied
    # over as well, since they're assumed to be images, data files etc.
    # Markdown files, however, get added to the document and also stored into
    # `mdpages`, to be used later.
    mdpages = String[]
    for (root, dirs, files) in walkdir(source)
        for dir in dirs
            d = normpath(joinpath(build, relpath(root, source), dir))
            isdir(d) || mkdir(d)
        end
        for file in files
            src = normpath(joinpath(root, file))
            dst = normpath(joinpath(build, relpath(root, source), file))

            if workdir == :build
                # set working directory to be the same as `build`
                wd = normpath(joinpath(build, relpath(root, source)))
            elseif workdir isa Symbol
                # Maybe allow `:src` and `:root` as well?
                throw(ArgumentError("Unrecognized working directory option '$workdir'"))
            else
                wd = normpath(joinpath(doc.user.root, workdir))
            end

            if endswith(file, ".md")
                push!(mdpages, Utilities.srcpath(source, root, file))
                Documents.addpage!(doc, src, dst, wd)
            else
                cp(src, dst; force = true)
            end
        end
    end

    #= essentially a copy of above, but reading the extra repositories
       and normalizing the paths suitably.
       only `.md` files are copied. There is/was a clash in the other
       files (icons, graphics, ...)
    =#
    for extra in ["Hecke", "Nemo", "AbstractAlgebra"]
        for (root, dirs, files) in walkdir(normpath(joinpath(dirname(pathof(getproperty(Main, Symbol(extra)))), "..", "docs", "src")))
            for dir in dirs
                d = normpath(joinpath(build, relpath(root, source), dir))
                isdir(d) || mkdir(d)
            end
            for file in files
                src = normpath(joinpath(root, file))
                dst = normpath(joinpath(build, relpath(root, source), file))

                if workdir == :build
                    # set working directory to be the same as `build`
                    wd = normpath(joinpath(build, relpath(root, source)))
                elseif workdir isa Symbol
                    # Maybe allow `:src` and `:root` as well?
                    throw(ArgumentError("Unrecognized working directory option '$workdir'"))
                else
                    wd = normpath(joinpath(doc.user.root, workdir))
                end

                if endswith(file, ".md")
                    push!(mdpages, Utilities.srcpath(source, root, file))
                    Documents.addpage!(doc, src, dst, wd, extra)
                end
            end
        end
    end


    # If the user hasn't specified the page list, then we'll just default to a
    # flat list of all the markdown files we found, sorted by the filesystem
    # path (it will group them by subdirectory, among others).
    userpages = isempty(doc.user.pages) ? sort(mdpages, lt=lt_page) : doc.user.pages

    # Populating the .navtree and .navlist.
    # We need the for loop because we can't assign to the fields of the immutable
    # doc.internal.
    for navnode in walk_navpages(userpages, nothing, doc)
        push!(doc.internal.navtree, navnode)
    end

    # Finally we populate the .next and .prev fields of the navnodes that point
    # to actual pages.
    local prev::Union{Documents.NavNode, Nothing} = nothing
    for navnode in doc.internal.navlist
        navnode.prev = prev
        if prev !== nothing
            prev.next = navnode
        end
        prev = navnode
    end
end

# new signature: trunc is added to sanitize paths (truncate)
function Documenter.Documents.addpage!(doc::Document, src::AbstractString, dst::AbstractString, wd::AbstractString, trunc::AbstractString)
    if occursin(trunc, dst)
      dst = sanitize_dst(dst, trunc)
    end
    page = Page(src, dst, wd)
    # page's identifier is the path relative to the `doc.user.source` directory
    name = normpath(relpath(src, doc.user.source))
    if occursin(trunc, name)
      name = sanitize(name, trunc)
    end
    doc.blueprint.pages[name] = page
end

# Overwrite printing to make the header not full of redundant nonsense
# Turns
#   Hecke.Order - Method
# into
#   Order - Method

# To remove the '-'
# Documenter.Utilities.print_signature(io::IO, signature)        = print(io, signature)

# To remove the "Method", "Type", "Module" use the following
# Documenter.Utilities.doccat(b::Base.Docs.Binding, ::Type)  = ""
# doccat(::Type)     = ""
# doccat(::Module)   = ""

# Remove the module prefix
Base.print(io::IO, b::Base.Docs.Binding) = print(io, b.var)

# Sanitize paths

function sanitize(a::AbstractString, n::AbstractString)
  b = splitpath(replace(a, Regex(".*$(n)") => n))
  return joinpath(b[1], b[3:end]...)
end

function sanitize_dst(a::AbstractString, n::AbstractString)
  b = splitpath(replace(a, Regex(".*$(n)") => n))
  return joinpath("build", b[1], b[3:end]...)
end

bla = normpath(joinpath(dirname(pathof(Hecke)), "..", "docs", "src"))
const hecke = sanitize(bla, "Hecke")

bla = normpath(joinpath(dirname(pathof(Nemo)), "..", "docs", "src"))
const nemo = sanitize(bla, "Nemo")

bla = normpath(joinpath(dirname(pathof(AbstractAlgebra)), "..", "docs", "src"))
const aa = sanitize(bla, "AbstractAlgebra")

bib = CitationBibliography(joinpath(@__DIR__, "oscar_references.bib"), sorting = :nyt)

function doit(strict::Bool = true, local_build::Bool = false)

  s = prod(eachline(joinpath(Oscar.oscardir, "docs", "doc.main")))
  s = replace(s, r"\$\(aa\)" => "$aa")
  s = replace(s, r"\$\(nemo\)" => "$nemo")
  s = replace(s, r"\$\(hecke\)" => "$hecke")


  doc = Meta.eval(Meta.parse(s))

  cd(joinpath(Oscar.oscardir, "docs")) do

  makedocs(bib,
         format   = Documenter.HTML(prettyurls = !local_build),
#         format   = Documenter.HTML(),
#         format   = Markdown(),
         sitename = "Oscar.jl",
         modules = [Oscar, Hecke, Nemo, AbstractAlgebra, Singular],
         clean = true,
         doctest = true,
         strict = strict,
         checkdocs = :none,
         pages    = doc)
end

end
