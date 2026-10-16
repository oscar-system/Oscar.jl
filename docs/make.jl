#
# Build and deploy the manual. Run from the root of the repository as
#
#     julia --project=. docs/make.jl
#
# For a local build use `Oscar.@build_doc` instead.
#
using Oscar

# Put Documenter and the OscarDocs driver on the load path; Oscar itself comes
# from the active project.
Oscar.docs_env()

using Documenter, OscarDocs

OscarDocs.build(Oscar; warnonly=false, local_build=false, doctest=false, open_browser=false)

should_push_preview = true
if get(ENV, "GITHUB_ACTOR", "") == "dependabot[bot]"
  # skip preview for dependabot PRs as they fail due to lack of permissions
  should_push_preview = false
end

deploydocs(
   repo   = "github.com/oscar-system/Oscar.jl.git",
   deploy_repo = "github.com/oscar-system/OscarDocumentation",
   push_preview = should_push_preview,
   forcepush = true,
)
