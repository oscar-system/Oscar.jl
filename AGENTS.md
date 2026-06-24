# Information for AI agents

## Module Organization
- `src/` - Main Julia sources. Subdirectories are organized by mathematical
  area or integration layer, for example `AlgebraicGeometry/`, `Groups/`,
  `NumberTheory/`, `PolyhedralGeometry/`, `Rings/`, `GAP/`, and `Polymake/`.
- `experimental/` - Julia sources for experimental components whose APIs are
  still allowed to change. Most subdirectories correspond to one experimental
  feature or package-style component.
- `gap/` - GAP-language sources and bundled GAP package code used by OSCAR's
  GAP integration.
- `test/` - Julia test suite. The directory layout mostly mirrors `src/`;
  add new tests to an existing file unless explicitly instructed otherwise.
- `docs/` - Documentation and user manual sources, including developer
  documentation under `docs/src/DeveloperDocumentation/`.
- `data/` - Static data files used by OSCAR functionality and tests, for
  example databases of solids and surfaces.
- `examples/` - Example scripts and demonstrations.
- `etc/` - Auxiliary developer and maintenance scripts.
- `Project.toml` and `Manifest.toml` - Julia package environment metadata.

## For all changes

- Add the AI tool as a Git co-author on all commits created by that tool.
- Whenever a pull request is opened, you MUST disclose that the pull request was written with the assistance of generative AI.

## Specific instructions for particular changes

### Doctests

#### Writing doctests

If you are asked to write new doctests, first review
`docs/src/DeveloperDocumentation/documentation.md` for best practices and for
the documented `build_doc`, `Oscar.doctest`, and `Oscar.doctest_fix` helpers.

#### Verifying doctests
If you have changed any `jldoctest` code blocks you should take
the following steps to verify your work:
- Review `docs/src/DeveloperDocumentation/documentation.md`. In particular,
  determine if any of the changed doctests require filters, labels, setup code,
  `build_doc`, or one of the targeted `Oscar.doctest` helpers.
- Run the doctests to verify that your change works:
    - To run all doctests: `julia --project=. -e 'using Oscar; Oscar.doctest()'`.
    - To run doctests for one function:
      `julia --project=. -e 'using Oscar; Oscar.doctest(function_name)'`.
    - To run doctests for files matching a path substring:
      `julia --project=. -e 'using Oscar; Oscar.doctest("/Rings/")'`.
    - IMPORTANT: The doctests may take up to 15 minutes. Do NOT terminate the
      doctests before completion. Do NOT use a timeout for doctests.
    - If you are ChatGPT, you may have to increase `yield_time_ms`.

Follow these steps for EVERY change you make in a doctest.

### Test changes

- If you are adding a new test, add it to an existing test file. Do not create
  a new test file unless explicitly instructed.
- Write one comment at the top of the test to explain what is being tested.
  Otherwise keep comments minimal.

### Writing code

After writing code, look up the docstring for each function you used. If there
are recommendations or additional considerations that apply to these functions,
make sure to take them into account.

Honor the style guide in `docs/src/DeveloperDocumentation/styleguide.md` for
code formatting, naming conventions, use of Unicode, and related conventions.

## Commit messages and pull requests

When writing commit messages, the title ("Brief summary") should not start with
a "component: ..." prefix. Instead it typically is a verb ("Add/Fix/Change/...").
In the body of the commit message, provide a brief prose summary
of the purpose of the changes made. Do not specifically mention added tests, comments,
documentation, etc., unless this is the main purpose of the change. Do not mention
the test plan, unless it differs from what you were instructed to do in AGENTS.md.
If your change fixes one or more issues, use the syntax "Fixes #" at the end of the
commit message, but do not include it in the title.

When referencing external GitHub PRs or issues, use proper GitHub interlinking
format, for example `owner/repo#123` for PRs/issues. When fixing CI failures,
include the link to the specific CI failure in the commit message.

When creating pull requests:
1. If the pull request consists of one commit only, use the body of the commit for the body of the pull request.
2. If there are multiple commits in the pull request, follow the same guidelines for the pull request as for the commit body.
3. Make sure that the base commit of the pull request is recent (within the past two days) - if not rebase your changes first.
4. You MUST disclose that the pull request was written with the assistance of generative AI.
