excluded = String[
                  # Excluded for easier testing
                  #"number-theory", # OK! + TODO check cohenlenstra
                  #"groups", # OK!
                  "algebraic-geometry",
                  #"polyhedral-geometry", # OK! but volume conflicts with mixedsubdiv from tropical
                  #"introduction", # OK!
                  #"specialized/boehm-breuer-git-fans", # OK!
                  #"specialized/markwig-ristau-schleis-faithful-tropicalization", # OK!
                  #"specialized/holt-ren-tropical-geometry", # OK!
                  #"specialized/bies-turner-string-theory-applications", # OK!
                  #"specialized/aga-boehm-hoffmann-markwig-traore", # OK!
                  #"specialized/joswig-kastner-lorenz-confirmable-workflows", # OK!
                  #"specialized/eder-mohr-ideal-theoretic", # OK!
                  #"specialized/kuehne-schroeter-matroids", # OK!
                  #"specialized/rose-sturmfels-telen-tropical-implicitization", # OK! + broken...
                  #"specialized/fang-fourier-monomial-bases", # OK!
                  #"specialized/brandhorst-zach-fibration-hopping", # OK! + TODO check long
                  #"specialized/breuer-nebe-parker-orthogonal-discriminants", # OK! + TODO waiting for fix
                  #"specialized/decker-schmitt-invariant-theory", # OK! + FIXME
                  #"specialized/bies-kastner-toric-geometry", # OK!
                 ]  

broken = [
                  # Something broken in Oscar
                  "specialized/breuer-nebe-parker-orthogonal-discriminants/expl_syl.jlcon",

                  # more control chars??
                  "specialized/boehm-breuer-git-fans/explG25_1.jlcon",

                  # weird errors with last line that has no output:
                  #"specialized/bies-turner-string-theory-applications/SU5.jlcon",
                  #"specialized/boehm-breuer-git-fans/explG25_8.jlcon",
                  #"specialized/aga-boehm-hoffmann-markwig-traore/graphname.jlcon",
                  #"specialized/rose-sturmfels-telen-tropical-implicitization/gen_impl.jlcon",
                  #"specialized/rose-sturmfels-telen-tropical-implicitization/hyperdet.jlcon",
                  #"specialized/rose-sturmfels-telen-tropical-implicitization/pol_from_surface.jlcon",
                  #"specialized/rose-sturmfels-telen-tropical-implicitization/chow_fan.jlcon",
                  #"specialized/rose-sturmfels-telen-tropical-implicitization/chow_transl.jlcon",

                  # TODO: need to fix column width for output?
                  "specialized/markwig-ristau-schleis-faithful-tropicalization/eliminate_xz.jlcon",

                  # non-stable:
                  "specialized/joswig-kastner-lorenz-confirmable-workflows/versioninfo.jlcon",

                  # output changed in oscar master? TODO check + adapt
                  "specialized/decker-schmitt-invariant-theory/gleason.jlcon",

                  # output changes?
                  "specialized/decker-schmitt-invariant-theory/cox_ring.jlcon",

                  # backtrace err
                  "introduction/introduction/julia.jlcon",
                  "introduction/introduction/julia2.jlcon",
                  "introduction/introduction/julia3.jlcon",
               ]
skipped = [
                  # sometimes very slow: 4000-30000s
                  "specialized/brandhorst-zach-fibration-hopping/vinberg_2.jlcon",
                  # very slow: 24000s
                  "cornerstones/number-theory/cohenlenstra.jlcon",
                  # ultra slow: time unknown
                  "specialized/markwig-ristau-schleis-faithful-tropicalization/eliminate_yz.jlcon",

                  # somewhat slow (~300s)
                  "cornerstones/polyhedral-geometry/ch-benchmark.jlcon",
                  "specialized/brandhorst-zach-fibration-hopping/vinberg_3.jlcon",
                ]

dispsize = (40, 130)

using Oscar

using REPL
using Random
using Pkg
import Random.Xoshiro
using Test
isdefined(Main, :FakeTerminals) || include(joinpath(pkgdir(REPL),"test","FakeTerminals.jl"))


function normalize_repl_output(s::AbstractString)
  result = string(s)
  lafter = length(result)
  lbefore = lafter+1
  while lafter < lbefore
    lbefore = lafter
    result = replace(result, r"^       "m => "")
    result = strip(result)
    result = replace(result, r"julia>$"s => "")
    result = replace(result, r"^\s*(?:#\d+)?(.* \(generic function with ).*\n"m => s"\1\n")
    result = replace(result, r"^\s*[0-9\.]+ seconds \(.* allocations: .*\)$"m => "<timing>\n")
    # this removes the package version slug, filename and linenumber
    result = replace(result, r" @ \w* ?~/\.julia/packages/(?:Nemo|Hecke|AbstractAlgebra|Polymake)/\K[\w\d]+/.*\.jl:\d+"m => "")
    lafter = length(result)
  end
  return strip(result)
end


function sanitize_input(s::AbstractString)
  result = s
  result = normalize_repl_output(result)
  return result
end

function sanitize_output(s::AbstractString)
  result = s
  # println("length before: ", length(result))
  lafter = length(result)
  lbefore = lafter+1
  while lafter < lbefore
    lbefore = lafter
    result = replace(result, r"\(Main\.__\d+\) julia>"m => "julia>")
    result = replace(result, r"\r\e\[\d+[A-Z]"=>"")
    result = replace(result, r"\e\[\?\d+[a-z]"=>"")
    result = replace(result, r"\r\r\n"=>"")
    # due to not interpreting control characters we need to remove duplicate prompts
    result = replace(result, r"(julia> ){2,}" => "julia> ")
    # replace weird duplicated comments
    result = replace(result, r"^(#.*)\n\1$"m => s"\1")
    # replace weird duplicated lines
    result = replace(result, r"^(julia> .*)\1$"m => s"\1")
    # remove end marker
    result = replace(result, "\n\njulia> println(\"\\nEND_BLOCK\");" => "")
    lafter = length(result)
  end
  result = replace(result, r"julia> visualize\(PC\)julia> visualize\(PC\)" => "julia> visualize(PC)")
  result = normalize_repl_output(result)
  # println("length after: ", length(result))
  return result
end

#function set_task_rngstate!(state::Xoshiro)
#  Random.setstate!(Random.default_rng(), state.s0, state.s1, state.s2, state.s3, state.s4)
#end
#function get_task_rngstate!(state::Xoshiro)
#  t = current_task()
#  Random.setstate!(state, t.rngState0, t.rngState1, t.rngState2, t.rngState3, t.rngState4)
#end

function get_preamble(modstr::String)
  return """
    REPL.activate($modstr);
    using Base.MainInclude: ans
    using Oscar;
    eval(Oscar.doctestsetup());
    println("\\nEND_PREAMBLE");
    """
end

struct MockREPLHelper
  mockdule::Module
  stdin_write::IO
  out_stream::IOContext
  input::Pipe
  output::Pipe
  err::Pipe
  repltask

  function MockREPLHelper(prefix::AbstractString)
    sym = Symbol("__", lstrip(string(gensym()), '#'))
    mockdule = Module(sym)
    # make it accessible
    setproperty!(Main, sym, mockdule)
    Core.eval(mockdule, :(eval(x) = Core.eval($(mockdule), x)))
    Core.eval(mockdule, :(include(x) = Base.include($(mockdule), abspath(x))))
    #global replmockdule[] = mockdule

    input = Pipe()
    output = Pipe()
    err = Pipe()
    Base.link_pipe!(input, reader_supports_async=true, writer_supports_async=true)
    Base.link_pipe!(output, reader_supports_async=true, writer_supports_async=true)
    Base.link_pipe!(err, reader_supports_async=true, writer_supports_async=true)
    stdin_write = input.in
    out_stream = IOContext(output.in, :displaysize=>dispsize)
    options = REPL.Options(confirm_exit=false, hascolor=false)
    repl = REPL.LineEditREPL(FakeTerminals.FakeTerminal(input.out, out_stream, err.in, options.hascolor), options.hascolor, false)
    repl.options = options
    Base.active_repl = repl
    repltask = @async begin
      #set_task_rngstate!(rng)
      REPL.run_repl(repl)
      #get_task_rngstate!(rng)
    end
    preout = redirect_stdout(out_stream) do
      preamble_task = @async write(stdin_write, get_preamble(string(mockdule)))
      wait(preamble_task)
      readuntil(output.out, "\nEND_PREAMBLE\n")
    end
    # crude error checking, but using `err` doesn't seem to work here
    if occursin("ERROR", preout)
      error("Error in preamble for $prefix:\n$preout")
    end

    return new(mockdule, stdin_write, out_stream, input, output, err, repltask)
  end
end

function close_repl(mockrepl::MockREPLHelper)
  input_task = @async write(mockrepl.stdin_write, "\x04")
  wait(input_task)
  wait(mockrepl.repltask)
  close(mockrepl.output)
end

function run_repl_string(mockrepl::MockREPLHelper, s::AbstractString; jlcon_mode=true)
  input_string = s
  if jlcon_mode
    input_string = "\e[200~$s\e[201~"
  end
  REPL.activate(mockrepl.mockdule)
  result = redirect_stdout(mockrepl.out_stream) do
    input_task = @async begin
                          write(mockrepl.stdin_write, input_string)
                          write(mockrepl.stdin_write, """\nprintln("\\nEND_BLOCK");\n""")
                        end
    wait(input_task)
    readuntil(mockrepl.output.out, "\nEND_BLOCK")
  end
  REPL.activate(Main)
  return sanitize_output(result)
end

function test_chapter(chapter::String="")
  println("Running OscarBook-Examples\n")

  # add overlay project for plots
  custom_load_path = []
  copy!(custom_load_path, Base.DEFAULT_LOAD_PATH)
  curdir = pwd()
  act_proj = dirname(Base.active_project())
  try
    plots = mktempdir()
    Pkg.activate(plots; io=devnull)
    Pkg.add("Plots"; io=devnull)
    Pkg.activate("$act_proj"; io=devnull)
    pushfirst!(custom_load_path, plots)

    oefile = joinpath(Oscar.oscardir, "test/book/ordered_examples.json")
    ordered_examples = load(oefile)
    if length(chapter) > 0
      ordered_examples = Dict("$chapter" => ordered_examples[chapter])
    end
    withenv("LINES" => dispsize[1], "COLUMNS" => dispsize[2], "DISPLAY" => "") do
      for (chapter, example_list) in ordered_examples
        cd(curdir)
        @testset "$chapter" verbose=true begin
          println(chapter)
          mockrepl = MockREPLHelper(chapter)
          println("  created mockrepl: $(mockrepl.mockdule)")

          copy!(LOAD_PATH, custom_load_path)
          auxmain = joinpath(Oscar.oscardir, "test/book", chapter, "auxiliary_code", "main.jl")
          if isfile(auxmain)
            # add overlay project for aux file
            # and run it from temp dir
            temp = mktempdir()
            Pkg.activate(temp; io=devnull)
            Pkg.develop(path=act_proj; io=devnull)
            pushfirst!(LOAD_PATH, "$act_proj")
            cp(auxmain,joinpath(temp, "main.jl"))
            cd(temp)
            run_repl_string(mockrepl, """include("$(joinpath(temp,"main.jl"))")\n""")
            LOAD_PATH[1] = temp
            Pkg.activate("$act_proj"; io=devnull)
            println("      done with aux")
          end

          run_repl_string(mockrepl, """Oscar.randseed!(42);""")
          for example in example_list
            full_file = joinpath(chapter, example)
            exclude = filter(s->occursin(s, full_file), excluded)
            if length(exclude) == 0
              println("    $example$(full_file in skipped ? ": skip" :
                                     full_file in broken ? ": broken" : "")")
              filetype = endswith(example, "jlcon") ? :jlcon :
                         endswith(example, "jl") ? :jl : :unknown
              content = read(joinpath(Oscar.oscardir, "test/book", full_file), String)
              if filetype == :jlcon && !occursin("julia> ", content)
                filetype = :jl
                @warn "possibly wrong file type: $full_file"
              end
              if full_file in skipped
                @test run_repl_string(mockrepl, content) isa AbstractString skip=true
              elseif filetype == :jlcon
                content = sanitize_input(content)
                computed = run_repl_string(mockrepl, content)
                @test normalize_repl_output(content) == computed broken=(full_file in broken)
              elseif filetype == :jl
                if occursin("# output\n", content)
                  (code, res) = split(content, "# output\n"; limit=2)
                  # TODO do we want to compare with `res` ?
                  @test run_repl_string(mockrepl, code; jlcon_mode=false) isa AbstractString broken=(full_file in broken)
                else
                  @test run_repl_string(mockrepl, content; jlcon_mode=false) isa AbstractString broken=(full_file in broken)
                end
              else
                @warn "unknown file type: $full_file"
              end
            # else
            #   println("   "*example*" excluded")
            end
          end
          println("  closing mockrepl: $(mockrepl.mockdule)")
          close_repl(mockrepl)
        end
      end
    end
  finally
    Main.REPL.activate(Main)
    Pkg.activate("$act_proj"; io=devnull)
    cd(curdir)
    copy!(LOAD_PATH, Base.DEFAULT_LOAD_PATH)
  end
  nothing
end
