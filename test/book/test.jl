using Oscar

using DeepDiffs
using REPL
using Pkg
using Test
isdefined(Main, :FakeTerminals) || include(joinpath(pkgdir(REPL),"test","FakeTerminals.jl"))

@testset "OscarBookExamples" verbose=true begin

  broken = [
            # non-stable:
            "specialized/joswig-kastner-lorenz-confirmable-workflows/versioninfo.jlcon",

            # backtrace err
            "introduction/introduction/julia.jlcon",
            "introduction/introduction/julia2.jlcon",
            "introduction/introduction/julia3.jlcon",
           ]
  skipped = [
             # Something broken, probably needs some updated GAP packages
             "specialized/breuer-nebe-parker-orthogonal-discriminants/expl_syl.jlcon",

             # these are skipped because they slow down the tests too much:

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


  # this is run on both sample and output
  function normalize_repl_output(s::AbstractString)
    result = string(s)
    lafter = length(result)
    lbefore = lafter+1
    while lafter < lbefore
      lbefore = lafter
      result = replace(result, r"^       "m => "")
      result = strip(result)
      result = replace(result, r"julia>$"s => "")
      # canonicalize numbered anonymous functions
      result = replace(result, r"^\s*(?:#\d+)?(.* \(generic function with ).*\n"m => s"\1\n")
      # remove timings
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

  # this is only applied to the output
  function sanitize_output(s::AbstractString)
    result = s
    lafter = length(result)
    lbefore = lafter+1
    while lafter < lbefore
      lbefore = lafter
      # remove control characters:
      #   kill line including previous line (move_up)
      result = replace(result, r"^.*\n.*?\e\[1A\r\e\[0K"m=>"")
      #   kill line
      result = replace(result, r"^.*(?<!\e\[1A)\r\e\[0K"m=>"")
      #   move right
      result = replace(result, r"\r\e\[\d+C"=>"")
      #   empty line with \r
      result = replace(result, r"^\r\r\n"m=>"")

      # remove mockdule from prompt-prefix
      result = replace(result, r"\(Main\.__\d+\) julia>"m => "julia>")

      # replace weird duplicated non-comments
      result = replace(result, r"^(#.*)\n\1$"m => s"\1")

      # remove end marker
      result = replace(result, "\n\njulia> println(\"\\nEND_BLOCK\");" => "")
      lafter = length(result)
    end
    result = normalize_repl_output(result)
    return result
  end

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
      # make it accessible from Main
      setproperty!(Main, sym, mockdule)
      Core.eval(mockdule, :(eval(x) = Core.eval($(mockdule), x)))
      Core.eval(mockdule, :(include(x) = Base.include($(mockdule), abspath(x))))

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
        REPL.run_repl(repl)
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
    # add overlay project for plots
    custom_load_path = []
    copy!(custom_load_path, Base.DEFAULT_LOAD_PATH)
    curdir = pwd()
    act_proj = dirname(Base.active_project())
    osc_proj = dirname(Base.identify_package_env("Oscar")[2])
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

            # clear verbosity levels before each chapter
            empty!(Hecke.VERBOSE_LOOKUP)

            copy!(LOAD_PATH, custom_load_path)
            auxmain = joinpath(Oscar.oscardir, "test/book", chapter, "auxiliary_code", "main.jl")
            if isfile(auxmain)
              # add overlay project for aux file
              # and run it from temp dir
              temp = mktempdir()
              Pkg.activate(temp; io=devnull)
              pushfirst!(LOAD_PATH, "$osc_proj")
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
                res = @test normalize_repl_output(content) == computed broken=(full_file in broken)
                if res isa Test.Fail
                  println(deepdiff(normalize_repl_output(content),computed))
                end
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
            end
            println("  closing mockrepl: $(mockrepl.mockdule)")
            close_repl(mockrepl)
          end
        end
      end
    finally
      # restore some state
      Main.REPL.activate(Main)
      Pkg.activate("$act_proj"; io=devnull)
      cd(curdir)
      copy!(LOAD_PATH, Base.DEFAULT_LOAD_PATH)
    end
    nothing
  end

  # test all chapters
  test_chapter()

end
