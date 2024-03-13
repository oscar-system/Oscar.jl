const excluded = [
                  # Excluded due to performance issues
                  "markwig-ristau-schleis-faithful-tropicalization/eliminate_xz",
                  "markwig-ristau-schleis-faithful-tropicalization/eliminate_yz",
                  "number-theory/cohenlenstra.jlcon",
                  "bies-turner-string-theory-applications/SU5.jlcon",
                  "brandhorst-zach-fibration-hopping/vinberg_2.jlcon",
                  "brandhorst-zach-fibration-hopping/vinberg_3.jlcon",
                  "cornerstones/polyhedral-geometry/ch-benchmark.jlcon",
                  "number-theory/unit_plot.jlcon",                                            
                  "number-theory/intro_plot_lattice.jlcon",                                   
                  "number-theory/intro5_0.jlcon",                                             
                  # Excluded for easier testing
                  "specialized",
                  "number-theory",
                  "groups",
                  "algebraic-geometry",
                  "introduction",
                 ]  

using REPL
include(joinpath(pkgdir(REPL),"test","FakeTerminals.jl"))


function sanitize_output(s::AbstractString)
  result = s
  println("length before: ", length(result))
  lafter = length(result)
  lbefore = lafter+1
  while lafter < lbefore
    lbefore = lafter
    result = replace(result, r"\r\e\[\d+[A-Z]"=>"")
    result = replace(result, r"\e\[\?\d+[a-z]"=>"")
    result = replace(result, r"\r\r\n"=>"")
    result = replace(result, r"julia> julia> julia>" => "julia>")
    lafter = length(result)
  end
  println("length after: ", length(result))
  return result[1:length(result)-32]
end

function run_repl_string(s::AbstractString; jlcon_mode=true)
  input_string = s
  if jlcon_mode
    input_string = "\e[200~$s\e[201~\n"
  end
  input = Pipe()
  err = Pipe()
  Base.link_pipe!(input, reader_supports_async=true, writer_supports_async=true)
  outputbuf = IOBuffer()
  Base.link_pipe!(err, reader_supports_async=true, writer_supports_async=true)
  stdin_write = input.in
  options = REPL.Options(confirm_exit=false, hascolor=false)
  # options.extra_keymap = REPL.LineEdit.escape_defaults
  repl = REPL.LineEditREPL(FakeTerminals.FakeTerminal(input.out, outputbuf, err.in, options.hascolor), options.hascolor, false)
  repl.options = options
  # repl.specialdisplay = REPL.REPLDisplay(repl)
  repltask = @async begin
    REPL.run_repl(repl)
  end
  input_task = @async write(stdin_write, input_string)
  wait(input_task)
  input_task = @async write(stdin_write, "\x04")
  wait(input_task)
  wait(repltask)
  result = String(take!(outputbuf))
  return sanitize_output(result)
end


oefile = joinpath(Oscar.oscardir, "test/book/ordered_examples.json")
ordered_examples = load(oefile)
for (chapter, example_list) in ordered_examples
  println(chapter)
  for example in example_list
    full_file = joinpath(chapter, example)
    exclude = filter(s->occursin(s, full_file), excluded)
    if length(exclude) == 0
      println("   "*example)
      if occursin("jlcon", example)
        content = read(joinpath(Oscar.oscardir, "test/book", full_file), String)
        computed = run_repl_string(content)
        @test content == computed
      elseif occursin("jl", example)
        content = read(joinpath(Oscar.oscardir, "test/book", full_file), String)
        run_repl_string(content; jlcon_mode=false)
      end
    # else
    #   println("   "*example*" excluded")
    end
  end
end

