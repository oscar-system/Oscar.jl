# A closure capturing a variable that is assigned more than once forces Julia
# to allocate a `Core.Box` for it, which defeats type inference at every use.
# See docs/src/DeveloperDocumentation/closure_boxes.md for the background and
# for how to avoid this. `Test.detect_closure_boxes` exists since Julia 1.14.
if isdefined(Test, :detect_closure_boxes)
  @testset "Closure boxes" begin
    boxes = Test.detect_closure_boxes(Oscar)
    for (m, vars) in boxes
      println("Boxed variable(s) ", join(vars, ", "), " in ", m)
    end
    @test length(boxes) == 0
  end
end
