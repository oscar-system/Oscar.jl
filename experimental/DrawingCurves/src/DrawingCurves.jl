
include("IsotopyGraph.jl")
include("backend_tikz.jl")



################################################################################
################################################################################
##
## msolve
## 
################################################################################
################################################################################

# Use real_solutions to compute roots of a univariate polynomial
function _real_roots(f; selected_precision::Int)
   (result,_) = Oscar.real_solutions(ideal(f); precision=selected_precision)
   # println("res: ",result)
   result = [r[1] for r in result]
   sort!(result)
   return result
end

function analyse_singularity(f, projy, xpt, ypt; selected_precision::Int)
   Rxy = parent(f)
   y = gens(Rxy)[2]
   x = gens(Rxy)[1]
   type = :singularity
   if(xpt[1] == xpt[2])
      # Exact case
      fpy = derivative(f, y)
      fpx = derivative(f, x)
      issingular = is_zero(Oscar.evaluate(fpy, [xpt[1], ypt[1]])) && is_zero(Oscar.evaluate(fpx, [xpt[1], ypt[1]]))
      if !issingular
         type = :ytangent
      end
      pts = _real_roots(projy(Oscar.evaluate(f, [Rxy(xpt[1]), y])); selected_precision)
      singindex = findall(p->p==ypt[1], pts)
      @assert length(singindex)==1 "Something went wrong for exact solution"
      return (pts, singindex[1], type)
   else
      xbefore = xpt[1]
      ptsbefore = _real_roots(projy(Oscar.evaluate(f, [Rxy(xbefore), y])); selected_precision)
      xafter = xpt[2]
      ptsafter = _real_roots(projy(Oscar.evaluate(f, [Rxy(xafter), y])); selected_precision)

      # Just guessing some precision, not optimal...
      yinterval = [ypt[1]-1//selected_precision, ypt[2]+1//selected_precision]
      result = ptsbefore
      diff = length(ptsbefore)-length(ptsafter)
      if diff != 0
         type = :ytangent
         if abs(diff) != 2
            type = :broken
         end
         # @assert abs(diff) == 2 "Please choose proper generic coordinate change $diff"
         if diff < 0
            result = ptsafter
         end
      end
      
      # Just verify that everything fits with float
      singy = Float64(result[findfirst(r->yinterval[1]<=r && r<=yinterval[2], result)])
      singindices = Int[]
      outr = typeof(result)()
      singfound = -1
      for i in 1:length(result)
         pt = result[i]
         if yinterval[1]<=pt && pt<=yinterval[2]
            if singfound == -1
               push!(outr, pt)
               singfound = i
            end
            # @assert Float64(pt) == singy "Please increase precision $singy $(Float64(pt))"
            if Float64(pt) != singy
               type = :broken
            end
            push!(singindices, i)
         else
            push!(outr, pt)
         end
      end

      if diff != 0
         ls = length(singindices)
         if ls != 2
            type = :broken
         end
         # @assert length(singindices) == 2 "Please choose proper generic coordinate change $ls"
      end

      return (outr, singfound, type)
   end
end


function msolve_sings(I ; info_level::Int=0, precision::Int=128, interval::Bool=true)
   # Just to remember how to get more debug info
   # (sings, rs2) = real_solutions(ideal(f) + ideal(derivative(f, y)); info_level=2, precision=selected_precision, interval=true)
   (sings, _) = Oscar.real_solutions(I; info_level, precision, interval)
   sort!(sings, lt=((a,b)->isless(a[1],b[1])))
   return sings
end


function isotopy_graph_from_msolve(IG::_IsotopyGraph, f_in, random_transform, selected_precision::Int)
   Rxy = parent(f_in)
   @assert nvars(Rxy)==2 "Need curve in affine plane"

   # This matrix is a generic rotation by definition
   f = f_in((random_transform*gens(Rxy))...)

   (x,y) = gens(Rxy)
   K = base_ring(Rxy)
   Ry, t = polynomial_ring(K, [:y])
   projy = hom(Rxy, Ry, [0, t[1]])
   sings = msolve_sings(ideal([f, derivative(f,y)]); precision=selected_precision, interval=true)
   println("There are ", length(sings), " points to investigate")
   # scale = 10/(max(xmax-xmin, ymax-ymin))
   svecs = []
   sindices = []
   stypes = Symbol[]

   # Compute points above and below the singularities.
   i = 0
   for (xpt, ypt) in sings
      i += 1
      # println(Float64(xpt[1])," ",Float64(xpt[2])," ",Float64(ypt[1])," ",Float64(ypt[2]))
      (svec, si, type) = analyse_singularity(f, projy, xpt, ypt; selected_precision)
      println("$i $(Float64(xpt[1]))")
      println("$i Was this singular? $type")
      if type == :broken
         return false
      end
      push!(svecs, svec)
      push!(sindices, si)
      # println("si: $si $(length(svec))")
      push!(stypes, type)
   end

   singcoords = [[sings[i][1][1], sings[i][2][1]] for i in 1:length(svecs)]
   _assemble_isotopy_graph!(IG, f, singcoords, stypes, svecs, sindices, random_transform, selected_precision, _real_roots)

   return true
end


################################################################################
################################################################################
##
## Main interface
## 
################################################################################
################################################################################

isotopy_graph_from_curve_inner(IG::_IsotopyGraph, f_in::QQMPolyRingElem, random_transform, selected_precision::Int) = isotopy_graph_from_msolve(IG, f_in, random_transform, selected_precision)::Bool

@doc raw"""
    draw_curve_tikz(f_in; filename::String="curve.tikz", selected_precision::Int=128, graph::Bool=false, custom_edge_plot=nothing)

Takes a polynomial in two variables and constructs a plot of the resulting real
algebraic curve.
"""
function draw_curve_tikz(f_in; filename::String="curve.tikz", selected_precision::Int=128, graph::Bool=false, custom_edge_plot=nothing)
   # Oscar.@check !isfile(filename) "Output file exists"
   IG = _IsotopyGraph()
   # Small rotation
   base_transform = matrix(QQ, [2052//2055 -111//2055; 111//2055 2052//2055])
   # Small shearing
   random_transform = matrix(QQ, base_transform)
   # base_transform *= matrix(QQ, [1 1//16; 0 1])
   result = isotopy_graph_from_curve_inner(IG, f_in, random_transform, selected_precision)
   println("result is $result")
   while !result
      # println("Turning!")
      random_transform *= base_transform
      result = isotopy_graph_from_curve_inner(IG, f_in, random_transform, selected_precision)
      println("result is $result")
   end

   scale = get_scale(IG, random_transform)
   io = open(filename, "w")
   if graph
      draw_graph_tikz(IG, io)
   else
      draw_curve_tikz(IG, scale, io; custom_edge_plot)
   end
   close(io)
end

export draw_curve_tikz
