import Voronoi
import Covering
import Draw

using Random
using Printf

#using Algencan

function unpack(x, n)
   r = x[1]

   points = Array{Tuple{Real, Real}}([])

   for i in 1:n
       append!(points, [(x[2i], x[2i + 1])])
   end

   return r, points
end


#function randf(start, finish, n)
#   v = rand(n)
#   return map(x -> start + x*(finish-start), v)
#end

#Em todos os testes, vamos considerar o retangulo de vertices (0,0); (3,0); (3,1); (0,1).

confs = Array{Real}[]

#########
# TESTE 1
#########

push!(confs, [0.3,
1.11, 0.4,
1.44, 0.25,
1.4, 0.77])

#########
# TESTE 2
#########

push!(confs, [0.3,
1.11, 0.4,
1.44, 0.25,
1.4, 0.77,
0.01, 0.01,
2.0, 0.8,
0.2, 0.73,
2.65, 0.51001])

#########
# TESTE 3
#########

push!(confs, [0.3,
1.11, 0.698,
0.887, 0.554,
1.4, 0.77,
1.554, 0.1,
1.8, 0.7,
0.99, 0.177,
1.805, 0.44])

#########
# TESTE 4
#########

push!(confs, [0.2,
1.11, 0.698,
1.05, 0.45,
1.322, 0.66,
1.39, 0.209,
2.99, 0.99,
1.12, 0.22,
1.505, 0.44,
0.01, 0.99,
0.01, 0.01,
2.99, 0.01])

#########
# TESTE 5
#########

push!(confs, [0.2777,
0.01, 0.499,
2.99, 0.499,
2.99, 0.99,
1.5, 0.99,
1.5, 0.01,
0.01, 0.99,
0.01, 0.01,
2.99, 0.01])

#########
# TESTE 6
#########

push!(confs, [0.3,
0.28, 0.8,
0.28, 0.6,
0.28, 0.4,
0.28, 0.2,
0.8, 0.8,
0.8, 0.6,
0.8, 0.4,
0.8, 0.2,
1.32, 0.8,
1.32, 0.6,
1.32, 0.4,
1.32, 0.2,
1.84, 0.8,
1.84, 0.6,
1.84, 0.4,
1.84, 0.2,
2.36, 0.8,
2.36, 0.6,
2.36, 0.4,
2.36, 0.2,
2.84, 0.8,
2.84, 0.6,
2.84, 0.4,
2.84, 0.2])

#########
# TESTE 7
#########

push!(confs, [0.5,
0.0, 0.5,
3.0, 0.5])

#########
# TESTE 8
#########

push!(confs, [0.5,
0.5, 0.5,
2.5, 0.5,
1.5, 0.5])

#########
# TESTE 9
#########

push!(confs, [0.5,
0.0, 1.0,
0.0, 0.0,
3.0, 0.0,
3.0, 1.0])

#########
# TESTE 10
#########

push!(confs, [0.5,
0.0, 1.0,
0.0, 0.0,
3.0, 0.0,
3.0, 1.0,
1.0, 0.5,
2.0, 0.5])

function gradientDescent()
   WIDTH = 3
   HEIGHT = 1

   command = nothing
   for (i, x) in enumerate(confs)
      n = Int((length(x) - 1)/2)
      r, points = unpack(x, n)

      println("=============================")
      println("CASO ", i)
      println("r: ", r)
      println("pontos:")
      for p in points
          println("(", p[1], ", ", p[2], ")")
      end
      println()

      iterations = 1
      fortunes = 1

      time = @elapsed V, f, ∇f, ξ = Voronoi.Optimization.init(points, WIDTH, HEIGHT)


      μ = 0.5
      α = 0.01
      ϵ = 10^-4
      λᵣ = 0.5
      λₛ = 1

      ρ = 0.1

      times = time

      V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
      Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
      P = Covering.voronoiDiagramToPartition(V, r)
      #@printf("f = %.10g\t\t|∇f| = %.10g\t\ttime: %s seconds\n", f, ξ, time)
      area, gᵣ, gₛ = Covering.areaAndGradient(P)
      #println("área coberta: ", area)
      @printf("area coberta= %.8g\ng = %.10g\n∂g/∂r = %.8g\n", area, (WIDTH*HEIGHT - area), gᵣ)
      for i in 1:length(gₛ)
          x, y = gₛ[i]
          @printf("∂g/∂c_%d = (%.8g, %.8g)\n", i, x, y)
      end
      @printf("\n")

      Draw.init(WIDTH, HEIGHT)
      Draw.voronoiDiagram(V)
      Draw.coveringPartition(P)
      Draw.commit()

      times += 0#time
      fortunes += 0#forts
      iterations += 1

      if command != "a"
          #println("Press Return for a step or enter \"a\" to animate until the end.")
          command = readline(stdin)
      else
          sleep(0.00001) # make sure stuff is drawn
      end
      println("\n")

      #println("n: ", n, "\t\titer: ", iterations, "\t\tfortunes: ", fortunes, "\t\ttime: ", times)
  end
end

gradientDescent()
