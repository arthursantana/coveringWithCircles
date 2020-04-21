import Voronoi
import Covering
import Draw
import Optimization

using Random
using Printf


function randf(start, finish, n)
   v = rand(n)
   return map(x -> start + x*(finish-start), v)
end

function gradientDescent()
   WIDTH = 300.0
   HEIGHT = 100.0

   #seed = 100
   #println("Random seed: ", seed)
   #Random.seed!(seed)

   command = nothing
   for n in 7:100
      if sqrt(n) == floor(sqrt(n))
         continue
      end
      #points = convert(Array{Tuple{Real, Real}}, collect(zip(randf(1, WIDTH-1, n), randf(1, HEIGHT-1, n))))
      #points = convert(Array{Tuple{Real, Real}}, [(111, 40), (144, 25), (140, 77)])
      points = convert(Array{Tuple{Real, Real}}, [(111, 40), (144, 25), (140, 77), (1, 1), (200, 80), (20, 73), (265, 51.001)])

      iterations = 1
      fortunes = 1

      time = @elapsed V, f, ∇f, ξ = Optimization.init(points, WIDTH, HEIGHT)


      μ = 0.5
      α = 0.01
      ϵ = 10^-4
      λᵣ = 0.5
      λₛ = 1

      ρ = 0.1

      r = WIDTH/(sqrt(n))
      r = 30

      times = time

      while ξ > ϵ
          V = Voronoi.Fortune.compute(points)
          Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
          P = Covering.voronoiDiagramToPartition(V, r)
          @printf("f = %.10g\t\t|∇f| = %.10g\t\ttime: %s seconds\n", f, ξ, time)
          area, gᵣ, gₛ = Covering.areaAndGradient(P)
          @printf("r = %.8g\t\tarea coberta= %.8g\ng = %.10g\n∂g/∂r = %.8g\n", r, area, (WIDTH*HEIGHT - area)/10000, gᵣ/100)
          for i in 1:length(gₛ)
              x, y = gₛ[i]
              @printf("∂g/∂c_%d = (%.8g, %.8g)\n", i, x/100, y/100)
          end
          @printf("\n")

          Draw.init(WIDTH, HEIGHT)
          Draw.voronoiDiagram(V)
          Draw.coveringPartition(P)
          Draw.commit()

          f = r + ρ*(WIDTH*HEIGHT - area)

          ∇obj = fill((0.0, 0.0), n+1)
          for i in 1:n
              ∇obj[i] = (ρ*gₛ[i][1], ρ*gₛ[i][2])
          end
          ∇obj[n+1] = (1 - ρ*gᵣ, 0)
          ξ = Optimization.norm2(∇obj)

          for i in 1:n
              a = points[i][1] - λₛ*∇obj[i][1]
              b = points[i][2] - λₛ*∇obj[i][2]
              points[i] = (a, b)
          end
          r -= λᵣ*∇obj[n+1][1]

          ρ *= 1.005
          ρ += 0.0001

          #time = @elapsed V, f, ∇f, ξ, forts = Optimization.lineSearch(points, d, α, λ₀, μ, ϵ, V, f, ∇f, ξ)

          times += 0#time
          fortunes += 0#forts
          iterations += 1

          if command != "a"
              println("Press Return for a step or enter \"a\" to animate until the end.")
              command = readline(stdin)
          else
              sleep(0.00001) # make sure stuff is drawn
          end
      end

      println("n: ", n, "\t\titer: ", iterations, "\t\tfortunes: ", fortunes, "\t\ttime: ", times)
   end
end

gradientDescent()
