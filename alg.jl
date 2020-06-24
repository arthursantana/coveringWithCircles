import Voronoi
import Covering
import Draw

using Random
using Printf

include("AlgencanWrapper.jl")

function randf(start, finish, n)
   v = rand(n)
   return map(x -> start + x*(finish-start), v)
end

function coverWithCircles()
   WIDTH = 100.0
   HEIGHT = 100.0

   n = tryparse(Int64, ARGS[1])
   r = 10

   #Random.seed!(1)
   points = convert(Array{Tuple{Real, Real}}, collect(zip(randf(1, WIDTH-1, n), randf(1, HEIGHT-1, n))))

   function pack(r, points)
       x = zeros(2n + 1)
       x[1] = r
       for (i, p) in enumerate(points)
           x[2i] = p[1]
           x[2i + 1] = p[2]
       end

       return x
   end

   function unpack(x)
       r = x[1]

       points = Array{Tuple{Real, Real}}([])

       for i in 1:n
           append!(points, [(x[2i], x[2i + 1])])
       end

       return r, points
   end

   function draw(r, points)
       println("# r: ", r)
       println("# points: ", points)
       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
       P = Covering.voronoiDiagramToPartition(V, r)
       area, gᵣ, gₛ = Covering.areaAndGradient(P)
       println("área coberta: ", area)

       Draw.init(WIDTH, HEIGHT)
       Draw.voronoiDiagram(V)
       Draw.coveringPartition(P)
       Draw.commit()
   end

   function f(x)
       r, points = unpack(x)

       println("\nEvaluating function")
       draw(r, points)
       return r
   end

   function ∇f(x)
       grad = zeros(2n + 1)
       grad[1] = 1
       return grad
   end

   function ∇∇f(x)
       return [], [], []
   end

   function c(ind, x)
       r, points = unpack(x)

       println("\nEvaluating constraints\nr = ", r)
       println("points = ", points)
       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
       P = Covering.voronoiDiagramToPartition(V, r)
       covered_area, gᵣ, gₛ = Covering.areaAndGradient(P)

       println("c = ", WIDTH*HEIGHT - covered_area)
       return WIDTH*HEIGHT - covered_area
   end

   function ∇c(ind, x)
       r, points = unpack(x)

       println("\nEvaluating constraint gradient\nr = ", r)
       println("points = ", points)
       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
       P = Covering.voronoiDiagramToPartition(V, r)
       area, gᵣ, gₛ = Covering.areaAndGradient(P)

       println("∂c/∂r: ", gᵣ)
       println("∂c/∂points: ", gₛ)

       return collect(0:n-1), pack(gᵣ, gₛ)
   end

   draw(r, points)
   readline(stdin)

   l = zeros(2n + 1)
   l[1] = 1e-10
   u = zeros(2n + 1)
   u[1] = 500.0
   for i in 2:2n+1
       u[i] = 100.0
   end

   x, fx = AlgencanWrapper.optimize(n = 2n + 1, m = 1,

                                    f = f, g = ∇f, h = ∇∇f,
                                    equatn = [1],
                                    c = c,# jac = ∇c,

                                    x = pack(r, points),
                                    l = l,
                                    u = u,

                                    nvparam = 1,
                                    vparam = [
                                              "ITERATIONS-OUTPUT-DETAIL 1",
                                             ],
                                    )

   r, points = unpack(x)
   draw(r, points)
   println("\n\n=================\nRESULTADO")
   println("=================\n")
   r /= 100
   println("r = ", r)
   println("s = ", 1/r)
   readline(stdin)
end

coverWithCircles()
