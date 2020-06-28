drawing = (haskey(ENV, "DRAW") && ENV["DRAW"] != "" && ENV["DRAW"] != "0")

using Random
using Printf

include("Algencan/AlgencanWrapper.jl")

import Voronoi
import Covering

function coverWithCircles(n, WIDTH, HEIGHT, r, points)
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
       if !drawing
           return
       end

       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
       P = Covering.voronoiDiagramToPartition(V, r)
       area, gᵣ, gₛ = Covering.areaAndGradient(P)

       Draw.init(WIDTH, HEIGHT)
       Draw.voronoiDiagram(V)
       Draw.coveringPartition(P)
       Draw.commit()
       print("")
   end

   function f(x)
       r, points = unpack(x)

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

       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
       P = Covering.voronoiDiagramToPartition(V, r)
       covered_area, gᵣ, gₛ = Covering.areaAndGradient(P)

       return WIDTH*HEIGHT - covered_area
   end

   function ∇c(ind, x)
       r, points = unpack(x)

       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
       P = Covering.voronoiDiagramToPartition(V, r)
       area, gᵣ, gₛ = Covering.areaAndGradient(P)

       return collect(1:2n+1), pack(gᵣ, gₛ)
   end

   l = zeros(2n + 1)
   l[1] = 1
   u = zeros(2n + 1)
   u[1] = 500.0
   for i in 2:2n+1
       u[i] = 100.0
   end

   x, fx = AlgencanWrapper.optimize(n = 2n + 1, m = 1,

                                    f = f, g = ∇f, h = ∇∇f,
                                    equatn = [1],
                                    c = c, jac = ∇c,

                                    x = pack(r, points),
                                    l = l,
                                    u = u,

                                    nvparam = 1,
                                    vparam = [
                                              "ITERATIONS-OUTPUT-DETAIL 1",
                                             ],
                                    epsopt = 1.0e-4,
                                    epsfeas = 1.0e-4
                                    )

   r, points = unpack(x)
   draw(r, points)
   println("\n\n=================\nRESULTADO\n=================\n")
   r /= 100
   println("r = ", r)
   println("s = ", 1/r)
   println("pontos: ")
   for p in points
       println("(", p[1], ", ", p[2], ")")
   end
   if drawing
       println("Digite ENTER para sair.")
       readline(stdin)
   end
end








if length(ARGS) == 0
    println("USO: cover.sh N_DE_PONTOS [DESENHAR]")
else
    if drawing
        import Draw
    end

    WIDTH = 100.0
    HEIGHT = 100.0

    n = tryparse(Int64, ARGS[1])
    r = 10

    function randf(start, finish, n)
        v = rand(n)
        return map(x -> start + x*(finish-start), v)
    end

    #Random.seed!(1)
    points = convert(Array{Tuple{Real, Real}}, collect(zip(randf(1, WIDTH-1, n), randf(1, HEIGHT-1, n))))

    coverWithCircles(n, WIDTH, HEIGHT, r, points)
end

