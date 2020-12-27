drawing = (haskey(ENV, "DRAW") && ENV["DRAW"] != "" && ENV["DRAW"] != "0")
drawingAll = false

using Random
using Printf

include("Algencan/AlgencanWrapper.jl")

import Voronoi
import Covering

if drawing
    import Draw

    function draw(WIDTH, HEIGHT, r, points, n)
        V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
        Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
        P = Covering.voronoiDiagramToPartition(V, r)
        area, gᵣ, gₛ = Covering.areaAndGradient(P)

        Draw.init(WIDTH, HEIGHT, false)
        Draw.voronoiDiagram(V)
        Draw.coveringPartition(P)
        Draw.savefig("bateria/" * string(n) * ".png")
    end
end

function randf(start, finish, n)
    v = rand(n)
    return map(x -> start + x*(finish-start), v)
end

#function avoidRepeats(points)
#    println("LÁ VAMOS NÓS")
#    println("points antes: ", points)
#    repeat = []
#    for i in 1:length(points)
#        for j in 1:length(points)
#            if j == i
#                continue
#            end
#
#            if points[i][1] == points[j][1] && points[i][2] == points[j][2]
#                push!(repeat, points[j])
#                deleteat!(points, j)
#            end
#        end
#    end
#
#    println("points depois: ", points)
#    println("repeats depois: ", repeat)
#    return points, repeat
#end
#
function coverWithCircles(n, WIDTH, HEIGHT, r, points)
   W = WIDTH
   H = HEIGHT

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

   function f(x)
       r, points = unpack(x)

       if drawing && drawingAll
           draw(W, H, r, points)
       end

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

       #points, repeats = avoidRepeats(points)

       #println("CALCULANDO C DE ", r, ", ", points)

       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
       P = Covering.voronoiDiagramToPartition(V, r)
       covered_area, gᵣ, gₛ = Covering.areaAndGradient(P)

       return WIDTH*HEIGHT - covered_area
   end

   function ∇c(ind, x)
       r, points = unpack(x)

       #points, repeats = avoidRepeats(points)

       #println("CALCULANDO C DE ", r, ", ", points)

       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
       P = Covering.voronoiDiagramToPartition(V, r)
       area, gᵣ, gₛ = Covering.areaAndGradient(P)

       #for i in repeats
       #    push!(gₛ, (0, 0))
       #end

       return collect(1:2n+1), pack(gᵣ, gₛ)
   end

   l = zeros(2n + 1)
   l[1] = 1e-6

   u = zeros(2n + 1)
   u[1] = 5.0
   for i in 2:2n+1
       u[i] = 1.0
   end

   try
       x, fx = AlgencanWrapper.optimize(n = 2n + 1, m = 1,

                                        f = f, g = ∇f, h = ∇∇f,
                                        equatn = [1],
                                        c = c, jac = ∇c,

                                        x = pack(r, points),
                                        l = l,
                                        u = u,

                                        nvparam = 1,
                                        vparam = [
                                                  "ITERATIONS-OUTPUT-DETAIL 10",
                                                 ],
                                        #checkder = 1,
                                        epsopt = 1.0e-8,
                                        epsfeas = 1.0e-8
                                       )
       r, points = unpack(x)

       return r, points
   catch exc
       println("Exception: ", exc)
       return 1e10, []
   end
end









function main()
    n = tryparse(Int64, ARGS[1])
    WIDTH = 1.0
    HEIGHT = 1.0

    rmin = 1e10
    pointsmin = []

    for k in 1:100
        points = convert(Array{Tuple{Real, Real}}, collect(zip(randf(1, WIDTH-1, n), randf(1, HEIGHT-1, n))))
        r = randf(0.5, 1.5, 1)[1]

        r, points = coverWithCircles(n, WIDTH, HEIGHT, r, points)

        if r < rmin
            rmin = r
            pointsmin = points
        end
    end

    sleep(5)
    println("\n\n=================\nRESULTADO (n = " * string(n) * ")\n=================\n")
    println("r = ", rmin)
    println("s = ", 1/rmin)
    println("pontos: ")
    for p in pointsmin
        println("(", p[1], ", ", p[2], ")")
    end
    println("")

    if drawing
        draw(WIDTH, HEIGHT, rmin, pointsmin, n)
        #readline(stdin)
    end
end

main()
