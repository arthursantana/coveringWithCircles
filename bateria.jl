using Random
using Printf

include("Algencan/AlgencanWrapper.jl")
include("sh.jl")

import Voronoi
import Covering

drawing = (haskey(ENV, "DRAW") && ENV["DRAW"] != "" && ENV["DRAW"] != "0")
drawingAll = (haskey(ENV, "DRAWALL") && ENV["DRAWALL"] != "" && ENV["DRAWALL"] != "0")

WIDTH = 1.0
HEIGHT = 1.0

A = Covering.Section(Array{Tuple{Real, Real}, 1}([
                                 (0, 0.125),
                                 (0.75, 0.5),
                                 (1 ,0.875),
                                 (0.25, 0.5),
                                ]))
total_area, _, _ = Covering.areaAndGradient(A, 0)

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
#function draw(r, points)
#   if !drawing
#       return
#   end
#
#   V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
#   Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
#   W = SH.intersect(V, A, r)
#
#   Draw.init(WIDTH, HEIGHT)
#   Draw.coveringPartition(W, "xkcd:pink")
#   Draw.coveringSection(A, 0, "xkcd:black")
#   #Draw.commit()
#   Draw.savefig("bateria/" * string(n) * ".png")
#   println("...")
#end
#
#function draw(WIDTH, HEIGHT, r, points, n)
#    V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
#    Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
#    P = Covering.voronoiDiagramToPartition(V, r)
#    area, gᵣ, gₛ = Covering.areaAndGradient(P)
#
#    Draw.init(WIDTH, HEIGHT)
#    Draw.voronoiDiagram(V)
#    Draw.coveringPartition(P)
#    Draw.savefig("img/" * string(n) * ".png")
#end
function draw(WIDTH, HEIGHT, r, points, n)
   if !drawing
       return
   end

   if length(points) == 0
       return
   end

   V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
   Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
   W = SH.intersect(V, A, r)

   Draw.init(WIDTH, HEIGHT)
   Draw.coveringPartition(W, nothing)
   Draw.coveringSection(A, 0, "xkcd:black")
   Draw.circles(r, points)
   Draw.commit()
   Draw.savefig("bateria/" * string(n) * ".png")
   println("...")
end


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

   function f(x)
       r, points = unpack(x)

       return r
   end

   function ∇f(x)
       grad = zeros(2n + 1)
       grad[1] = 1
       #println("GRAD: ", grad)
       return grad
   end

   function ∇∇f(x)
       return [], [], []
   end

   function c(ind, x)
       r, points = unpack(x)

       #points, repeats = avoidRepeats(points)

       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
       W = SH.intersect(V, A, r)
       covered_area, gᵣ, gₛ = Covering.areaAndGradient(W)

       return total_area - covered_area
   end

   function ∇c(ind, x)
       r, points = unpack(x)

       #points, repeats = avoidRepeats(points)

       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
       W = SH.intersect(V, A, r)
       covered_area, gᵣ, gₛ = Covering.areaAndGradient(W)

       #for i in repeats
       #    push!(gₛ, (0, 0))
       #end
       #println("GRAD C.: ", pack(gᵣ, gₛ))

       return collect(1:2n+1), pack(gᵣ, gₛ)
   end

   l = zeros(2n + 1)
   l[1] = 1e-6

   u = zeros(2n + 1)
   u[1] = 5.0
   for i in 2:2n+1
       u[i] = 1.0
   end

   #draw(r, points)
   #println("Aperte Enter pra continuar")
   #readline(stdin)

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
       #draw(r, points)

       println("\n\n=================\nRESULTADO\n=================\n")
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

       return r, points
   catch exc
       println("Exception: ", exc)
       return 1e10, []
   end
end






function randf(start, finish, n)
    v = rand(n)
    return map(x -> start + x*(finish-start), v)
end

if drawing
    import Draw
end


function main()
    n = tryparse(Int64, ARGS[1])
    WIDTH = 1.0
    HEIGHT = 1.0

    rmin = 1e10
    pointsmin = []

    total_area, _, _ = Covering.areaAndGradient(A, 0)
    ϵ = 10e-4

    for k in 1:1000
        Random.seed!(1)
        points = convert(Array{Tuple{Real, Real}}, collect(zip(randf(1, WIDTH-1, n), randf(1, HEIGHT-1, n))))

        function insideHalfSpace(point, p, q)
            pq = Voronoi.Geometry.subVector(q, p)
            v = Voronoi.Geometry.rotateVectorCCW(pq)

            function dot(a, b)
                return a[1]*b[1] + a[2]*b[2]
            end

            u = Voronoi.Geometry.subVector(point, p)

            return (dot(u, v) >= 0)
        end

        function pointIsInsideConvexPolygon(point, polygon)
            p = polygon.borderHead

            first = true
            while p != polygon.borderHead || first
                first = false
                q = p.next
                if !insideHalfSpace(point, p.origin, q.origin)
                    return false
                end

                p = p.next
            end

            return true
        end

        for (i, p) in enumerate(points)
            while !pointIsInsideConvexPolygon(points[i], A)
                println("Ponto $i ainda está fora. Re-sorteando.")
                points[i] = (randf(1, WIDTH-1, 1)[1], randf(1, HEIGHT-1, 1)[1])
            end
        end

        r = randf(0.05, 0.15, 1)[1]
        r = 0.1

        try
            r, points = coverWithCircles(n, WIDTH, HEIGHT, r, points)
            #@info "points: $(points)"
        catch exc
            println("Exception: ", exc)
            return 1e10, []
        end

        #draw(WIDTH, HEIGHT, r, points, n)

        @info "n = $(n); k = $(k); rmin = $(rmin); r = $(r)"
        if r < rmin
           V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
           Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
           W = SH.intersect(V, A, r)
           covered_area, gᵣ, gₛ = Covering.areaAndGradient(W)

           if total_area - covered_area < ϵ
               rmin = r
               pointsmin = points
           else
               @warn "UNFEASIBLE"
           end
        end
    end

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
