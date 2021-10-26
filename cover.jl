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

A = Array{Covering.Section, 1}([
        Covering.Section(Array{Tuple{Real, Real}, 1}([
            (0.2, 0.5),
            (0.8, 0.5),
            (0.8, 0.75),
            (0.2, 0.75),
        ])),
        Covering.Section(Array{Tuple{Real, Real}, 1}([
            (0.4, 0.5),
            (0.4, 0.25),
            (0.6, 0.25),
            (0.6, 0.5),
        ]))
    ])

#A = Array{Covering.Section, 1}([
#        Covering.Section(Array{Tuple{Real, Real}, 1}([
#            (0, 0.125),
#            (0.75, 0.5),
#            (1 ,0.875),
#            (0.25, 0.5),
#        ])),
#        Covering.Section(Array{Tuple{Real, Real}, 1}([
#            (0, 1),
#            (0, 0.5),
#            (0.5, 1),
#        ]))
#    ])

function f()
    total_area = 0
    for Aⱼ in A
        t, _, _ = Covering.areaAndGradient(Aⱼ, 0)
        total_area += t
    end
    return total_area
end
total_area = f()

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

       Draw.init(WIDTH, HEIGHT)
       Draw.circles(r, points)

       for Aⱼ in A
           W = SH.intersect(V, Aⱼ, r)
           Draw.coveringPartition(W, "xkcd:pink")
           Draw.coveringSection(Aⱼ, 0, "xkcd:black")
       end

       Draw.commit()
       print("")
   end

   function f(x)
       r, points = unpack(x)

       if drawing && drawingAll
           draw(r, points)
       end

       println("R: ", r)
       return r
   end

   function ∇f(x)
       grad = zeros(2n + 1)
       grad[1] = 1
       println("GRAD: ", grad)
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

       uncovered_area = total_area
       for Aⱼ in A
           W = SH.intersect(V, Aⱼ, r)
           area, _, _ = Covering.areaAndGradient(W)
           uncovered_area -= area
       end

       println("C: ", uncovered_area)
       return uncovered_area
   end

   function ∇c(ind, x)
       r, points = unpack(x)

       #points, repeats = avoidRepeats(points)

       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))

       gᵣ = 0
       gₛ = fill((0.0, 0.0), n)
       for Aⱼ in A
           W = SH.intersect(V, Aⱼ, r)
           _, gᵣⱼ, gₛⱼ = Covering.areaAndGradient(W)
           gᵣ += gᵣⱼ
           for i in 1:n
               gₛ[i] = (gₛ[i][1] + gₛⱼ[i][1], gₛ[i][2] + gₛⱼ[i][2])
           end
       end

       #for i in repeats
       #    push!(gₛ, (0, 0))
       #end
       println("GRAD C.: ", pack(gᵣ, gₛ))
       #println("Aperte Enter pra continuar")
       #readline(stdin)

       return collect(1:2n+1), pack(gᵣ, gₛ)
   end

   l = zeros(2n + 1)
   l[1] = 1e-14

   u = zeros(2n + 1)
   u[1] = 10e20
   for i in 2:2n+1
       u[i] = 1.0
   end

   draw(r, points)
   #println("Aperte Enter pra continuar")
   #readline(stdin)

   try
       x, fx = AlgencanWrapper.optimize(n = 2n + 1, m = 1,

                                        f = f, g = ∇f,# h = ∇∇f,
                                        equatn = [1],
                                        c = c, jac = ∇c,

                                        x = pack(r, points),
                                        l = l,
                                        u = u,

                                        nvparam = 2,
                                        vparam = [
                                                  "ITERATIONS-OUTPUT-DETAIL 10",
                                                  "PENALTY-PARAMETER-INITIAL-VALUE 1000"
                                                 ],
                                        #checkder = 1,
                                        epsopt = 1.0e-6,
                                        epsfeas = 1.0e-6
                                        )

       r, points = unpack(x)
       draw(r, points)

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
   catch exc
       println("Exception: ", exc)
       return 1e10, []
   end
end








if length(ARGS) == 0
    println("USO: cover.sh N_DE_PONTOS [DESENHAR] [DESENHAR_TODA_AVALIAÇÃO]")
else
    if drawing
        import Draw
    end

    n = tryparse(Int64, ARGS[1])

    function randf(start, finish, n)
        v = rand(n)
        return map(x -> start + x*(finish-start), v)
    end

    Random.seed!(1)
    points = convert(Array{Tuple{Real, Real}}, collect(zip(randf(1, WIDTH-1, n), randf(1, HEIGHT-1, n))))

    r = (0.5 + (rand(1)[1]))/n

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

    function pointIsInsideSomeConvexPolygon(point, polygons)
        for polygon in polygons
            if pointIsInsideConvexPolygon(point, polygon)
                return true
            end
        end

        return false
    end

    for (i, p) in enumerate(points)
        while !pointIsInsideSomeConvexPolygon(points[i], A)
            #println("Ponto $i ainda está fora. Re-sorteando.")
            points[i] = (randf(1, WIDTH-1, 1)[1], randf(1, HEIGHT-1, 1)[1])
        end
    end

    try
        coverWithCircles(n, WIDTH, HEIGHT, r, points)
    catch e
        @error "Deu pau:" exception=(e, catch_backtrace())
    end
end

