using Random
using Printf

include("Algencan/AlgencanWrapper.jl")
include("sh.jl")

import Voronoi
import Covering

drawing = (haskey(ENV, "DRAW") && ENV["DRAW"] != "" && ENV["DRAW"] != "0")
drawingAll = (haskey(ENV, "DRAWALL") && ENV["DRAWALL"] != "" && ENV["DRAWALL"] != "0")

TRIES = 100
WIDTH = 1.0
HEIGHT = 1.0
INSTANCE = 1

As = [
      Array{Covering.Section, 1}([
                                  Covering.Section(Array{Tuple{Real, Real}, 1}([
                                                                                (0.1, 0.1),
                                                                                (0.9, 0.1),
                                                                                (0.9, 0.9),
                                                                                (0.1, 0.9),
                                                                               ])),
                                 ]),
      Array{Covering.Section, 1}([
                                  Covering.Section(Array{Tuple{Real, Real}, 1}([
                                                                                (0.5, 0.5),
                                                                                (0.4, 0.4),
                                                                                (0.5, 0.1),
                                                                                (0.6, 0.4),
                                                                               ])),
                                  Covering.Section(Array{Tuple{Real, Real}, 1}([
                                                                                (0.5, 0.5),
                                                                                (0.6, 0.4),
                                                                                (0.9, 0.5),
                                                                                (0.6, 0.6),
                                                                               ])),
                                  Covering.Section(Array{Tuple{Real, Real}, 1}([
                                                                                (0.5, 0.5),
                                                                                (0.6, 0.6),
                                                                                (0.5, 0.9),
                                                                                (0.4, 0.6),
                                                                               ])),
                                  Covering.Section(Array{Tuple{Real, Real}, 1}([
                                                                                (0.5, 0.5),
                                                                                (0.4, 0.6),
                                                                                (0.1, 0.5),
                                                                                (0.4, 0.4),
                                                                               ])),
                                 ]),
      Array{Covering.Section, 1}([
                                  Covering.Section(Array{Tuple{Real, Real}, 1}([
                                                                                (0, 0.125),
                                                                                (0.75, 0.5),
                                                                                (1 ,0.875),
                                                                                (0.25, 0.5),
                                                                               ])),
                                  Covering.Section(Array{Tuple{Real, Real}, 1}([
                                                                                (0, 1),
                                                                                (0, 0.5),
                                                                                (0.5, 1),
                                                                               ]))
                                 ]),
      Array{Covering.Section, 1}([
                                  Covering.Section(Array{Tuple{Real, Real}, 1}([
                                                                                (0.1, 0.1),
                                                                                (0.2, 0.2),
                                                                                (0.2, 0.8),
                                                                                (0.1, 0.9),
                                                                               ])),
                                  Covering.Section(Array{Tuple{Real, Real}, 1}([
                                                                                (0.1, 0.1),
                                                                                (0.9, 0.1),
                                                                                (0.8, 0.2),
                                                                                (0.2, 0.2),
                                                                               ])),
                                  Covering.Section(Array{Tuple{Real, Real}, 1}([
                                                                                (0.9, 0.1),
                                                                                (0.9, 0.9),
                                                                                (0.8, 0.8),
                                                                                (0.8, 0.2),
                                                                               ])),
                                  Covering.Section(Array{Tuple{Real, Real}, 1}([
                                                                                (0.9, 0.9),
                                                                                (0.1, 0.9),
                                                                                (0.2, 0.8),
                                                                                (0.8, 0.8),
                                                                               ])),
                                 ]),
     ]

if haskey(ENV, "INSTANCE")
    INSTANCE = parse(Int, ENV["INSTANCE"])
end
A = As[INSTANCE]

function calculateArea()
    total_area = 0
    for Aⱼ in A
        t, _, _ = Covering.areaAndGradient(Aⱼ, 0)
        total_area += t
    end
    return total_area
end
total_area = calculateArea()

function draw(WIDTH, HEIGHT, r, points, n)
   if !drawing
       return
   end

   if length(points) == 0
       return
   end

   V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
   Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))

   Draw.init(WIDTH, HEIGHT)

   for Aⱼ in A
       W = SH.intersect(V, Aⱼ, r)
       #Draw.coveringPartition(W, "xkcd:pink")
       Draw.coveringPartition(W, nothing)
       Draw.filledPolygon(Aⱼ, "blue")
   end

   Draw.circles(r, points)
   Draw.commit()
   Draw.savefig("bateria/" * string(INSTANCE) * "/" * string(n))
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

   function draw(r, points) if !drawing
           return
       end

       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))

       Draw.init(WIDTH, HEIGHT)

       for Aⱼ in A
           W = SH.intersect(V, Aⱼ, r)
           Draw.coveringPartition(W, "xkcd:pink")
           Draw.filledPolygon(Aⱼ, "xkcd:black")
       end

       Draw.circles(r, points)
       Draw.commit()
       print("")
   end

   function f(x)
       r, points = unpack(x)

       if drawing && drawingAll
           draw(r, points)
       end

       #println("R: ", r)
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

       uncovered_area = total_area
       for Aⱼ in A
           W = SH.intersect(V, Aⱼ, r)
           area, _, _ = Covering.areaAndGradient(W)
           uncovered_area -= area
       end

       #println("C: ", uncovered_area)
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
       #println("GRAD C.: ", pack(gᵣ, gₛ))
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

   #draw(r, points)
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
       #draw(r, points)

       println("\n\n=================\nRESULTADO\n=================\n")
       println("r = ", r)
       println("s = ", 1/r)
       println("pontos: ")
       for p in points
           println("(", p[1], ", ", p[2], ")")
       end
       #if drawing
       #    println("Digite ENTER para sair.")
       #    readline(stdin)
       #end

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

    #total_area, _, _ = Covering.areaAndGradient(A, 0)
    ϵ = 10e-4
    
    Random.seed!(1)

    for k in 1:TRIES
        #println("Tentativa ", k)
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

        r = (0.5 + (rand(1)[1]))/n

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

           uncovered_area = total_area
           for Aⱼ in A
               W = SH.intersect(V, Aⱼ, r)
               area, _, _ = Covering.areaAndGradient(W)
               uncovered_area -= area
           end

           if uncovered_area < ϵ
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
