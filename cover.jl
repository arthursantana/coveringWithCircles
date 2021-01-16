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
       W = SH.intersect(V, A, r)

       Draw.init(WIDTH, HEIGHT)
       Draw.coveringPartition(W, "xkcd:pink")
       Draw.coveringSection(A, 0, "xkcd:black")
       Draw.commit()
       print("")
   end

   function f(x)
       r, points = unpack(x)

       if drawing && drawingAll
           draw(r, points)
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

       V = Voronoi.Fortune.compute(points, WIDTH, HEIGHT)
       Voronoi.Intersect.intersect(V, Voronoi.Intersect.Rectangle(WIDTH, HEIGHT))
       W = SH.intersect(V, A, r)
       covered_area, gᵣ, gₛ = Covering.areaAndGradient(W)

       return WIDTH*HEIGHT - covered_area
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
    println("USO: cover.sh N_DE_PONTOS [DESENHAR]")
else
    if drawing
        import Draw
    end

    n = tryparse(Int64, ARGS[1])
    r = 0.1

    function randf(start, finish, n)
        v = rand(n)
        return map(x -> start + x*(finish-start), v)
    end

    #Random.seed!(1)
    points = convert(Array{Tuple{Real, Real}}, collect(zip(randf(1, WIDTH-1, n), randf(1, HEIGHT-1, n))))

    try
        coverWithCircles(n, WIDTH, HEIGHT, r, points)
    catch e
        @error "Deu pau:" exception=(e, catch_backtrace())
    end
end

