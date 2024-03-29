module Draw

import Voronoi
import Covering

TIKZ = true
MATPLOTLIB = false

if MATPLOTLIB
    using PyPlot
    using PyCall

    const patch = PyNULL()
    function __init__()
        copy!(patch, pyimport("matplotlib.patches"))
    end

    plt = PyPlot
end

tikz_text = ""
ax = nothing
WIDTH = nothing
HEIGHT = nothing
FRAME = 0.1

if TIKZ
    COLOR_BLACK = "black"
    COLOR_SILVER = "silver"
    COLOR_AZURE = "azure"
    COLOR_ORANGERED = "orangered"
    COLOR_TOMATO = "blue"
    COLOR_MAGENTA = "magenta"
    COLOR_GOLD = "gold"
    COLOR_GRAY = "gray"
    COLOR_TAN = "tan"
    COLOR_RED = "red"
    COLOR_GREEN = "green"
else
    COLOR_BLACK = "xkcd:black"
    COLOR_SILVER = "xkcd:silver"
    COLOR_AZURE = "xkcd:azure"
    COLOR_ORANGERED = "xkcd:orangered"
    COLOR_TOMATO = "xkcd:tomato"
    COLOR_MAGENTA = "xkcd:magenta"
    COLOR_GOLD = "xkcd:gold"
    COLOR_GRAY = "xkcd:gray"
    COLOR_TAN = "xkcd:tan"
    COLOR_RED = "xkcd:red"
    COLOR_GREEN = "xkcd:green"
end

colors = ["white"]
#colors = ["xkcd:grey", "xkcd:crimson", "xkcd:gold", "xkcd:green", "xkcd:azure", "xkcd:beige", "xkcd:silver", "xkcd:lavender", "xkcd:lightgreen", "xkcd:magenta", "xkcd:ivory", "xkcd:maroon", "xkcd:orange", "xkcd:orangered", "xkcd:orchid", "xkcd:pink", "xkcd:plum", "xkcd:gold", "xkcd:salmon", "xkcd:sienna", "xkcd:lime", "xkcd:tan", "xkcd:tomato", "xkcd:violet", "xkcd:wheat", "xkcd:indigo", "xkcd:yellowgreen", "xkcd:chocolate", "xkcd:coral", "xkcd:brown"]


function clear(title::String)
   global plt, ax
   global WIDTH, HEIGHT, FRAME

   plt.cla()
   plt.title(title)
   if WIDTH > 0
      ax.set_aspect("equal")
      ax.set_xlim([0 - FRAME, WIDTH + FRAME])
      ax.set_ylim([0 - FRAME, HEIGHT + FRAME])
      ax.grid(false)
      ax.get_xaxis().set_visible(false)
      ax.get_yaxis().set_visible(false)
   else
      ax.grid(true)
      ax.get_xaxis().set_visible(true)
      ax.get_yaxis().set_visible(true)
   end
end

function init(w, h, ion=true)
    if TIKZ
        global tikz_text

        tikz_text = "\\begin{tikzpicture}[thick, scale=3]\n"
    end

    if MATPLOTLIB
        global plt, ax
        global WIDTH, HEIGHT

        plt.pygui(true)
        if ion
            plt.ion()
        end

        plt.clf()
        #plt.xkcd() # uncomment for generalized wobbliness
        ax = plt.gca() # get current axes

        WIDTH = w
        HEIGHT = h
        clear("")
    end
end


function legend(xlabel, ylabel)
   global plt

   plt.xlabel(xlabel)
   plt.ylabel(ylabel)
   plt.legend()
end

function commit()
    if TIKZ
        global tikz_text

        tikz_text *= "\\end{tikzpicture}\n\n"
    end

    if MATPLOTLIB
        global plt

        plt.draw()
    end
end

function savefig(filename)
    if TIKZ
        global tikz_text

        filename *= ".tex"
        println("Escrevendo Tikz no arquivo ", filename)

        open(filename, "w") do io
            write(io, tikz_text)
        end;
    end

    if MATPLOTLIB
        global plt

        filename *= ".png"

        plt.savefig(filename)
    end
end

function arc(p::Tuple{Real, Real}, r, θ1, θ2, color)
    if TIKZ
        global tikz_text

        t1 = 2π*θ1/360
        x = p[1] + r*cos(t1)
        y = p[2] + r*sin(t1)

        if θ2 < θ1
            θ2 += 360
        end

        tikz_text *= "\\draw[$color] ($(x), $(y)) arc ($(θ1):$(θ2):$(r)cm);\n"
    end

    if MATPLOTLIB
        global ax

        ax.add_artist(patch.Arc((p[1], p[2]), 2r, 2r, 0, θ1, θ2, color=color, linewidth=3, zorder=1))
    end
end

function circle(p::Tuple{Real, Real}, color, fill, r, l, zorder=1)
    if TIKZ
        global tikz_text

        tikz_text *= "\\draw[$color] ($(p[1]), $(p[2])) circle ($(r)cm);\n"
        if fill
            tikz_text *= "\\fill[$color] ($(p[1]), $(p[2])) circle ($(r)cm);\n"
        end
    end

    if MATPLOTLIB
        global ax

        ax.add_artist(patch.Circle((p[1], p[2]), color=color, radius=r, fill=fill, zorder=zorder, linewidth=l))
    end
end

function point(p::Tuple{Real, Real}, color, fill)
    circle(p, color, fill, 0.005, 1)
end

function circle(p::Tuple{Real, Real}, color, r)
   circle(p, color, false, r, 1)
end

function startFill(color)
    if TIKZ
        global tikz_text

        tikz_text *= "\\filldraw[color=$color, fill=$color" * "!10] "
    end
end

function fillVertex(p::Tuple{Real, Real})
    if TIKZ
        global tikz_text

        tikz_text *= "($(p[1]), $(p[2])) -- "
    end
end

function endFill()
    if TIKZ
        global tikz_text

        tikz_text *= "cycle;\n"
    end
end

function line(p1::Tuple{Real, Real}, p2::Tuple{Real, Real}, color)
    if TIKZ
        global tikz_text

        tikz_text *= "\\draw[$color] ($(p1[1]), $(p1[2])) -- ($(p2[1]), $(p2[2]));\n"
    end

    if MATPLOTLIB
        global plt

        plt.plot([p1[1], p2[1]], [p1[2], p2[2]], color=color, linestyle="-", linewidth=3, zorder=1)
    end
end

function thinLine(p1::Tuple{Real, Real}, p2::Tuple{Real, Real}, color)
   global plt

   plt.plot([p1[1], p2[1]], [p1[2], p2[2]], color=color, linestyle="-", linewidth=1, zorder=1)
end

function plot(f, color, start, finish)
   global plt

   x = range(start, stop=finish, length=1000)
   y = map(f, x)

   plt.plot(x, y, color=color, linewidth=3, zorder=2)
end

function plot(x, y, color, label, style)
   global plt

   plt.plot(x, y, color=color, label=label, linestyle=style, linewidth=3, zorder=2)
end



function fortuneIteration(V::Voronoi.Diagram.DCEL, T::Voronoi.BeachLine.BST, Q::Voronoi.EventQueue.Heap, points::Array{Tuple{Real, Real}}, ly::Real)
	#clear("Computing Voronoi Voronoi.Diagram using Fortune's Algorithm")
	clear("")

   # draw points
   for p in points
      point(p, COLOR_BLACK, (ly <= p[2]))
      if ly < p[2]
         f = Voronoi.Geometry.parabola(p, ly)
         plot(f, COLOR_SILVER, 0, WIDTH)
      end
   end

   # calculate parabolas and breakpoints
	beachLineFoci = Voronoi.BeachLine.beachLine(T, ly)

   # draw beachline
   start = (0, HEIGHT)
   i = 2
   while i <= size(beachLineFoci)[1] + 2
      if i <= size(beachLineFoci)[1]
         finish = beachLineFoci[i]
      else
         finish = (WIDTH, start[2])
      end

      p = beachLineFoci[i-1]
		f = Voronoi.Geometry.parabola(p, ly)

		if f == nothing # point is over the sweep line
         if start == nothing # special case where there the first couple of points are on the same y coordinate
         else
            Draw.line((p[1], ly), (p[1], start[2]), COLOR_AZURE)
         end
		else
         if 0 <= start[1]
            st = start[1]
         else
            st = 0
         end

         if finish[1] <= WIDTH
            fn = finish[1]
         else
            fn = WIDTH
         end

			plot(f, COLOR_AZURE, st, fn)
		end

      start = finish
      i += 2
	end

   # draw circumcircles
   i = 1
   while i < Q.pos
      if isa(Q.data[i], Voronoi.EventQueue.CircleEvent) && !(Q.data[i].removed)
         b = Q.data[i].disappearingArc
         a = b.prev
         c = b.next
         O = Q.data[i].center
         r = O[2] - Q.data[i].coordinates[2]
         circle(O, COLOR_ORANGERED, r)
         thinLine(O, b.region.generator, COLOR_ORANGERED)
         point(O, COLOR_MAGENTA, true)
         point((O[1], O[2] - r), COLOR_ORANGERED, true)
      end

      i += 1
   end

   # draw diagram edges
   for he in V.halfEdges
      if he.origin != nothing && he.twin.origin != nothing
         Draw.line(he.origin, he.twin.origin, COLOR_BLACK)
      end
   end

   # draw sweepline
	line((0, ly), (WIDTH, ly), COLOR_GOLD)
	commit()
end

function voronoiDiagram(V::Voronoi.Diagram.DCEL)
    global WIDTH, HEIGHT

	clear("")

    regions = Voronoi.Diagram.regionBorders(V)

    # draw regions
    i = 1
    for region in regions
        if size(region[1])[1] > 0
            ax.fill(region[1], region[2], colors[(i % size(colors)[1])+ 1])
        end
        i += 1
    end

    # draw diagram edges
    for he in V.halfEdges
        if he.origin != nothing
            if he.twin != nothing
                Draw.line(he.origin, he.twin.origin, COLOR_SILVER)
            elseif he.next != nothing
                Draw.line(he.origin, he.next.origin, COLOR_SILVER)
            end
        end
    end

    # draw generators
    for r in V.regions
        point(r.generator, COLOR_BLACK, true)
    end
end

function crossingPoint(center, dir, r)
    vector = (dir[1] - center[1], dir[2] - center[2])
    n = sqrt(vector[1]^2 + vector[2]^2) 
    unitVector = (vector[1]/n, vector[2]/n)
    crossingPoint = (center[1] + r * unitVector[1], center[2] + r * unitVector[2])
    
    return crossingPoint
end

function circles(r::Real, points)
    for center in points
        circle(center, COLOR_GRAY, false, r, 1)
        point(center, COLOR_BLACK, true)
    end
end

function filledPolygon(section::Covering.Section, color=COLOR_TAN)
    if TIKZ
        if section.borderHead == nothing
            return
        end

        startFill(color)

        el = section.borderHead.next
        fillVertex(el.origin)
        while (el != section.borderHead)
            el = el.next
            fillVertex(el.origin)
        end

        endFill()
    else
        return coverginSection(section, 0, color)
    end
end

function coveringSection(section::Covering.Section, r::Real, color=COLOR_TAN)
    if section.borderHead == nothing
        return
    end

    center = section.center

    #circle(center, COLOR_GRAY, false, r, 1)

    el = section.borderHead.next
    i = 0
    if el isa Covering.Arc
        arc(section.center, r,
            Covering.radianToDegrees(Covering.segmentAngle(center, el.origin)),
            Covering.radianToDegrees(Covering.segmentAngle(center, el.next.origin)),
            COLOR_TOMATO)
        #cp = crossingPoint(center, el.origin, r)
        #line(cp, center, COLOR_GREEN)
        #line(cp, el.origin, "xkcd:pale green")
    else
        if color !== nothing
            line(el.origin, el.next.origin, color)
        end
        #line(el.origin, center, COLOR_GREEN)
    end

    while (el != section.borderHead)
        el = el.next
        i += 1
        if el isa Covering.Arc
            arc(section.center, r,
                Covering.radianToDegrees(Covering.segmentAngle(center, el.origin)),
                Covering.radianToDegrees(Covering.segmentAngle(center, el.next.origin)),
                COLOR_TOMATO)

            #cp = crossingPoint(center, el.origin, r)
            #line(cp, center, COLOR_GREEN)
            #line(cp, el.origin, COLOR_GREENreen")
        else
            if color !== nothing
                line(el.origin, el.next.origin, color)
            end
            #line(el.origin, center, COLOR_GREEN)
        end
        for p in Covering.Covering.segmentArcIntersections(Voronoi.Diagram.HalfEdge(el.origin, false, nothing, 
                                                                          Voronoi.Diagram.HalfEdge(el.next.origin, false, nothing, nothing, nothing, center),
                                                                          nothing, center), center, r)
            #point(p, COLOR_RED, true)
        end
    end
end

function coveringPartition(P::Covering.Partition, color=COLOR_TAN)
    for section in P.sections
        Draw.coveringSection(section, P.r, color)
    end
end

end # module
