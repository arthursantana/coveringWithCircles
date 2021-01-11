module SH

import Covering
import Voronoi

function intersectWithHalfSpace(W::Covering.Section, p::Tuple{Real, Real}, q::Tuple{Real, Real})
    pq = Voronoi.Geometry.subVector(q, p)
    v = Voronoi.Geometry.rotateVectorCCW(pq)

    function isVisible(point)
        function dot(a, b)
            return a[1]*b[1] + a[2]*b[2]
        end

        u = Voronoi.Geometry.subVector(point, p)

        return (dot(u, v) >= 0)
    end

    function intersectionPoint(a::Tuple{Real, Real}, b::Tuple{Real, Real})
        p1 = p
        dir1 = pq

        p2 = a
        dir2 = Voronoi.Geometry.subVector(b, a)

        vec = [p2[1] - p1[1], p2[2] - p1[2]]
        M = [dir1[1] dir2[1]; dir1[2] dir2[2]]
        res = M\vec
        return (p1[1] + res[1]*dir1[1], p1[2] + res[1]*dir1[2])
    end

    newW = Covering.Section((0, 0), nothing)
    global current = nothing
    function putPoint(p::Tuple{Real, Real})
        global current
        if current == nothing
            newW.borderHead = Covering.Edge(p, nothing)
            current = newW.borderHead
        else
            current.next = Covering.Edge(p, nothing)
            current = current.next
        end
    end

    w = W.borderHead
    ranOnce = false

    while w != W.borderHead || !ranOnce
        if w == nothing
            break
        end
        a = w.origin
        b = w.next.origin

        if isVisible(a)
            if isVisible(b)
                putPoint(b)
            else
                putPoint(intersectionPoint(a, b))
            end
        elseif isVisible(b)
            putPoint(intersectionPoint(a, b))
            putPoint(b)
        else
            # print none
        end

        w = w.next
        ranOnce = true
    end

    if current != nothing
        current.next = newW.borderHead
    end

    return newW
end

function intersectOne(Vᵢ::Voronoi.Diagram.Region, A::Covering.Section)
    W = Covering.Section(A)

    W = intersectWithHalfSpace(W, Vᵢ.borderHead.origin, Vᵢ.borderHead.next.origin)

    v = Vᵢ.borderHead.next
    while v != nothing && v != Vᵢ.borderHead
        W = intersectWithHalfSpace(W, v.origin, v.next.origin)
        v = v.next
    end

    return W
end

function intersectWithCircle(W::Covering.Section, center::Tuple{Real, Real}, r::Real)
    function isVisible(point)
        dist = Voronoi.Geometry.distanceSquared(center, point)

        return (dist <= r^2)
    end

    newW = Covering.Section(center, nothing)
    global current = nothing
    function putPoint(p::Tuple{Real, Real}, isArc=false)
        global current
        if current == nothing
            if isArc
                newW.borderHead = Covering.Arc(p, nothing)
            else
                newW.borderHead = Covering.Edge(p, nothing)
            end
            current = newW.borderHead
        else
            if isArc
                current.next = Covering.Arc(p, nothing)
            else
                current.next = Covering.Edge(p, nothing)
            end
            current = current.next
        end
    end

    w = W.borderHead
    ranOnce = false

    while w != W.borderHead || !ranOnce
        if w == nothing
            break
        end
        a = w.origin
        b = w.next.origin

        if isVisible(a)
            if isVisible(b)
                putPoint(b)
            else
                putPoint(Covering.segmentArcIntersections(a, b, center, r)[1], true) # put arc
            end
        elseif isVisible(b)
            putPoint(Covering.segmentArcIntersections(a, b, center, r)[1])
            putPoint(b)
        else
            intersections = Covering.segmentArcIntersections(a, b, center, r)

            # arcs with no intersection or just tangential are ignored
            if length(intersections) == 2
                putPoint(intersections[1])
                putPoint(intersections[2], true) # put arc
            end
        end

        w = w.next
        ranOnce = true
    end

    if current != nothing
        current.next = newW.borderHead
    end

    return newW
end

function intersect(V::Voronoi.Diagram.DCEL, A::Covering.Section, r::Real)
    n = length(V.regions)
    sections = Array{Covering.Section, 1}(undef, n)

    for i in 1:n
        sections[i] = intersectWithCircle(intersectOne(V.regions[i], A), V.regions[i].generator, r)
    end

    return Covering.Partition(sections, r)
end

end
