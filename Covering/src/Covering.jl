module Covering

import Voronoi

abstract type Element end

mutable struct Edge <: Element
   origin::Tuple{Real, Real}
   next::Union{Element, Nothing}
end

mutable struct Arc <: Element
   origin::Tuple{Real, Real}
   next::Union{Element, Nothing}
end

mutable struct Section
   center::Union{Tuple{Real, Real}, Nothing}
   borderHead::Union{Element, Nothing}
end

function segmentArcIntersections(segment::Voronoi.Diagram.HalfEdge, center::Tuple{Real, Real}, r::Real)
    function dot(a, b)
        return a[1]*b[1] + a[2]*b[2]
    end

    p = segment.origin
    q = segment.next.origin

    PQ = Voronoi.Geometry.subVector(q, p)
    CP = Voronoi.Geometry.subVector(p, center)

    a = dot(PQ, PQ)
    b = 2*dot(CP, PQ)
    c = dot(CP, CP) - r^2

    Δ = b^2 - 4a*c # discriminant

    if Δ < 0 # no intersection
        return []
    elseif Δ == 0 # one intersection
        α = (-b - sqrt(Δ))/2a

        if 0 <= α <= 1
            return [(p[1] + α*PQ[1], p[2] + α*PQ[2])]
        else
            return []
        end
    else # two intersections
        α1 = (-b - sqrt(Δ))/2a
        α2 = (-b + sqrt(Δ))/2a

        if (α1 < α2)
            if 0 <= α1 <= α2 <= 1
                return [(p[1] + α1*PQ[1], p[2] + α1*PQ[2]), (p[1] + α2*PQ[1], p[2] + α2*PQ[2])]
            elseif 0 <= α1 <= 1
                return [(p[1] + α1*PQ[1], p[2] + α1*PQ[2])]
            elseif 0 <= α2 <= 1
                return [(p[1] + α2*PQ[1], p[2] + α2*PQ[2])]
            else
                return []
            end
        else
            if 0 <= α2 <= α1 <= 1
                return [(p[1] + α2*PQ[1], p[2] + α2*PQ[2]), (p[1] + α1*PQ[1], p[2] + α1*PQ[2])]
            elseif 0 <= α1 <= 1
                return [(p[1] + α1*PQ[1], p[2] + α1*PQ[2])]
            elseif 0 <= α2 <= 1
                return [(p[1] + α2*PQ[1], p[2] + α2*PQ[2])]
            else
                return []
            end
        end
    end
end

function intersectedSegment(segment::Voronoi.Diagram.HalfEdge, center::Tuple{Real, Real}, r::Real)
    p = segment.origin
    if Voronoi.Geometry.distanceSquared(p, center) < r^2
        inside = true
        el = head = Edge(segment.origin, nothing)
    else
        inside = false
        el = head = Arc(segment.origin, nothing)
    end

    for i in segmentArcIntersections(segment, center, r)
        if inside
            el.next = Arc(i, nothing);
        else
            el.next = Edge(i, nothing);
        end

        el = el.next
        inside = !inside
    end

    return head, el
end

function voronoiRegionToCoveringSection(region::Voronoi.Diagram.Region, r::Real)
    he = region.borderHead.next
    head, tail = intersectedSegment(he, region.generator, r)

    el = head
    while (he != region.borderHead)
        he = he.next
        tail.next, tail = intersectedSegment(he, region.generator, r)
    end
    tail.next = head

    return Section(region.generator, head)
end

mutable struct Partition
    sections::Array{Covering.Section}
    r::Real
end

function voronoiDiagramToPartition(V::Voronoi.Diagram.DCEL, r::Real)
    n = length(V.regions)
    sections = Array{Covering.Section, 1}(undef, n)

    for i in 1:n
        sections[i] = Covering.voronoiRegionToCoveringSection(V.regions[i], r)
    end

    return Partition(sections, r)
end

end # module
