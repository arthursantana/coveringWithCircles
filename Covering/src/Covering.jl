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

function distanceOutsideCircle(a, center, r)
    d = Voronoi.Geometry.distanceSquared(a, center)
    return d - r^2
end

function intersectedSegment(segment::Voronoi.Diagram.HalfEdge, center::Tuple{Real, Real}, r::Real)
    p = segment.origin
    q = segment.next.origin

    docp = distanceOutsideCircle(p, center, r)
    docq = distanceOutsideCircle(q, center, r)

    if docp <= 0 && docq <= 0
        el = head = Edge(p, nothing)
    else
        intersections = segmentArcIntersections(segment, center, r)

        if docq >= 0 && docp >= 0 && length(intersections) == 1
            el = head = Arc(p, nothing)
        else
            if docp < 0
                inside = true
                el = head = Edge(p, nothing)
            else
                inside = false
                el = head = Arc(p, nothing)
            end

            for i in intersections
                if inside
                    el.next = Arc(i, nothing);
                else
                    el.next = Edge(i, nothing);
                end

                el = el.next
                inside = !inside
            end
        end
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

function radianToDegrees(rad)
    return rad*180/π
end

function segmentAngle(p, q)
    num = q[2] - p[2]
    den = q[1] - p[1]
    if num == 0 && den > 0
        angle = 0
    elseif num == 0 && den < 0
        angle = π
    elseif den == 0 && num > 0
        angle = pi/2
    elseif den == 0 && num < 0
        angle = 3pi/2
    elseif num > 0 && den > 0
        angle = atan(num/den)
    elseif num > 0 && den < 0
        angle = π - atan(-num/den)
    elseif num < 0 && den < 0
        angle = π + atan(num/den)
    else # num < 0 && den > 0
        angle = -atan(-num/den)
    end

    return angle
end

function areaAndGradient(section::Section, r::Real)
    totalArea = 0
    el = section.borderHead
    gᵣ = 0
    gₛx = gₛy = 0.0

    firstTime = true # gimme do whiles, Julia, dammit!
    while (el != section.borderHead || firstTime)
        firstTime = false
        if (el isa Arc)
            θ1 = segmentAngle(section.center, el.origin)
            θ2 = segmentAngle(section.center, el.next.origin)
            θ = θ2 - θ1
            if θ < 0
                θ += 2π
            end

            totalArea += θ*(r^2)/2 # (θ/(2π))*π*(r^2)
            gᵣ += -r*θ # (θ/(2π))*2π*r
            gₛx += r*(sin(θ1) - sin(θ2))
            gₛy += r*(cos(θ2) - cos(θ1))
        else # isa Edge
            point = el
            shoelace = 0
            shoelace += el.origin[1]*el.next.origin[2] - el.origin[2]*el.next.origin[1]
            shoelace += el.next.origin[1]*section.center[2] - el.next.origin[2]*section.center[1]
            shoelace += section.center[1]*el.origin[2] - section.center[2]*el.origin[1]

            totalArea += shoelace/2
        end

        el = el.next
    end

    return totalArea, gᵣ, (gₛx, gₛy)
end

function areaAndGradient(P::Partition)
    n = length(P.sections)

    totalArea = gᵣ = 0
    gₛ = fill((0.0, 0.0), n)

    for i in 1:n
        a, g_r, g_s_i = areaAndGradient(P.sections[i], P.r)
        totalArea += a
        gᵣ += g_r
        gₛ[i] = g_s_i
    end

    return totalArea, gᵣ, gₛ
end

end # module
