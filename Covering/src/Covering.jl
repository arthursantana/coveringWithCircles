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
   center::Tuple{Real, Real}
   borderHead::Element
end

function regionToCoveringSection(region::Voronoi.Diagram.Region, r::Real)
    he = region.borderHead
    el = head = Edge(he.origin, nothing)
    he = he.next

    while (he != region.borderHead)
        el.next = Edge(he.origin, nothing)
        el = el.next
        he = he.next
    end

    el.next = Edge(he.origin, head)

    return Section(region.generator, head)
end

end # module
