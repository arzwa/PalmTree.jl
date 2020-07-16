module PalmTree

using Parameters
using Luxor
using Printf
import Luxor: RGB
import NewickTree: isleaf, id, distance, children
export TreeLayout, drawtree, nodemap

"""
    TreeLayout(node; dims=(400,300), cladogram=false)

Define a tree layout for a tree with given root node. Methods that should
be implemented in order to get be able to get a TreeLayout are `distance` (unless
`cladogram = true`), `children`, `isleaf` and `id`.
"""
struct TreeLayout{I,T}
    edges ::Vector{Tuple{I,I}}
    coord ::Dict{I,Point}
    nodes ::Dict{I,T}
    leaves::Vector{I}
    dims  ::Tuple
end

Base.getindex(tl::TreeLayout, i) = tl.coord[i]

function TreeLayout(n::T; dims=(400,300), cladogram=false) where T
    if !hasmethod(distance, Tuple{T}) && !cladogram
        @warn "No distance method for time $T, setting `cladogram = true`"
        cladogram = true
    end
    I = typeof(id(n))
    coord = Dict{I,Point}()
    nodes = Dict{I,T}()
    edges = Tuple{I,I}[]
    leaves = I[]
    yleaf = -1.
    function walk(n, x)
        d = cladogram ? 1. : distance(n)
        d = isnan(d) ? 0. : d
        nodes[id(n)] = n
        if isleaf(n)
            yleaf += 1.
            currx = x + d
            coord[id(n)] = Point(currx, yleaf)
            push!(leaves, id(n))
            return yleaf
        else
            childy = []
            currx = x + d
            for c in children(n)
                push!(childy, walk(c, currx))
                push!(edges, (id(n), id(c)))
            end
            y = sum(childy)/length(childy)
            coord[id(n)] = Point(x + d, y)
            return y
        end
    end
    walk(n, 0.)
    tl = TreeLayout(edges, coord, nodes, leaves, dims)
    scale!(tl, dims...)
    cladogram && cladogram!(tl)
    tl
end

function scale!(tl::TreeLayout, width, height)
    @unpack Δx, Δy = bbox(tl)
    xscale = width /Δx
    yscale = height/Δy
    for (k, p) in tl.coord
        tl.coord[k] *= (xscale, yscale)
    end
end

function cladogram!(tl::TreeLayout)
    @unpack leaves = tl
    @unpack p2 = bbox(tl)
    for i in leaves
        tl.coord[i] = Point(p2.x, tl.coord[i].y)
    end
end

function bbox(tl::TreeLayout)
    maxx, minx, maxy, miny = 0., 0., 0., 0.
    for p in values(tl.coord)
        @unpack x, y = p
        maxx = x > maxx ? x : maxx
        maxy = y > maxy ? y : maxy
        minx = x < minx ? x : minx
        miny = y < miny ? y : miny
    end
    (p1=Point(minx, miny), p2=Point(maxx, maxy),
     Δx=maxx-minx, Δy=maxy-miny)
end

function drawtree(tl::TreeLayout; color=(x)->RGB(), bblend=true)
    @unpack coord, edges, nodes = tl
    for (a,b) in edges
        corner = Point(coord[a].x, coord[b].y)
        poly([coord[a], corner, coord[b]])
        bblend ?
            setblend(blend(coord[a], corner, color(nodes[a]),color(nodes[b]))) :
            sethue(color(nodes[b]))
        Luxor.strokepath()
    end
end

function colorbar(pos, len, cs; δ=0.1)
    init = deepcopy(pos)
    Δy = δ*len
    for v=0.0:δ:1.0-δ
        ca = get(cs, 1 - v)
        cb = get(cs, 1 - v + δ)
        bl = Luxor.blend(pos, pos + (0, Δy), cb, ca)
        Luxor.line(pos, pos + (0, Δy))
        Luxor.setblend(bl)
        Luxor.strokepath()
        pos += (0, Δy)
    end
    return (p1=init, p2=pos)
end

function nodemap(tl::TreeLayout, nodes, f::Function)
    for n in nodes
        f(n, tl[id(n)])
    end
end

nodemap(tl::TreeLayout, f::Function) = nodemap(tl, values(tl.nodes), f)

end
