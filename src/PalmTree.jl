module PalmTree

using Parameters
using Luxor
using Printf
import Luxor: RGB
import NewickTree: isleaf, id, distance, children
export TreeLayout, drawtree, nodemap

"""
    TreeLayout(node)

Define a tree layout for a tree with given root node. Methods that should
be implemented in order to get be able to get a TreeLayout are `distance`,
`children`, `isleaf` and `id`.

```julia
using NewickTree, PalmTree, Luxor

tr = readnw("((A:77.3,B:77.3):16.4,(C:46.9,(D:23.4,E:23.4):23.4):46.9);")
tl = TreeLayout(tr[1])
@svg begin
    origin(Point(10,10))
    drawtree(tl)
    setfont("Noto mono", 16)
    nodemap(tl, (k, p)->
        settext(" "*tr.leaves[k], p, valign="center"), keys(tr.leaves))
end
```

```julia
r = backtrack(wm, ccd)
function leafname(n)
   γ, _ = parse.(Int16, split(n, "."))
   return 0 < γ <= length(ccd.leaves)  ? " "*ccd.leaves[γ] : ""
end
tl = TreeLayout(r)
cladogram!(tl)

@svg begin
   origin(Point(10,10))
   drawtree(tl, color=(n)->startswith(n, "0") ? RGB(0.8,0.8,0.8) : RGB())
   nodemap(tl, (k, p)->settext(leafname(k), p, valign="center"), tl.leaves)
end 800 320 "/tmp/tree.svg"
"""
struct TreeLayout{I}
    edges ::Vector{Tuple{I,I}}
    coord ::Dict{I,Point}
    leaves::Vector{I}
end

Base.getindex(tl::TreeLayout, i) = tl.coord[i]

function TreeLayout(n; dim=(400,300), cladogram=false)
    I = typeof(id(n))
    coord = Dict{I,Point}()
    edges = Tuple{I,I}[]
    leaves = I[]
    yleaf = -1.
    function walk(n, x)
        if isleaf(n)
            yleaf += 1.
            currx = x + distance(n)
            coord[id(n)] = Point(currx, yleaf)
            push!(leaves, id(n))
            return yleaf
        else
            childy = []
            currx = x + distance(n)
            for c in children(n)
                push!(childy, walk(c, currx))
                push!(edges, (id(n), id(c)))
            end
            y = sum(childy)/length(childy)
            coord[id(n)] = Point(x + distance(n), y)
            return y
        end
    end
    walk(n, 0.)
    tl = TreeLayout(edges, coord, leaves)
    scale!(tl, dim...)
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

function drawtree(tl::TreeLayout; color=(x)->RGB())
    @unpack coord, edges = tl
    for (a,b) in edges
        corner = Point(coord[a].x, coord[b].y)
        poly([coord[a], corner, coord[b]])
        setblend(blend(coord[a], corner, color(a), color(b)))
        Luxor.strokepath()
    end
end

function colorbar(pos, len; δ=0.1)
    init = deepcopy(pos)
    Δy = δ*len
    for v=0.0:δ:1.0-δ
        ca = get(ColorSchemes.viridis, 1 - v)
        cb = get(ColorSchemes.viridis, 1 - v + δ)
        bl = Luxor.blend(pos, pos + (0, Δy), cb, ca)
        Luxor.line(pos, pos + (0, Δy))
        Luxor.setblend(bl)
        Luxor.strokepath()
        pos += (0, Δy)
    end
    return (p1=init, p2=pos)
end

function nodemap(tl::TreeLayout, f, ks=keys(tl.coord))
    for k in ks
        f(k, tl[k])
    end
end

end
